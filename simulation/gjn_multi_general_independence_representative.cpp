// Simulation to record the joint and marginal occurance happend to test the asymptotic independence
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <ctime>
#include <chrono>
#include <string>
#include <sstream>
#include <random>
#include <iomanip>
#include <numeric>
#include <cctype>
using namespace std;
using namespace chrono;

const double inf = numeric_limits<double>::infinity();
// Declare or define all needed global variables
double T = 0;
double last_T = 0;
int K;
int obs = 3000;
vector<double> gengammaE;
vector<double> gengammaS;
vector<double> rho;
vector<double> nominalArrivalRate;
vector<vector<double>> transitionMatrix;
vector<double> externalArrivalRate;
vector<double> serviceRate;
vector<double> serviceTimes;
vector<int> threshold;
vector<vector<double>> timeRecord;
vector<vector<int>> occurRecord;
vector<vector<int>> occurRecord_joint;
vector<int> N, M;
vector<vector<vector<double>>> result;
default_random_engine generator(random_device{}());

// Utility: inverse CDF of geometric distribution
int geom_ppf(double p, double prob, int loc) {
    return static_cast<int>(ceil(log(1.0 - p) / log(1.0 - prob))) + loc;
}

void initializeOccurRecords() {
    occurRecord.resize(K);
    for (int k = 0; k < K; ++k) {
        occurRecord[k].resize(threshold[k] + 1, 0);
    }

    int rows = threshold[0] + 1;
    int cols = threshold[1] + 1;
    occurRecord_joint.resize(rows, vector<int>(cols, 0));
}

class SERVER;  // Forward declare only SERVER (optional, see below)

class CUSTOMER {
private:
    int _category;
    double _arrivalTime, _commitedTime, _departureTime;
    int _serverCategory;
    SERVER* _serverPointer;

public:
    CUSTOMER(int category)
        : _category(category),
          _arrivalTime(-1),
          _commitedTime(-1),
          _departureTime(-1),
          _serverCategory(-1),
          _serverPointer(nullptr) {}

    int getType() const { return _category; }
    double getArrival() const { return _arrivalTime; }
    double setArrival(double t) { return _arrivalTime = t; }
    double getDeparture() const { return _departureTime; }
    double setDeparture(double t) { return _departureTime = t; }
    double getCommited() const { return _commitedTime; }
    double setCommited(double t) { return _commitedTime = t; }

    int getServer() const { return _serverCategory; }
    void setServer(int cat) { _serverCategory = cat; }

    SERVER* getServerPointer() const { return _serverPointer; }
    void setServerPointer(SERVER* ptr) { _serverPointer = ptr; }
};

class FIFOQUEUE {
private:
    vector<shared_ptr<CUSTOMER>> _queue;
    double _servedTime;
    int _size;
    SERVER* _server;

public:
    FIFOQUEUE(double meanServiceTime)
        : _servedTime(meanServiceTime), _size(0), _server(nullptr) {}

    // string to_string() const {
    //     stringstream ss;
    //     ss << "FIFOQUEUE:size " << _queue.size();
    //     return ss.str();
    // }

    bool isEmpty() const { return _size == 0; }
    //const vector<shared_ptr<CUSTOMER>>& getItems() const { return _queue; }
    int getSize() const { return _size; }

    void joinQueue(const shared_ptr<CUSTOMER>& customer) {
        _queue.push_back(customer);
        _size += 1;
    }

    shared_ptr<CUSTOMER> commitQueue() {
        if (_queue.empty()) throw runtime_error("Attempted to pop from an empty queue");
        auto customer = _queue.front();
        _queue.erase(_queue.begin());
        _size -= 1;
        return customer;
    }

    SERVER* getServer() const { return _server; }
    void setServer(SERVER* server) { _server = server; }
};

class SERVER {
private:
    int _category, _status, _customerCategory;
    double _currentRate;
    CUSTOMER* _currentCustomer;

    vector<double> _serviceRates, _meanServiceTime;
    vector<shared_ptr<FIFOQUEUE>> _queues;

public:
    SERVER(int category, int status = 0)
        : _category(category), _status(status),
          _customerCategory(-1), _currentRate(-1),
          _currentCustomer(nullptr) {}

    int getType() const { return _category; }
    int getStatus() const { return _status; }

    const vector<double>& getRates() const { return _serviceRates; }
    const vector<double>& getMeans() const { return _meanServiceTime; }
    const vector<shared_ptr<FIFOQUEUE>>& getQueues() const { return _queues; }

    void setQueue(shared_ptr<FIFOQUEUE> queue, int queueType) {
        _queues[queueType] = queue;
    }

    vector<double> setRates(const vector<double>& serviceRates) {
        _serviceRates = serviceRates;
        int l = _serviceRates.size();
        _queues.resize(l);
        _meanServiceTime.resize(l);
        for (int i = 0; i < l; ++i) {
            if (_serviceRates[i] != 0) {
                _queues[i] = make_shared<FIFOQUEUE>(1.0 / _serviceRates[i]);
                _meanServiceTime[i] = 1.0 / _serviceRates[i];
            } else {
                _queues[i] = make_shared<FIFOQUEUE>(0.0);
                _meanServiceTime[i] = inf;
            }
        }
        return _serviceRates;
    }

    int getCustomer() const { return _customerCategory; }

    CUSTOMER* commitNow(CUSTOMER* customer, double currentTime) {
        _status = 1;
        _customerCategory = customer->getType();
        _currentRate = _serviceRates[0];
        _currentCustomer = customer;
        customer->setCommited(currentTime);
        customer->setServer(_category);
        customer->setServerPointer(this);
        return customer;
    }

    CUSTOMER* releaseNow(double currentTime) {
        _status = 0;
        _customerCategory = -1;
        _currentRate = -1;
        CUSTOMER* customer = _currentCustomer;
        _currentCustomer = nullptr;
        customer->setDeparture(currentTime);
        customer->setServerPointer(nullptr);
        return customer;
    }
};

vector<shared_ptr<FIFOQUEUE>> Queues;
vector<shared_ptr<SERVER>> Servers;
vector<pair<shared_ptr<CUSTOMER>, double>> q;

// gamma sampling
double gamma_sample(double shape, double scale) {
    gamma_distribution<double> distribution(shape, scale);
    return distribution(generator);
}

// Arrival distribution
double arrival(double interArrivalRate, int type) {
    if (interArrivalRate == 0.0) {
        return inf;
    }
    return gamma_sample(gengammaE[type], 1.0 / (interArrivalRate * gengammaE[type]));
}

// Internal transition type
int internalTransitionType(int k) {
    vector<double>& probs = transitionMatrix[k];
    vector<double> cumArrival(probs.size());
    partial_sum(probs.begin(), probs.end(), cumArrival.begin());
    double random_number = uniform_real_distribution<double>(0.0, 1.0)(generator);

    for (int i = 0; i < K; ++i) {
        if (random_number <= cumArrival[i]) {
            return i;
        }
    }

    if (random_number > cumArrival[K - 1]) {
        return K;  // exit
    }

    return -1;  // default fallback
}

// Service distribution
double service(double interServiceTime, int type) {
    return gamma_sample(gengammaS[type], interServiceTime / gengammaS[type]);
}

// Next event
int nextEvent() {
    int index = -1;
    double minval = inf;
    for (int i = 0; i < 3 * K; ++i) {
        if (q[i].second < minval && q[i].second >= 0.0) {
            minval = q[i].second;
            index = i;
        }
    }
    return index;
}

// Event complete
void eventComplete(int index) {
    last_T = T;
    T = q[index].second;
    q[index].first = nullptr;  // Assuming placeholder for no CUSTOMER
    q[index].second = -1;
}

// Clock tick (external arrival)
void clockTick(int k) {
    if (q[k].second == -1) {
        double arrivalTime = T + arrival(externalArrivalRate[k], k);
        q[k].second = arrivalTime;
        shared_ptr<CUSTOMER> cust = make_shared<CUSTOMER>(k);
        cust->setArrival(arrivalTime);
        q[k].first = cust;  // <-- fix: also set the pointer
    }
}

void timeCheckInd(double timeHorizon) {
    if (T > timeHorizon / 4.0) {
        double timecut = (3.0 * timeHorizon) / (4.0 * obs);

        if (floor(T / timecut) != floor(last_T / timecut) && T > timeHorizon / 4.0) {
            for (int k = 0; k < K; ++k) {
                int diff = N[k] - M[k];
                int idx = std::min(diff, threshold[k]);
                occurRecord[k][idx] += 1;
            }

            int diff0 = N[0] - M[0];
            int diff1 = N[1] - M[1];
            int idx0 = std::min(diff0, threshold[0]);
            int idx1 = std::min(diff1, threshold[1]);

            if (diff0 >= threshold[0] && diff1 >= threshold[1]) {
                occurRecord_joint[threshold[0]][threshold[1]] += 1;
            } else if (diff0 >= threshold[0]) {
                occurRecord_joint[threshold[0]][idx1] += 1;
            } else if (diff1 >= threshold[1]) {
                occurRecord_joint[idx0][threshold[1]] += 1;
            } else {
                occurRecord_joint[idx0][idx1] += 1;
            }
        }
    }
}

// Route external customer
shared_ptr<SERVER> customerRouteExternal(const CUSTOMER& customer) {
    int customerType = customer.getType();
    int serverIndex = customerType;
    return Servers[serverIndex];
}

// Route internal transition
shared_ptr<SERVER> customerRouteTransit(int k) {
    int customerType = internalTransitionType(k);
    if (customerType == K) return nullptr;
    return Servers[customerType];
}

// Updated simulate function with null pointer checks

tuple<double, vector<int>, vector<int>, vector<vector<double>>, vector<vector<vector<double>>>>
simulate(double timeHorizon) {
    while (T < timeHorizon) {
        int nextEventIndex = nextEvent();

        if (nextEventIndex == -1) {
            vector<double> firstArrivalTime(K);
            for (int k = 0; k < K; ++k) {
                clockTick(k);
                firstArrivalTime[k] = q[k].second;
            }
            int firstEventIndex = distance(firstArrivalTime.begin(), min_element(firstArrivalTime.begin(), firstArrivalTime.end()));
            auto customer = make_shared<CUSTOMER>(firstEventIndex);
            q[firstEventIndex].first = customer;
            customer->setArrival(q[firstEventIndex].second);

            auto server = customerRouteExternal(*customer);
            auto queue = server->getQueues()[0];
            queue->joinQueue(customer);

            eventComplete(firstEventIndex);
            timeCheckInd(timeHorizon);
            N[firstEventIndex]++;

            queue->commitQueue();
            server->commitNow(customer.get(), T);
            double serviceTime = service(server->getMeans()[0], server->getType());
            q[firstEventIndex + K].first = customer;
            q[firstEventIndex + K].second = T + serviceTime;

            clockTick(firstEventIndex);
        }

        else if (nextEventIndex < K) {
            if (!q[nextEventIndex].first) {
                cerr << "[Error] Null pointer at external arrival index " << nextEventIndex << endl;
                exit(1);
            }
            auto customer = q[nextEventIndex].first;
            int customerType = customer->getType();

            auto server = customerRouteExternal(*customer);
            auto queue = server->getQueues()[0];
            queue->joinQueue(customer);

            eventComplete(nextEventIndex);
            timeCheckInd(timeHorizon);
            N[customerType]++;

            if (server->getStatus() == 0) {
                queue->commitQueue();
                server->commitNow(customer.get(), T);
                double serviceTime = service(server->getMeans()[0], server->getType());
                q[customerType + K].first = customer;
                q[customerType + K].second = T + serviceTime;
            }

            clockTick(customerType);
        }

        else if (nextEventIndex < 2 * K) {
            if (!q[nextEventIndex].first) {
                cerr << "[Error] Null pointer at service completion index " << nextEventIndex << endl;
                exit(1);
            }
            CUSTOMER* customer = q[nextEventIndex].first.get();
            double departureTime = q[nextEventIndex].second;
            int customerType = customer->getType();

            auto server = customer->getServerPointer();
            if (!server) {
                cerr << "[Error] Null SERVER pointer from customer." << endl;
                exit(1);
            }

            server->releaseNow(departureTime);
            auto queues = server->getQueues();

            auto nextServer = customerRouteTransit(customerType);
            if (nextServer != nullptr) {
                int reenterIndex = nextServer->getType() + 2 * K;
                int nextType = nextServer->getType();

                if (q[reenterIndex].second == -1) {
                    auto reCustomer = make_shared<CUSTOMER>(nextType);
                    q[reenterIndex].first = reCustomer;
                    q[reenterIndex].second = departureTime;
                    reCustomer->setArrival(departureTime);
                }
            }

            eventComplete(nextEventIndex);
            timeCheckInd(timeHorizon);
            M[server->getType()]++;

            if (!queues[0]->isEmpty()) {
                shared_ptr<CUSTOMER> newCustomer = queues[0]->commitQueue();  // use shared_ptr
                server->commitNow(newCustomer.get(), T);                      // pass raw pointer if needed
                double serviceTime = service(server->getMeans()[0], server->getType());
                q[server->getType() + K].first = newCustomer;                 // store shared_ptr directly
                q[server->getType() + K].second = T + serviceTime;
            }
        }

        else {
            if (!q[nextEventIndex].first) {
                cerr << "[Error] Null pointer at internal reentry index " << nextEventIndex << endl;
                exit(1);
            }
            auto customer = q[nextEventIndex].first;
            int customerType = customer->getType();
            auto server = Servers[customerType];
            auto queue = server->getQueues()[0];
            queue->joinQueue(customer);

            eventComplete(nextEventIndex);
            timeCheckInd(timeHorizon);
            N[customerType]++;

            if (server->getStatus() == 0) {
                queue->commitQueue();
                server->commitNow(customer.get(), T);
                double serviceTime = service(server->getMeans()[0], server->getType());
                q[server->getType() + K].first = customer;
                q[server->getType() + K].second = T + serviceTime;
            }
        }
    }

    return make_tuple(T, N, M, timeRecord, result);
}


void writeVector2D(const vector<vector<int>>& vec, const string& filename) {
    ofstream file(filename, ios::app);  // append mode
    file << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        file << "[";
        for (size_t j = 0; j < vec[i].size(); ++j) {
            file << vec[i][j];
            if (j < vec[i].size() - 1) file << ", ";
        }
        file << "]";
        if (i < vec.size() - 1) file << ", ";
    }
    file << "]\n";
    file.close();
}

string buildFilename(const string& prefix,
                     const vector<double>& rho,
                     const vector<double>& gengammaE,
                     const vector<double>& gengammaS) {
    ostringstream ss;
    ss << prefix << "[";
    for (size_t i = 0; i < rho.size(); ++i) {
        ss << fixed << setprecision(2) << rho[i];
        if (i < rho.size() - 1) ss << ", ";
    }
    ss << "]_[";
    for (size_t i = 0; i < gengammaE.size(); ++i) {
        ss << fixed << setprecision(2) << gengammaE[i];
        if (i < gengammaE.size() - 1) ss << ", ";
    }
    ss << "]_[";
    for (size_t i = 0; i < gengammaS.size(); ++i) {
        ss << fixed << setprecision(2) << gengammaS[i];
        if (i < gengammaS.size() - 1) ss << ", ";
    }
    ss << "].txt";
    return ss.str();
}

// Main write block
#include <iomanip>  // for std::fixed, std::setprecision
#include <sstream>  // for std::ostringstream


// Main block
int main() {
    cout << "Start time:" << endl;
    auto start = system_clock::now();
    time_t start_time = system_clock::to_time_t(start);
    cout << ctime(&start_time) << endl;

    result.clear();
    rho.resize(2);
    rho[0] = 0.92;
    rho[1] = 0.98;

    double time = 4000000000;
    double Horizon_t = time;

    gengammaE = {0.2, 0.45};
    gengammaS = {0.5, 0.4};

    cout << "starting rho: ";
    for (auto v : rho) cout << v << " ";
    cout << endl;

    cout << "Horizon_t: " << Horizon_t << endl;

    nominalArrivalRate = {1, 1};
    transitionMatrix = {{0.3, 0.6}, {0.4, 0.2}};
    K = nominalArrivalRate.size();

    externalArrivalRate.resize(K, 0.0);
    for (int i = 0; i < K; ++i) {
        externalArrivalRate[i] = nominalArrivalRate[i];
        for (int j = 0; j < K; ++j) {
            externalArrivalRate[i] -= transitionMatrix[j][i] * nominalArrivalRate[j];
        }
    }

    cout << "externalArrivalRate: ";
    for (auto x : externalArrivalRate) cout << x << " ";
    cout << endl;

    serviceRate.resize(K);
    for (int i = 0; i < K; ++i) {
        serviceRate[i] = nominalArrivalRate[i] / rho[i];
    }

    cout << "serviceRate: ";
    for (auto x : serviceRate) cout << x << " ";
    cout << endl;

    cout << "rho: ";
    for (auto x : rho) cout << x << " ";
    cout << endl;

    cout << "gengammaE: ";
    for (auto x : gengammaE) cout << x << " ";
    cout << endl;

    cout << "gengammaS: ";
    for (auto x : gengammaS) cout << x << " ";
    cout << endl;

    serviceTimes.resize(K);
    threshold.resize(K);
    for (int i = 0; i < K; ++i) {
        serviceTimes[i] = 1.0 / serviceRate[i];
    }

    threshold[0] = max(350, geom_ppf(0.9999, 1 - rho[0], -1));
    threshold[1] = max(700, geom_ppf(0.99999, 1 - rho[1], -1));

    // timeRecord.resize(K);
    // for (int k = 0; k < K; ++k) {
    //     timeRecord[k].resize(threshold[k] + 1, 0.0);
    // }

    initializeOccurRecords();

    N.assign(K, 0);
    M.assign(K, 0);

    q.resize(3 * K, make_pair(nullptr, -1.0));
    for (int i = 0; i < 3 * K; ++i) {
        q[i].first = nullptr;
        q[i].second = -1;
    }

    Queues.clear();
    for (int i = 0; i < K; ++i) {
        Queues.push_back(make_shared<FIFOQUEUE>(serviceTimes[i]));
    }

    Servers.clear();
    for (int k = 0; k < K; ++k) {
        Servers.push_back(make_shared<SERVER>(k));
        Servers[k]->setRates({serviceRate[k]});
        Queues[k]->setServer(Servers[k].get());
        Servers[k]->setQueue(Queues[k], 0);
    }

    auto roundResult = simulate(Horizon_t);
    result = get<4>(roundResult);

    // Write occurRecord and occurRecord_joint to match Python output
    string fname1 = buildFilename("repre_ind_ind_gen_queue_", rho, gengammaE, gengammaS);       //marginal
    string fname2 = buildFilename("repre_joint_ind_gen_queue_", rho, gengammaE, gengammaS);     //joint

    writeVector2D(occurRecord, fname1);
    writeVector2D(occurRecord_joint, fname2);

    return 0;
}