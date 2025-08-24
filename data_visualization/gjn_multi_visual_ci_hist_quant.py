## This is a visualization of GJN queuing system 
# %%
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import statistics 
from scipy.stats import geom
from scipy.stats import expon
from scipy import stats
from datetime import datetime

inf = 1000000000000000000000


###################################################################################################
# %%
def plot_confidence_interval(x, values, z=1.96, color='blue', horizontal_line_width=0.25, marker='.'):
    mean = statistics.mean(values)
    stdev = statistics.stdev(values)
    confidence_interval = z * stdev / math.sqrt(len(values))

    left = x - horizontal_line_width 
    top = mean + confidence_interval
    right = x + horizontal_line_width 
    bottom = mean - confidence_interval
    plt.plot([x, x], [top, bottom], color=color)
    plt.plot([left, right], [top, top], color=color)
    plt.plot([left, right], [bottom, bottom], color=color)
    plt.plot(x, mean, marker, color='red')
    return mean, confidence_interval

def plot_confidence_interval_label(x, values, z=1.96, color='blue', horizontal_line_width=0.25, marker='.', label='ci'):
    mean = statistics.mean(values)
    stdev = statistics.stdev(values)
    confidence_interval = z * stdev / math.sqrt(len(values))

    left = x - horizontal_line_width / 2
    top = mean + confidence_interval
    right = x + horizontal_line_width / 2
    bottom = mean - confidence_interval
    plt.plot([x, x], [top, bottom], color=color)
    plt.plot([left, right], [top, top], color=color, label=label)
    plt.plot([left, right], [bottom, bottom], color=color)
    plt.plot(x, mean, marker, color='red', label='mean of estimates')

    return mean, confidence_interval


###################################################################################################
# Static parameter setting
###################################################################################################
# %%
rho = [0, 0]
######################### MAKE SURE TO UPDATE ##########################

# # copy data obstained from simulation here 
# result = 
 
time = 1000000000
rho[0] = 0.95
rho[1] =  0.9975
obs = 20

# gengammaE = [0.75, 0.8]
# gengammaS = [0.95, 0.6]

# gengammaE = [0.4, 0.5]
# gengammaS = [0.65, 0.5]

gengammaE = [0.2, 0.45]
gengammaS = [0.5, 0.4]

#########################################################################
Horizon_t = time
print('starting rho: ', rho)
print('Horizon_t:', Horizon_t)
nominalArrivalRate = [1, 1]
transitionMatrix = [[0.3, 0.6], [0.4, 0.2]]

K = len(nominalArrivalRate)
externalArrivalRate = np.matmul(np.identity(K) - np.transpose(transitionMatrix), nominalArrivalRate) 
print('externalArrivalRate: ', externalArrivalRate)

serviceRate = [nominalArrivalRate[i]/rho[i] for i in range(K)]
#serviceRate = [nominalArrivalRate[i]+r for i in range(K)]
print('serviceRate' + str(serviceRate))
#rho = [nominalArrivalRate[k]/(serviceRate[k]) for k in range(K)]
print("rho: ", rho)
print('gengammaE: ', gengammaE)
print('gengammaS: ', gengammaS)
# static global variables -------------------------------------------------------------------------
serviceTimes = [1/serviceRate[i] for i in range(K)]  
threshold = [int(len(result[0][0]))-1, int(len(result[0][1]))-1]
timeRecord = [[0] * (threshold[k] +1) for k in range(K)]

###################################################################################################
# Simulation result -----------------------------------------------------------------------
###################################################################################################
#result = []

##########################################################################################################
# Theoretical estimation  ------------------------------------------------------------------
###################################################################################################
# %%
wMatrix = [[0]* K for x in range(K)]
d = [0 for x in range(K)]
d_a = [0 for x in range(K)]
wMatrix[0][0] = transitionMatrix[0][0]
wMatrix[1][0] = transitionMatrix[1][0]
wMatrix[0][1] = transitionMatrix[0][1]/(1-transitionMatrix[0][0])
wMatrix[1][1] = transitionMatrix[1][1] + (transitionMatrix[0][1] * transitionMatrix[1][0]) / (1-transitionMatrix[0][0])
c2e = [1/i for i in gengammaE]
c2s = [1/i for i in gengammaS]
d[0] = 0.5 * (externalArrivalRate[0]* c2e[0] + nominalArrivalRate[1] *(wMatrix[1][0]**2*c2s[1] + wMatrix[1][0]*(1-wMatrix[1][0]))+nominalArrivalRate[0]*(c2s[0]*(1-wMatrix[0][0])**2 + wMatrix[0][0]*(1 - wMatrix[0][0]))) / (1 - wMatrix[0][0])
d[1] = 0.5 * (externalArrivalRate[0] * (wMatrix[0][1]**2 * c2e[0] + wMatrix[0][1]*(1-wMatrix[0][1])) + externalArrivalRate[1]*c2e[1] + nominalArrivalRate[1] * ((1- wMatrix[1][1])**2 *c2s[1] + (1- wMatrix[1][1])*wMatrix[1][1])) / (1- wMatrix[1][1])
for k in range(K):
    d_a[k] = rho[k] * d[k]

print("d: ", d)
print('d/(mu-lambda)[MEAN]:', [d[k]/(serviceRate[k]-nominalArrivalRate[k]) for k in range(K)])


######################################### (scaled)######################################################
# one confidence interval for exponential distirbution
################################################################################################
# %%
timecut = Horizon_t/(10 * obs)
scale_avg_queue = [[], []]

for i in range(obs): #add serviceRate[k]
    for k in range(K):
        scale_avg_q = sum([result[i][k][j] * j* (serviceRate[k]-nominalArrivalRate[k]) for j in range(threshold[k]+1)]) / timecut
        scale_avg_queue[k].append(scale_avg_q)

value_ticks = []
for k in range(K):
    print('class: ', k+1)
    print(plot_confidence_interval(1+0.5* k, scale_avg_queue[k]))
    plt.plot(1.1+0.5* k, d[k], 'o', color='orange')
    plt.plot(1.1+0.5* k, d_a[k], 'o', color='green')   
    #plt.plot(1.1+0.5* k, mean_discrete[k]*(serviceRate[k]-nominalArrivalRate[k]), 'o', color='purple') 
    value_ticks.append(1 +0.5* k)

plt.plot(1.1+0.5* k, d[k], 'o', color='orange', label='d_limit')
plt.plot(1.1+0.5* k, d_a[k], 'o', color='green', label='d_adjusted')
plot_confidence_interval_label(1+0.5* k, scale_avg_queue[k])

plt.legend(loc='upper left')
plt.xticks(value_ticks, ['station-1','station-2'])
plt.ylabel('Mean queue length')
plt.title('Mean of scaled queue lengths')
  
plt.show()

########################################### (unscaled) ##############################################
#  # one confidence interval for exponential distirbution 
####################################################################################################
# %%
timecut = Horizon_t/(10 * obs)
avg_queue = [[], []]

for i in range(obs): #add serviceRate[k]
    for k in range(K):
        avg_q = sum([result[i][k][j] * j for j in range(threshold[k]+1)]) / timecut
        avg_queue[k].append(avg_q)

value_ticks = []
for k in range(K):
    print('class: ', k+1)
    print(plot_confidence_interval(1+0.5* k, avg_queue[k]))
    plt.plot(1.1+0.5* k, d[k]/(serviceRate[k]-nominalArrivalRate[k]), 'o', color='orange')
    value_ticks.append(1 +0.5* k)

plt.plot(1.1+0.5* k, d[k]/(serviceRate[k]-nominalArrivalRate[k]), 'o', color='orange', label='d_limit')
plot_confidence_interval_label(1+0.5* k, avg_queue[k])

plt.legend(loc='upper left')
plt.xticks(value_ticks, ['station-1','station-2'])
plt.ylabel('Mean queue length')
plt.title('Mean of scaled queue lengths')
  
plt.show()

# ##########################################################################################################
# # Histogram  ------------------------------------------------------------------
# ###################################################################################################
# %%
# exponential visualization
def plot_confidence_interval_label(x, values, z=1.96, color='blue', horizontal_line_width=0.5, marker='.'):
    mean = statistics.mean(values)
    stdev = statistics.stdev(values)
    confidence_interval = z * stdev / math.sqrt(len(values))

    left = x - horizontal_line_width / 2
    top = mean + confidence_interval
    right = x + horizontal_line_width / 2
    bottom = mean - confidence_interval
    plt.plot(x, mean, marker, color='red', label='Sim')
    plt.plot([left, right], [bottom, bottom], color=color, label='Sim CI')
    return mean, confidence_interval

output = [[], []]
for k in range(K):
    rg = threshold[k]
    output[k] = [[] for i in range(rg)]

    for r in range(rg):
        output[k][r] = [result[i][k][r]/timecut for i in range(obs)]

# find the last nonzero element
last_nonzero_r_list = []
for k in range(K):
    rg = threshold[k]
    last_nonzero_r = -1  # default in case all are zero

    for r in reversed(range(rg)):
        has_nonzero = False
        for i in range(obs):
            if result[i][k][r] != 0:
                has_nonzero = True
                break
        if has_nonzero:
            last_nonzero_r = r
            break

    last_nonzero_r_list.append(last_nonzero_r)

lambda_0 = [0, 0]
mean_discrete = [0, 0]
for k in range(K):
    rg = threshold[k]
    c = 0
    wid = 0.5
    lambda_0[k] = rho[k]/(d[k]* (1-expon.cdf(c*(1-rho[k]),loc = 0, scale=d[k])))
    pointsa = np.linspace(c, rg+c, rg+1)
    points2a = np.linspace(1+ c, rg+1 + c, rg+1)
    cdf1 = np.append(1-rho[k], lambda_0[k] * d[k] * np.delete(expon.cdf(points2a*(1-rho[k])* serviceRate[k],loc = 0, scale=d[k]) - expon.cdf(pointsa*(1-rho[k])* serviceRate[k],loc = 0, scale=d[k]),-1) )
    plt.bar(pointsa-c-wid/2, cdf1, width=wid, color='orange', alpha=0.8, edgecolor='orange',align='edge', label='M-Scale') 


    m=0
    for r in range(rg):
        r_mean = plot_confidence_interval(r, output[k][r])[0]
        m +=  r * r_mean
        #print('r=', r, 'mean=', r_mean)
    plot_confidence_interval_label(r, output[k][r])
    print('scaled queue mean:',  m)
    mean_discrete[k] = np.dot(points2a-1-c, cdf1)
    print('scaled exp mean: ', np.dot(points2a-1-c, cdf1))
    plt.xlabel('Queue length')
    plt.xticks(np.arange(0,rg,5))
    plt.ylabel('Percentage')
    #plt.title('Histogram: station-'+str(k+1) +' with rho=' + str(np.round(rho[k], 4)*100)+'%')
    plt.legend(loc='upper right')

    # if k == 0:
    #     plt.xlim(0, last_nonzero_r_list[0])
    # else:
    #     plt.xlim(0, last_nonzero_r_list[1])

    if k ==0:
        plt.xlim(-1,20)
        plt.ylim(0, )
    else:
        plt.xlim(-1,20)
        plt.ylim(0, )

    # if k ==0:
    #     plt.xlim(40,60)
    #     plt.ylim(0,0.01)
    # else:
    #     plt.xlim(90,110)
    #     plt.ylim(0,0.01)

    # if k ==0:
    #     plt.xlim(70,90)
    #     plt.ylim(0,0.01)
    # else:
    #     plt.xlim(70,90)
    #     plt.ylim(0,0.01)

    # if k ==0:
    #     plt.xlim(100,120)
    #     plt.ylim(0,0.002)
    # else:
    #     plt.xlim(400,420)
    #     plt.ylim(0,0.003)        
    plt.show()

