## This is visualization of GJN independence
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
from scipy.stats import chi2_contingency
import json
import pandas as pd
from scipy.stats import chi2

inf = 1000000000000000000000

################################################################################
# read data
#  %%
with open('/Users/sally/Documents/codes/GJN/simulated_data/repre_joint_ind_gen_queue_[0.80, 0.96]_[0.20, 0.45]_[0.50, 0.40].txt', 'r') as file:
    content1 = file.read()
result_joint = json.loads(content1)

with open('/Users/sally/Documents/codes/GJN/simulated_data/repre_ind_ind_gen_queue_[0.80, 0.96]_[0.20, 0.45]_[0.50, 0.40].txt', 'r') as file:
    content2 = file.read()
result = json.loads(content2)

################################################################
# one confidence interval for exponential distirbution
# %%
s1 = 10 # Number of rows 
s2 = 10 # Number of columns 
s = 10

prob = [0, 0]

result_joint = np.array(result_joint)
total_joint= result_joint.sum()  # Total sum

#change this to specify the point
K = 2
for i in range(6):
    for j in range(6):
        point = [j, i]
        for k in range(K):
            temp = np.array(result[k])
            total = temp.sum()
            prob[k] = temp[point[k]]/total
        prob_joint = result_joint[point[0]][point[1]]/total_joint
        print('p(z1=',point[0],')=', prob[0])
        print('p(z2=',point[1],')=', prob[1])

        print('p(z1=',point[0], ')*','p(z2=',point[1],')=', np.round(prob[0] * prob[1],5)*1000)
        print('p(z1=',point[0],',z2=',point[1],')=', np.round(prob_joint,5)*1000)
        print('----')
    print('***********')

################################################################################################
##################### chi-square test for whole batch ################################
################################################################################################
# %%

print('chi-square test for whole batch')
# Assuming `result_joint` is defined as a list of joint distributions
obs = len(result_joint)  # Number of elements in the result_joint
s= [len(result_joint), len(result_joint[0])]
result_joint = np.array(result_joint)

# Example Input: Joint and Marginal Distributions
row_marginals = result_joint.sum(axis=1)  # Row sums
col_marginals = result_joint.sum(axis=0)  # Column sums
total = result_joint.sum()  # Total sum

# Step 1: Compute Expected Frequencies
expected = np.outer(row_marginals, col_marginals) / total

# Step 2: Drop terms with zero observations or expected frequencies
# Only keep terms where both result_joint_sum and expected > 0
mask = (expected > 0)
filtered_result_joint = result_joint[mask]
filtered_expected = expected[mask]

# Step 3: Compute the Chi-Square Statistic
chi2_stat = np.sum((filtered_result_joint - filtered_expected) ** 2 / filtered_expected)

# Step 4: Degrees of Freedom
# Exclude rows and columns with zero sums
row_sums = result_joint.sum(axis=1)
col_sums = result_joint.sum(axis=0)

valid_rows = row_sums > 0
valid_cols = col_sums > 0

filtered_result_df = result_joint[valid_rows, :][:, valid_cols]
df = (filtered_result_df.shape[0] - 1) * (filtered_result_df.shape[1] - 1)
print(filtered_result_df.shape)
# Step 5: Compute p-value from Chi-Square distribution
p_value = 1 - chi2.cdf(chi2_stat, df)

# Output the Results
print(f"Chi-Square Statistic: {chi2_stat}")
print(f"Degrees of Freedom: {df}")
print(f"P-Value: {p_value}")



################################################################################################
################################# package chi-square test for whole batch ################################
################################################################################################
# %%
print('package chi-square test for whole batch')

# Assuming `result_joint` is defined as a list of joint distributions
obs = len(result_joint)  # Number of elements in the result_joint
s= [len(result_joint), len(result_joint[0])]
result_joint = np.array(result_joint)
# Step 2: Exclude rows and columns with all zeros
row_sums = result_joint.sum(axis=1)
col_sums = result_joint.sum(axis=0)

# Filter rows and columns where the sums are greater than 0
valid_rows = row_sums > 0
valid_cols = col_sums > 0

filtered_result_joint= result_joint[valid_rows, :][:, valid_cols]
print(filtered_result_joint.shape)
# Step 3: Check if the filtered table is still valid
if filtered_result_joint.size == 0 or filtered_result_joint.shape[0] < 2 or filtered_result_joint.shape[1] < 2:
    raise ValueError("The contingency table has insufficient data after filtering. Cannot perform Chi-Square test.")

# Step 4: Use chi2_contingency to perform the Chi-Square test
chi2_stat, p_value, df, expected = chi2_contingency(filtered_result_joint)

# Step 5: Output the Results
print(f"Chi-Square Statistic: {chi2_stat}")
print(f"Degrees of Freedom: {df}")
print(f"P-Value: {p_value}")
#print(f"Expected Frequencies:\n{expected}")



################################################################################################
########### likelyhood ratio test for whole batch ################################
################################################################################################
# %%
print('likelyhood ratio test for whole batch')

# Assuming `result_joint` is defined as a list of joint distributions
obs = len(result_joint)  # Number of elements in the result_joint
s= [len(result_joint), len(result_joint[0])]
result_joint = np.array(result_joint)

row_marginals = result_joint.sum(axis=1)  # Row sums
col_marginals = result_joint.sum(axis=0)  # Column sums
total = result_joint.sum()  # Total sum

# Step 1: Compute Expected Frequencies
expected = np.outer(row_marginals, col_marginals) / total

# Step 2: Drop terms with zero observations or expected frequencies
# Only keep terms where both result_joint_sum and expected > 0
mask = (result_joint > 0) & (expected > 0)
filtered_result_joint= result_joint[mask]
filtered_expected = expected[mask]

# Step 3: Compute the G-statistic
G_stat = 2 * np.sum(filtered_result_joint * np.log(filtered_result_joint / filtered_expected))

# Step 4: Degrees of Freedom
# Exclude rows and columns with zero sums
row_sums = result_joint.sum(axis=1)
col_sums = result_joint.sum(axis=0)

valid_rows = row_sums > 0
valid_cols = col_sums > 0

filtered_result_df = result_joint[valid_rows, :][:, valid_cols]
df = (filtered_result_df.shape[0] - 1) * (filtered_result_df.shape[1] - 1)
print(filtered_result_df.shape)
# Step 5: Compute p-value from Chi-Square distribution
p_value = 1 - chi2.cdf(G_stat, df)

# Output the Results
print(f"G-Statistic: {G_stat}")
print(f"Degrees of Freedom: {df}")
print(f"P-Value: {p_value}")
