# Minsoo Lee
# This program will conduct cross-dependence (CD) testing for environmental data across 43
# countries between 1990 to 2019


import pandas as pd
import numpy as np
from linearmodels.panel import PanelOLS
from scipy.stats import norm
import matplotlib.pyplot as plt

# Load data
df = pd.read_excel('Cabon_Emm_disaggregate.xlsx')

# Confirm N and T values
df['Year'] = pd.to_datetime(df['Year']).dt.year

# Two-way fixed effects model, consumption only as we are measuring carbon emissions
df = df.set_index(['Country.Name', 'Year']).sort_index() # set country name (i) and year (t) as index

# Estimate model
y = df['carbon_emissions_pc']

# These variables have been split into production and consumption side to avoid multicollinearity
X = df[[
    'gdp_pc',
    'gdpsq_pc',
    'fd_index',
    'coal_consumption_pc',
    'oil_consumption_pc',
    'gas_consumption_pc',
    'renewables_consumption_pc'
]] # Dummy variables are absorbed



# Fit model using OLS
model = PanelOLS(y, X, entity_effects=True, time_effects=True)
results = model.fit(cov_type='clustered', cluster_entity=True)

print(results.summary)

resid = results.resids.reset_index()
resid.columns = ['Country.Name', 'Year', 'resid']

# reshape residuals to wider form for correlation matrix
resid_wide = resid.pivot(index='Year', columns='Country.Name', values='resid').sort_index()

# correlation matrix across countries
corr_matrix = resid_wide.corr()

# Determine N countries and T years for calculation in Pesaran CD statistic
N = resid_wide.shape[1]
T = resid_wide.shape[0]

# upper triangle correlations only
upper_idx = np.triu_indices(N, k=1)
rho_ij = corr_matrix.to_numpy()[upper_idx]

# drop undefined values
# rho_ij = rho_ij[~np.isnan(rho_ij)]

# Pesaran CD statistic
CD = np.sqrt(2 * T / (N * (N - 1))) * np.sum(rho_ij) # i think this is correct?

# two-sided p-value to evaluate test statistic
p_value = 2 * (1 - norm.cdf(abs(CD)))

print(corr_matrix.iloc[:5, :5])

# plot a subset of the residuals
avg_resid = resid_wide.mean(axis=1)

plt.plot(avg_resid)
plt.title("Average residual across countries")
plt.xlabel("Year")
plt.ylabel("Residual")
plt.show()

print("N =", N)
print("T =", T)
print("CD statistic =", CD)
print("p-value =", p_value)












