# Minsoo Lee
# This program will conduct cross-dependence (CD) testing for environmental data across 43
# countries between 1990 to 2019, and tests the strength of cross dependence


import pandas as pd
import numpy as np
from linearmodels.panel import PanelOLS
from scipy.stats import norm
from statsmodels.stats.outliers_influence import variance_inflation_factor
import matplotlib.pyplot as plt

# Load data
df = pd.read_excel('Cabon_Emm_disaggregate.xlsx')

# Confirm N and T values
df['Year'] = pd.to_datetime(df['Year']).dt.year

# Two-way fixed effects model, consumption only as we are measuring carbon emissions
df = df.set_index(['Country.Name', 'Year']).sort_index() # set country name (i) and year (t) as index

# Estimate model
y = df['carbon_emissions_pc']

# Combined
X = df[[
    'gdp_pc',
    'gdpsq_pc',
    'fd_index',
    'coal_consumption_pc',
    'oil_consumption_pc',
    'gas_consumption_pc',
    'renewables_consumption_pc',
    'coal_production_pc',
    'oil_production_pc',
    'gas_production_pc',
    'solar_production_pc',
    'wind_production_pc',
    'hydro_production_pc',
    'other_renewables_production_pc' ]]

# VIF Test for multicollinearity

# Demean by entity and time
X_demeaned = X.copy()
X_demeaned = X_demeaned - X_demeaned.groupby(level=0).transform('mean') # entity
X_demeaned = X_demeaned - X_demeaned.groupby(level=1).transform('mean') # time

# Create dataframe, and calculate VIF for each regressor
vif_data = pd.DataFrame()
vif_data["variable"] = X_demeaned.columns
vif_data["VIF"] = [
    variance_inflation_factor(X_demeaned.values, i)
    for i in range(X_demeaned.shape[1])
]

print(vif_data)

# Results show low multicollinearity, but this might be due to high heterogeneity in the data (Australia produces more
# coal that it consumes, Japan consumes more coal than it produces)

# Fit model using OLS
model = PanelOLS(y, X, entity_effects=True, time_effects=True)
results = model.fit(cov_type='clustered', cluster_entity=True)

print(results.summary)

resid = results.resids.reset_index()
resid.columns = ['Country.Name', 'Year', 'resid']

# Reshape residuals to wider form for correlation matrix
resid_wide = resid.pivot(index='Year', columns='Country.Name', values='resid').sort_index()

# Correlation matrix across countries
corr_matrix = resid_wide.corr()

# Determine N countries and T years for calculation in Pesaran CD statistic
N = resid_wide.shape[1]
T = resid_wide.shape[0]

# Upper triangle correlations only
upper_idx = np.triu_indices(N, k=1)
rho_ij = corr_matrix.to_numpy()[upper_idx]

# drop undefined values
# rho_ij = rho_ij[~np.isnan(rho_ij)]

# Pesaran CD statistic
CD = np.sqrt(2 * T / (N * (N - 1))) * np.sum(rho_ij)

# Two-sided p-value to evaluate test statistic
p_value = 2 * (1 - norm.cdf(abs(CD)))

print(corr_matrix.iloc[:5, :5])

plt.figure(figsize=(10, 6))

for col in resid_wide.columns:
    plt.plot(resid_wide.index, resid_wide[col], label=col)

plt.title("Model Residuals across countries")
plt.xlabel("Year")
plt.ylabel("Residual")

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', ncol=2, fontsize=8)
plt.tight_layout()
plt.show()
print("N =", N)
print("T =", T)
print("CD statistic =", CD)
print("p-value =", p_value)

# Test for strength of cross-dependence (Bailey et al 2019)

# Full correlation matrix as numpy array
R = corr_matrix.to_numpy()

# Threshold choice
threshold = 2 * np.sqrt(np.log(N)) / np.sqrt(T)

# Construct thresholded correlation matrix represented by delta with diagonal of 1
delta = np.zeros((N, N))
np.fill_diagonal(delta, 1)

for i in range(N):
    for j in range(N):
        if i != j:
            if abs(R[i, j]) > threshold:
                delta[i, j] = R[i, j]
            else:
                delta[i, j] = 0

# Construct tau as vector of ones
tau = np.ones((N, 1))

# Quadratic form of tau used in computation of alpha statistic
quad_form = (tau.T @ delta @ tau).item()

# Calculate measure of strength of cross-dependence
alpha = np.log(quad_form) / (2 * np.log(N))

print("Estimated alpha =", alpha)
