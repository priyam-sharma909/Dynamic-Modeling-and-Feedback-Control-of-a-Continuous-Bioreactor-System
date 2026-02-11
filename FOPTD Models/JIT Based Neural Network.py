# JIT_NN_bioreactor.py
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor


# SSE, AIC, BIC

def compute_sse_aic_bic(err, n_params):
    r = err.reshape(-1, 1)
    n_samples, n_outputs = r.shape
    n_eff = n_samples * n_outputs
    sse = np.sum(r ** 2)
    sigma2 = sse / n_eff
    sigma2 = max(sigma2, 1e-12)
    ll = -0.5 * n_eff * (np.log(2 * np.pi * sigma2) + 1.0)
    aic = 2 * n_params - 2 * ll
    bic = n_params * np.log(n_eff) - 2 * ll
    return sse, aic, bic

# autocorrelation function

def compute_acf(e, max_lag=40):
    r = e - np.mean(e)
    denom = np.sum(r ** 2)
    if denom < 1e-12:
        return np.zeros(max_lag + 1)
    acf_vals = []
    for lag in range(max_lag + 1):
        if lag == 0:
            acf_vals.append(1.0)
        else:
            num = np.sum(r[lag:] * r[:-lag])
            acf_vals.append(num / denom)
    return np.array(acf_vals)

# Loading bioreactor data 
mat = loadmat("bioreactor_data.mat")
t_vec = mat["t_vec"].squeeze()           # (N,)
x_nl_noisy = mat["x_nl_noisy"]           # (N,3)
D0 = float(mat["D0"].squeeze())
deltaD = float(mat["deltaD"].squeeze())

# creating D_vec and deviation u 
t = t_vec
D_vec = D0 * np.ones_like(t)
D_vec[t >= 100] = D0 - 1 * deltaD
D_vec[t >= 200] = D0 - 2 * deltaD
D_vec[t >= 300] = D0 - 3 * deltaD
u = (D_vec - D0).reshape(-1, 1)          # (N,1)

# Choosing prediction target

# predict biomass X (state 1) index = 0. 
# for S index = 1
# for P index = 2
target_state_index = 0
y_all = x_nl_noisy[:, target_state_index]    # (N,)

X_all = np.column_stack((t_vec.reshape(-1,1), u))  # shape (N,2)


# Train/Validation split (70/30), shuffle

X_tr_raw, X_va_raw, y_tr_raw, y_va_raw = train_test_split(
    X_all, y_all, test_size=0.3, random_state=42, shuffle=True
)

# Scaling

sx = StandardScaler()
sy = StandardScaler()

X_tr = sx.fit_transform(X_tr_raw)
X_va = sx.transform(X_va_raw)

y_tr_s = sy.fit_transform(y_tr_raw.reshape(-1, 1)).ravel()   # scaled target for training local models

# Local KNN + 1-hidden-layer MLP 

k = 20     # number of neighbors
h = 1     # hidden layers in local MLP

nnbrs = NearestNeighbors(n_neighbors=k, metric="euclidean")
nnbrs.fit(X_tr)

y_va_pred = np.zeros_like(y_va_raw, dtype=float)

for i in range(len(X_va)):
    q = X_va[i].reshape(1, -1)
    _, idx = nnbrs.kneighbors(q, return_distance=True)
    idx = idx.flatten()
    X_loc = X_tr[idx]                   # (k, n_features)
    y_loc = y_tr_s[idx].ravel()         # (k,)

    mlp = MLPRegressor(
        hidden_layer_sizes=(h,),
        activation="relu",
        solver="adam",
        alpha=1e-2,
        max_iter=500,
        random_state=42
    )
    # train on local neighborhood (scaled y)

    mlp.fit(X_loc, y_loc)
    # predict for query 
    # mlp expects input in training feature scaling

    y_pred_s = mlp.predict(q).reshape(1, 1)   # scaled
    y_pred_orig = sy.inverse_transform(y_pred_s)  # back to original units
    y_va_pred[i] = y_pred_orig[0, 0]


# Metrics

mse = mean_squared_error(y_va_raw, y_va_pred)
r2 = r2_score(y_va_raw, y_va_pred)

print("\n JIT (NN) — Bioreactor (predict X) ====")
print(f"MSE  : {mse:.6e}")
print(f"R²   : {r2:.4f}")

# residuals
res = y_va_raw - y_va_pred

# number of parameters for a 1-hidden-layer MLP: n_feat*h + h + h*n_out + n_out
n_feat = X_all.shape[1]
n_out = 1
n_params = n_feat * h + h + h * n_out + n_out

sse, aic, bic = compute_sse_aic_bic(res, n_params)
print("\nJIT(NN) Information criteria (validation):")
print(f"  n_params: {n_params}")
print(f"  SSE: {sse:.6f}, AIC: {aic:.2f}, BIC: {bic:.2f}")


# Plots: actual vs predicted for validation samples 
plt.figure(figsize=(12,5))
plt.plot(y_va_raw, label="Actual X (validation)", linewidth=2)
plt.plot(y_va_pred, label="Predicted X (local MLP)", linewidth=2, alpha=0.8)
plt.title(" JIT (NN)— Biomass X")
plt.xlabel("Validation sample index")
plt.ylabel("Biomass X (g/L)")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(12,4))
plt.plot(res, linewidth=1.5)
plt.title("JIT (NN) — Residuals (validation)")
plt.xlabel("Validation sample index")
plt.ylabel("ResiduaL)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# Residual ACF

max_lag = 30
lags = np.arange(max_lag + 1)
acf_vals = compute_acf(res, max_lag=max_lag)

print("\n JIT(NN) Residual autocorrelation (first few lags):")
print(acf_vals[:5])

plt.figure(figsize=(8,4))
(markerline, stemlines, baseline) = plt.stem(lags, acf_vals)

plt.axhline(0, linewidth=1)
plt.title("JIT(NN) — Residual ACF")
plt.xlabel("Lag")
plt.ylabel("Autocorrelation")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()
