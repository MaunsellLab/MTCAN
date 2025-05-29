"""
 This script predicts responses CRF responses to single or dual stimuli using a weighted RF model
 and it incorporates attention. Tries to resolve response gain vs contrast gain effect. The prediction
 shows that depending on stimulus configuration, you can have either response or contrast gain.
"""
#####################################################################################

import numpy as np
import matplotlib.pyplot as plt

# Parameters
contrasts = np.logspace(-1, 3, 100)  # contrast range: 0.01 to 1
sigma = 7
L1 = 1.0  # fixed gains

# X-axis tick positions and labels
xticks = [0.1, 1.0, 10.0, 100.0]
xtick_labels = ['0.1', '1', '10', '100']

# Attention levels: w1 increases with attention
w1_levels = [1.0, 0.5, 0.25]  # simulates unattended, neutral, attended

# Single stimulus condition
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
for w1 in w1_levels:
    R_single = (contrasts * w1 * L1) / ((contrasts * w1) + sigma)
    plt.plot(contrasts, R_single, label=f'intensity={w1}')
plt.xticks(xticks, xtick_labels)
plt.xscale('log')
plt.xlabel('Contrast')
plt.ylabel('Response')
plt.title('Single Stimulus (Attention ~ Contrast Gain)')
plt.ylim(top=1)
plt.legend()
plt.grid(True, which='both', linestyle='--')

plt.subplot(1, 2, 2)

# Dual stimulus condition
# Stimulus 2 is fixed at  w=1, L=0.5
w2, L2 = 1, 0.1
w1 = 1.0
attentionLevels = [1.0, 2.0, 3.0]

for a1 in attentionLevels:
    numerator = (contrasts * w1 * a1 * L1) + (contrasts * 1 * w2 * L2)
    denominator = (contrasts * a1 * w1) + (contrasts * 1 * w2) + sigma
    R_dual = numerator / denominator
    plt.plot(contrasts, R_dual, label=f'intensity={a1}')
plt.xticks(xticks, xtick_labels)
plt.xscale('log')
plt.xlabel('Contrast')
plt.ylabel('Response')
plt.title('Dual Stimulus (Attention ~ Response Gain)')
plt.legend()
plt.ylim(top=1.0)
plt.grid(True, which='both', linestyle='--')

plt.tight_layout()
plt.show()


################################################ heuristic model one stim ##############################################

# Parameters
contrasts = np.logspace(-3, 2, 100)  # contrast range: 0.01 to 1
iNaught = 1.0
L1 = 1.0  # fixed gains
baseline = .05

# Attention levels: w1 increases with attention
w1_levels = [1.0, 2.0, 3.0]  # simulates unattended, neutral, attended

# Single stimulus condition
plt.figure(figsize=(12, 5))
plt.subplot()
for w1 in w1_levels:
    R_single = ((contrasts * w1 * L1) + (iNaught * baseline)) / ((contrasts * w1) + iNaught)
    plt.plot(contrasts, R_single, label=f'intensity={w1}')
plt.xscale('log')
plt.xlabel('Contrast')
plt.ylabel('Response')
plt.title('Single Stimulus (Attention ~ Contrast Gain)')
plt.ylim(top=1)
plt.legend()
plt.grid(True, which='both', linestyle='--')

plt.show()

################################################ heuristic model one stim ##############################################
############################################## attention on iNaught as well ############################################

# Parameters
contrasts = np.logspace(-3, 2, 100)  # contrast range: 0.01 to 1
iNaught = .1
L1 = 1.0  # fixed gains
baseline = .05

# Attention levels: w1 increases with attention
w1_levels = [1.0, 5.0, 10.0]  # simulates unattended, neutral, attended

# Single stimulus condition
plt.figure(figsize=(12, 5))
plt.subplot()
for w1 in w1_levels:
    R_single = ((contrasts * w1 * L1) + (iNaught * w1 * baseline)) / ((contrasts * w1) + (w1 * iNaught))
    plt.plot(contrasts, R_single, label=f'intensity={w1}')
plt.xscale('log')
plt.xlabel('Contrast')
plt.ylabel('Response')
plt.title('Single Stimulus (Attention ~ Contrast Gain)')
plt.ylim(top=1)
plt.legend()
plt.grid(True, which='both', linestyle='--')

plt.show()


#####
# Common fixed parameters
L1 = 1.0
contrast = 1.0  # 100% contrast
attention_conditions = [1, 2]  # attention weights
baseline_range = np.linspace(0, 1, 100)
i0_range = np.linspace(0.01, 2.0, 100)  # sweep i0

# Helper function to compute Rmax
def compute_rmax(w1, contrast, L1, i0, baseline):
    return ((contrast * w1 * L1) + (i0 * baseline)) / ((contrast * w1) + i0)

# Setup plots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.flatten()

# === Plot 1: Sweep baseline, low i0 ===
i0_low = 0.1
for w1 in attention_conditions:
    rmax_vals = compute_rmax(w1, contrast, L1, i0_low, baseline_range)
    axs[0].plot(baseline_range, rmax_vals, label=f'Attention = {w1}')
axs[0].set_title('Low i0, Varying Baseline')
axs[0].set_xlabel('Baseline (b)')
axs[0].set_ylabel('Rmax')
axs[0].legend()
axs[0].grid(True)

# === Plot 2: Sweep baseline, high i0 ===
i0_high = 2.0
for w1 in attention_conditions:
    rmax_vals = compute_rmax(w1, contrast, L1, i0_high, baseline_range)
    axs[1].plot(baseline_range, rmax_vals, label=f'Attention = {w1}')
axs[1].set_title('High i0, Varying Baseline')
axs[1].set_xlabel('Baseline (b)')
axs[1].set_ylabel('Rmax')
axs[1].legend()
axs[1].grid(True)

# === Plot 3: Sweep i0, low baseline ===
baseline_low = 0.1
for w1 in attention_conditions:
    rmax_vals = compute_rmax(w1, contrast, L1, i0_range, baseline_low)
    axs[2].plot(i0_range, rmax_vals, label=f'Attention = {w1}')
axs[2].set_title('Low Baseline, Varying i0')
axs[2].set_xlabel('i0')
axs[2].set_ylabel('Rmax')
axs[2].legend()
axs[2].grid(True)

# === Plot 4: Sweep i0, high baseline ===
baseline_high = 1.0
for w1 in attention_conditions:
    rmax_vals = compute_rmax(w1, contrast, L1, i0_range, baseline_high)
    axs[3].plot(i0_range, rmax_vals, label=f'Attention = {w1}')
axs[3].set_title('High Baseline, Varying i0')
axs[3].set_xlabel('i0')
axs[3].set_ylabel('Rmax')
axs[3].legend()
axs[3].grid(True)

plt.tight_layout()
plt.show()

#########################################
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Fixed parameters
L1 = 1.0
contrast_vals = np.linspace(0, 1, 500)
baseline_vals = np.linspace(0, 1, 100)
i0_vals = np.linspace(0.01, 2.0, 100)
attention_conditions = [1, 2]

# Define model
def response(c, w1, L1, i0, baseline):
    return ((c * w1 * L1) + (i0 * baseline)) / ((c * w1) + i0)

def compute_rmax(w1, L1, i0, b):
    return response(1.0, w1, L1, i0, b)

def compute_c50(w1, L1, i0, b):
    r = response(contrast_vals, w1, L1, i0, b)
    rmax = r[-1]
    half = rmax / 2
    interp = interp1d(r, contrast_vals, bounds_error=False, fill_value='extrapolate')
    c50 = interp(half)
    return np.clip(c50, 0, 1)

# Store results for sweeping baseline (fixed i0)
baseline_results = {}
fixed_i0 = 0.1
for attn in attention_conditions:
    rmax_list, c50_list = [], []
    for b in baseline_vals:
        rmax_list.append(compute_rmax(attn, L1, fixed_i0, b))
        c50_list.append(compute_c50(attn, L1, fixed_i0, b))
    baseline_results[attn] = (np.array(rmax_list), np.array(c50_list))

# Store results for sweeping i0 (fixed baseline)
i0_results = {}
fixed_baseline = 0.1
for attn in attention_conditions:
    rmax_list, c50_list = [], []
    for i0 in i0_vals:
        rmax_list.append(compute_rmax(attn, L1, i0, fixed_baseline))
        c50_list.append(compute_c50(attn, L1, i0, fixed_baseline))
    i0_results[attn] = (np.array(rmax_list), np.array(c50_list))

# Plotting both subplots
fig = plt.figure(figsize=(16, 7))

# Left: Sweep baseline
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
colors = ['blue', 'orange']
for i, attn in enumerate(attention_conditions):
    rmax, c50 = baseline_results[attn]
    ax1.plot(baseline_vals, rmax, c50, color=colors[i], label=f'Attention = {attn}', lw=2)
ax1.set_xlabel('Baseline (b)')
ax1.set_ylabel('Rmax')
ax1.set_zlabel('C50')
ax1.set_title('Sweep Baseline (fixed i0 = 0.1)')
ax1.view_init(elev=25, azim=135)
ax1.legend()

# Right: Sweep i0
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
for i, attn in enumerate(attention_conditions):
    rmax, c50 = i0_results[attn]
    ax2.plot(i0_vals, rmax, c50, color=colors[i], label=f'Attention = {attn}', lw=2)
ax2.set_xlabel('i0')
ax2.set_ylabel('Rmax')
ax2.set_zlabel('C50')
ax2.set_title('Sweep i0 (fixed baseline = 0.1)')
ax2.view_init(elev=25, azim=135)
ax2.legend()

plt.tight_layout()
plt.show()

########################## 4x4 grid
import numpy as np
import matplotlib.pyplot as plt

# Parameters
L1 = 1.0
contrast_vals = np.logspace(np.log10(0.01), 0, 200)
# contrast_vals = np.linspace(0, 1, 200)
attention_conditions = [1, 2]


# Sweep values
baseline_vals = np.linspace(0.1, 1, 4)     # 4 columns
i0_vals = np.linspace(0.1, 3, 4)           # 4 rows

# Model function
def response(c, w1, L1, i0, baseline):
    return ((c * w1 * L1) + (i0 * baseline)) / ((c * w1) + i0)

# Create figure and grid of subplots
fig, axes = plt.subplots(4, 4, figsize=(8, 8), sharex=True, sharey=True)

colors = ['blue', 'orange']

for row, i0 in enumerate(i0_vals):
    for col, baseline in enumerate(baseline_vals):
        ax = axes[row, col]
        for i, attn in enumerate(attention_conditions):
            crf = response(contrast_vals, attn, L1, i0, baseline)
            ax.plot(contrast_vals, crf, label=f'Attention={attn}', color=colors[i], lw=2)
            ax.set_xscale('log')
        ax.set_title(f'i0={i0:.2f}, b={baseline:.2f}', fontsize=10)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.1)
        if row == 3:
            ax.set_xlabel('Contrast')
        if col == 0:
            ax.set_ylabel('Response')

# Add single legend to the top-left plot
axes[0, 0].legend(loc='upper left', fontsize=8)
fig.suptitle('Contrast Response Functions under Varying i0 and Baseline', fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

########################################################################################################################
########### generate CRF using Naka-Rushton and then fit my model to all curves separately
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Naka-Rushton function with fixed baseline
def naka_rushton(c, Rmax, c50, n, baseline=0.1):
    return baseline + Rmax * (c ** n) / (c ** n + c50 ** n)


# Updated normalization model with L1 as a parameter
def norm_model(c, w1, i0, b, L1):
    return ((c * w1 * L1) + (i0 * b)) / ((c * w1) + i0)


# Fit helper function
def fit_model(c_vals, responses):
    bounds_lower = [0, 0, 0, 0]
    bounds_upper = [5, 5, 5, 5]
    initial_guess = [1, 0.5, 0.1, 1]
    popt, _ = curve_fit(norm_model, c_vals, responses, bounds=(bounds_lower, bounds_upper), p0=initial_guess)
    return popt, norm_model(c_vals, *popt)


# Generate contrast values
contrast_vals = np.logspace(np.log10(0.01), 1, 100)

# Parameters
baseline = 0.1
n = 1.1
Rmax_base = 1.0
Rmax_gain = 1.5
c50_base = 0.2
c50_gain = 0.1

# Generate data
resp_base = naka_rushton(contrast_vals, Rmax_base, c50_base, n, baseline)
resp_response_gain = naka_rushton(contrast_vals, Rmax_gain, c50_base, n, baseline)
resp_contrast_gain = naka_rushton(contrast_vals, Rmax_base, c50_gain, n, baseline)

# Fit each condition
params_base_rg, fit_base_rg = fit_model(contrast_vals, resp_base)
params_gain_rg, fit_gain_rg = fit_model(contrast_vals, resp_response_gain)
params_base_cg, fit_base_cg = fit_model(contrast_vals, resp_base)
params_gain_cg, fit_gain_cg = fit_model(contrast_vals, resp_contrast_gain)

# Plotting
fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

# Response Gain
axs[0].plot(contrast_vals, resp_base, 'ko', label='Base Data')
axs[0].plot(contrast_vals, fit_base_rg, 'k-', label='Base Fit')
axs[0].plot(contrast_vals, resp_response_gain, 'bo', label='Response Gain Data')
axs[0].plot(contrast_vals, fit_gain_rg, 'b-', label='Response Gain Fit')
axs[0].set_xscale('log')
axs[0].set_title('Response Gain')
axs[0].set_xlabel('Contrast')
axs[0].set_ylabel('Response')
axs[0].legend()
axs[0].grid(True)
axs[0].text(0.65, 0.4,
    f'Base Fit:\na1={params_base_rg[0]:.2f}, i0={params_base_rg[1]:.2f}, b={params_base_rg[2]:.2f}, L1={params_base_rg[3]:.2f}\n'
    f'Gain Fit:\na1={params_gain_rg[0]:.2f}, i0={params_gain_rg[1]:.2f}, b={params_gain_rg[2]:.2e}, L1={params_gain_rg[3]:.2f}',
    transform=axs[0].transAxes, fontsize=9, verticalalignment='center')

# Contrast Gain
axs[1].plot(contrast_vals, resp_base, 'ko', label='Base Data')
axs[1].plot(contrast_vals, fit_base_cg, 'k-', label='Base Fit')
axs[1].plot(contrast_vals, resp_contrast_gain, 'go', label='Contrast Gain Data')
axs[1].plot(contrast_vals, fit_gain_cg, 'g-', label='Contrast Gain Fit')
axs[1].set_xscale('log')
axs[1].set_title('Contrast Gain')
axs[1].set_xlabel('Contrast')
axs[1].legend()
axs[1].grid(True)
axs[1].text(0.02, 0.9,
    f'Base Fit:\na1={params_base_cg[0]:.2f}, i0={params_base_cg[1]:.2f}, b={params_base_cg[2]:.2f}, L1={params_base_cg[3]:.2f}\n'
    f'Gain Fit:\na1={params_gain_cg[0]:.2f}, i0={params_gain_cg[1]:.2f}, b={params_gain_cg[2]:.2e}, L1={params_gain_cg[3]:.2f}',
    transform=axs[1].transAxes, fontsize=9, verticalalignment='top')

plt.tight_layout()
plt.show()

########################################################################################################################

# Parameters
c = np.logspace(np.log10(0.01), np.log10(10.0), 100)  # contrast values from 0.01 to 1 in log space
w = 1
L = 1
b = 0.1
sigma = 0.2
i0 = 0.8
a = 2  # attention gain

# Model 1 – Full Attention (Stimulus + Baseline, Additive)
R1_no_attn = (c * w * L) / (c * w + sigma) + b
R1_attn = (a * c * w * L) / (a * c * w + sigma) + a * b

# Model 2 – Stimulus Attention Only (No Baseline Gain)
R2_no_attn = (c * w * L) / (c * w + sigma) + b
R2_attn = (a * c * w * L) / (a * c * w + sigma) + b

# Model 3 – Normalization Model with Asymmetric Attention
R3_no_attn = (c * w * L + i0 * b) / (c * w + i0)
R3_attn = (a * c * w * L + i0 * a * b) / (a * c * w + i0)

# Model 4 – Normalization Model with Stimulus Attention Only (no i0*b attention)
R4_no_attn = (c * w * L + i0 * b) / (c * w + i0)
R4_attn = (a * c * w * L + i0 * b) / (a * c * w + i0)

# Plotting
fig, axs = plt.subplots(2, 2, figsize=(8, 8))
titles = ["Model 1: Full Attention (Additive b)",
          "Model 2: Stimulus Only Attention",
          "Model 3: Asymmetric Attention on b",
          "Model 4: Stimulus Attention Only (no attention on b)"]

model_pairs = [(R1_no_attn, R1_attn),
               (R2_no_attn, R2_attn),
               (R3_no_attn, R3_attn),
               (R4_no_attn, R4_attn)]

for ax, (r_no_attn, r_attn), title in zip(axs.flat, model_pairs, titles):
    ax.plot(c, r_no_attn, label="No Attention", color="blue")
    ax.plot(c, r_attn, label="Attention", color="green", linestyle="--")
    ax.set_xscale('log')
    ax.set_title(title)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Contrast")
    ax.set_ylabel("Response")
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()

######## model 3 with 2 stimuli (preferred and non-preferred)

# Parameters
c = np.logspace(np.log10(0.01), np.log10(10.0), 100)  # contrast values from 0.01 to 1 in log space
L1 = 1.0      # Preferred stimulus
L2 = 0.2      # Non-preferred stimulus
w1 = 1.0
w2 = 1.0
b = 0.1
i0 = 0.3
a1_no_attn = 1.0
a1_attn = 2.0  # Attention on stimulus 1 only

# Model 3: (c1*w1*L1 + c2*w2*L2 + i0*b) / (c1*w1 + c2*w2 + i0)

# Without attention
R_no_attn = (c * w1 * L1 + c * w2 * L2 + i0 * b) / (c * w1 + c * w2 + i0)

# With attention on preferred stimulus (stimulus 1)
R_attn = (a1_attn * c * w1 * L1 + c * w2 * L2 + i0 * b) / (a1_attn * c * w1 + c * w2 + i0)

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(c, R_no_attn, label="No Attention", color='blue')
plt.plot(c, R_attn, label="Attention on Preferred Stimulus", color='orange', linestyle='--')
plt.xlabel("Contrast")
plt.ylabel("Response")
plt.title("Model 3: Two Stimuli (Preferred & Non-Preferred) in RF")
plt.legend()
plt.grid(True)
plt.xscale("log")
plt.ylim(0, 1.1)
plt.tight_layout()
plt.show()


####################### 3 summary figures (contrast, activity gain and multiple stimuli with attention to preferred)


# === Combined Figure: Contrast Gain, Response Gain, Two Stimuli, and Asymmetric Attention ===
# Fixed contrast vector
c = np.logspace(np.log10(0.01), np.log10(100.0), 100)

# Parameters for each subplot
# Subplot 1 (Model 4: Contrast Gain)
params_cg = {
    "w": 1.0,
    "L": 1.0,
    "b": 0.1,
    "i0": 0.5,
    "attn_levels": [1, 2]
}

# Subplot 2 (Model 3: Response Gain)
params_rg = {
    "w": 1.0,
    "L": 1.0,
    "b": 0.1,
    "i0": 0.3,
    "attn_levels": [1, 3]
}

# Subplot 3 (Model 3: Two Stimuli)
params_dual = {
    "L1": 1.0,
    "L2": 0.2,
    "w1": 1.0,
    "w2": 1.0,
    "b": 0.1,
    "i0": 0.3,
    "a1_no_attn": 1.0,
    "a1_attn": 2.0
}

# Create combined figure with 4 subplots
fig, axs = plt.subplots(1, 4, figsize=(24, 5), sharey=True)

# Subplot 1: Contrast Gain (Model 4)
for a in params_cg["attn_levels"]:
    R = (a * c * params_cg["w"] * params_cg["L"] + params_cg["i0"] * params_cg["b"]) / (a * c * params_cg["w"] + params_cg["i0"])
    axs[0].plot(c, R, label=f'Attention = {a}')
axs[0].set_xscale('log')
axs[0].set_xlabel("Contrast")
axs[0].set_ylabel("Response")
axs[0].set_title("Model 4: Contrast Gain")
axs[0].set_ylim(0, 1.1)
axs[0].grid(True)
axs[0].legend()

# Subplot 2: Response Gain (Model 3)
for a in params_rg["attn_levels"]:
    R = (a * c * params_rg["w"] * params_rg["L"] + params_rg["i0"] * a * params_rg["b"]) / (a * c * params_rg["w"] + params_rg["i0"])
    axs[1].plot(c, R, label=f'Attention = {a}')
axs[1].set_xscale('log')
axs[1].set_xlabel("Contrast")
axs[1].set_title("Model 3: Activity/Response Gain")
axs[1].set_ylim(0, 1.1)
axs[1].grid(True)
axs[1].legend()

# Subplot 3: Two Stimuli in RF (Model 3)
R_no_attn = (c * params_dual["w1"] * params_dual["L1"] + c * params_dual["w2"] * params_dual["L2"] + params_dual["i0"] * params_dual["b"]) / (c * params_dual["w1"] + c * params_dual["w2"] + params_dual["i0"])
R_attn = (params_dual["a1_attn"] * c * params_dual["w1"] * params_dual["L1"] + c * params_dual["w2"] * params_dual["L2"] + params_dual["i0"] * params_dual["b"]) / (params_dual["a1_attn"] * c * params_dual["w1"] + c * params_dual["w2"] + params_dual["i0"])
axs[2].plot(c, R_no_attn, label="No Attention", color='blue')
axs[2].plot(c, R_attn, label="Attention on Preferred", color='orange', linestyle='--')
axs[2].set_xscale('log')
axs[2].set_xlabel("Contrast")
axs[2].set_title("Model 3: Two Stimuli, Attention on Preferred")
axs[2].set_ylim(0, 1.1)
axs[2].grid(True)
axs[2].legend()

# Subplot 4: Asymmetric Attention on Stimulus vs i0 (Model 2)
a_stim = 2.0
a_i0 = 1
a_baseline = 1.0

# Baseline (no attention)
R_asym_no_attn = (a_baseline * c * params_rg["w"] * params_rg["L"] + a_baseline * params_rg["i0"] * params_rg["b"]) / (a_baseline * c * params_rg["w"] + a_baseline * params_rg["i0"])
# Asymmetric attention: different gain on stimulus vs i0
R_asym_attn = (a_stim * c * params_rg["w"] * params_rg["L"] + a_i0 * params_rg["i0"] * params_rg["b"]) / (a_stim * c * params_rg["w"] + a_i0 * params_rg["i0"])

axs[3].plot(c, R_asym_no_attn, label="No Attention", color='blue')
axs[3].plot(c, R_asym_attn, label="Asymmetric Attention", color='orange', linestyle='--')
axs[3].set_xscale('log')
axs[3].set_xlabel("Contrast")
axs[3].set_title("Asymmetric Attention: Stim vs i0")
axs[3].set_ylim(0, 1.1)
axs[3].grid(True)
axs[3].legend()

plt.tight_layout()
plt.show()
