for filler in range(1):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    from scipy.io import savemat

    # ================= USER SETTINGS =================
    N_TRIALS       = 500      # trials per (pair × condition) per loop
    DT_SECONDS     = 0.25     # 250 ms counting window

    # ---- Weight gains per loop (global, shared across all neurons) ----
    WGAIN_SHAPE    = 20.0     # Gamma(shape, scale=1/shape), mean=1, var=1/shape
    APPLY_W0_GAIN  = True     # If False, w0 is not scaled (g0=1)

    # ---- Trialwise gain (creates r_sc; weights remain fixed) ----
    SIGMA_GAIN     = 1.00     # log-SD of per-trial gain
    RHO_SHARED     = 1.00     # shared fraction in gain latent
    KAPPA          = 0.0      # 1.0 cancels gain in denom; <1 lets gain modulate rates

    # ---- BC gating and sign forcing ----
    BC_THRESH      = 0.65      # BC split (>=thr => positive shared sign; <thr => negative)
    FORCE_SIGN_WITH_BC = True  # enforce sign of r_sc by BC

    # ---- Reps / seeds ----
    N_REPS         = 200     # number of loops (each draws new g0..g4)
    RNG_SEED       = 123

    # ========== PREREQS you must already have in memory ==========
    # all_combs: list of (unit_i, unit_j)
    # masterGoodUnits: np.array of unit IDs (strings), aligned with rfWeightNeuronFits
    # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
    # bc_lookup[tuple(sorted((u1,u2)))] = BC float


    # ---------- paired-condition decoder (your indices 10..17, 19..26) ----------
    def decode_paired_condition(i):
        if 10 <= i <= 17:
            site = (i - 10) // 2
            periph_is_pref = (i - 10) % 2 == 1
            center_is_pref = False
        elif 19 <= i <= 26:
            site = (i - 19) // 2
            periph_is_pref = (i - 19) % 2 == 1
            center_is_pref = True
        else:
            raise ValueError(f"Condition {i} not in paired blocks (10..17 or 19..26).")
        return center_is_pref, periph_is_pref, site


    def sign_from_bc(bc, thresh=BC_THRESH):
        return +1.0 if (np.isfinite(bc) and bc >= thresh) else -1.0


    # ---------- Trialwise gain (creates r_sc without touching weights) ----------
    def draw_gain_latents(T, rho, sign, rng):
        """
        Draw one set of trialwise gain latents for a pair×condition.
        We'll reuse these SAME latents for baseline vs jittered.
        """
        eta1 = rng.normal(0, 1, size=T)
        eta2 = rng.normal(0, 1, size=T)
        xi   = rng.normal(0, 1, size=T)
        g1 = np.sqrt(1 - rho) * eta1 + np.sqrt(rho) * xi
        g2 = np.sqrt(1 - rho) * eta2 + sign * np.sqrt(rho) * xi
        return g1, g2


    def apply_gain_to_rate(base_rate, g, sigma_gain, kappa):
        # Multiplicative gain on rate; if kappa<1, rate varies with g; if kappa=1, no effect
        return base_rate * np.exp((1.0 - kappa) * sigma_gain * g)


    # ---------- Base rate & F from frozen weights ----------
    def rate_from_frozen_weights(Lc, Lpnp, w0_frozen, wsite_frozen, sigma, base):
        num = w0_frozen * Lc + wsite_frozen * Lpnp
        den = w0_frozen + wsite_frozen + max(1e-9, float(sigma))
        return num / den + base


    def weight_fraction_from_frozen(w0_frozen, wsite_frozen, sigma):
        den = w0_frozen + wsite_frozen + max(1e-9, float(sigma))
        f = (w0_frozen + wsite_frozen) / den
        return float(np.clip(f, 0.0, 1.0))


    # ---------- faster correlation ----------
    def fast_corr(x, y):
        sx = x.std(); sy = y.std()
        if sx == 0 or sy == 0:
            return np.nan
        return float(np.corrcoef(x, y)[0,1])


    # ---------- Precompute static per (pair, cond) info ----------
    def precompute_pair_cond_table(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup):
        params_by_unit = {unit: rfWeightNeuronFits[i] for i, unit in enumerate(masterGoodUnits)}
        paired_indices = list(range(10, 18)) + list(range(19, 27))
        rows = []
        for (u1, u2) in all_combs:
            if (u1 not in params_by_unit) or (u2 not in params_by_unit):
                continue
            p1 = params_by_unit[u1]; p2 = params_by_unit[u2]
            Lp1, Lnp1, w1_1, w2_1, w3_1, w4_1, sig1, base1 = p1
            Lp2, Lnp2, w1_2, w2_2, w3_2, w4_2, sig2, base2 = p2
            key = tuple(sorted((u1, u2)))
            bc  = bc_lookup.get(key, np.nan)
            sgn = sign_from_bc(bc)
            for cond_idx in paired_indices:
                center_is_pref, periph_is_pref, site = decode_paired_condition(cond_idx)
                wsite1 = [w1_1, w2_1, w3_1, w4_1][site]
                wsite2 = [w1_2, w2_2, w3_2, w4_2][site]
                Lc1   = Lp1 if center_is_pref else Lnp1
                Lpnp1 = Lp1 if periph_is_pref else Lnp1
                Lc2   = Lp2 if center_is_pref else Lnp2
                Lpnp2 = Lp2 if periph_is_pref else Lnp2
                rows.append((u1,u2,cond_idx,site,center_is_pref,periph_is_pref,
                             Lc1,Lpnp1,Lc2,Lpnp2,wsite1,wsite2,sig1,base1,sig2,base2,bc,sgn))
        cols = ["unit_i","unit_j","cond_index","site","center_is_pref","periph_is_pref",
                "Lc1","Lpnp1","Lc2","Lpnp2","wsite1","wsite2","sig1","base1","sig2","base2","bhattCoef","sign"]
        return pd.DataFrame(rows, columns=cols), params_by_unit


    # ---------- Per-loop weight gains (Gamma, mean=1) ----------
    def draw_weight_gains(rng, shape=WGAIN_SHAPE, apply_w0=APPLY_W0_GAIN):
        g0 = rng.gamma(shape, 1.0/shape) if apply_w0 else 1.0
        g1 = rng.gamma(shape, 1.0/shape)
        g2 = rng.gamma(shape, 1.0/shape)
        g3 = rng.gamma(shape, 1.0/shape)
        g4 = rng.gamma(shape, 1.0/shape)
        return g0, g1, g2, g3, g4


    def loop_noise_indices_from_gains(g0, g1, g2, g3, g4):
        lg = np.log([g0, g1, g2, g3, g4])        # log-gains; ~symmetric around 0 (mean≈0)
        noise_L2 = float(np.sqrt(np.mean(lg**2)))
        noise_L1 = float(np.mean(np.abs(lg)))
        # effective F jitter per site in paired conditions uses w0 + w_site
        lg_eff = lg[0] + lg[1:]                  # size 4: sites 1..4
        noise_eff_std  = float(np.std(lg_eff))
        noise_eff_mean = float(np.mean(lg_eff))
        signed_rms = float(np.sign(noise_eff_mean) * np.sqrt(np.mean(lg_eff ** 2)))
        return dict(noise_L2=noise_L2,
                    noise_L1=noise_L1,
                    noise_eff_std=noise_eff_std,
                    noise_eff_mean=noise_eff_mean,
                    signed_rms=signed_rms)


    # ---------- Simulate one rep (Baseline vs Jittered using SAME gain latents) ----------
    def simulate_one_rep(pre_df, params_by_unit, rng):
        # Global gains for this loop
        g0, g1, g2, g3, g4 = draw_weight_gains(rng, WGAIN_SHAPE, APPLY_W0_GAIN)
        lg0, lg1, lg2, lg3, lg4 = np.log([g0, g1, g2, g3, g4])
        noise_idx = loop_noise_indices_from_gains(g0, g1, g2, g3, g4)

        out = {
            "unit_i": [], "unit_j": [], "cond_index": [],
            "site": [], "center_is_pref": [], "periph_is_pref": [],
            "bhattCoef": [], "sign": [],
            # F and corr (baseline vs jittered) + deltas
            "F_base": [], "F_jitter": [], "F_delta": [],
            "corr_base": [], "corr_jitter": [], "corr_delta": [],
            # cov/var (baseline and jittered)
            "mean1_base": [], "mean2_base": [], "var1_base": [], "var2_base": [], "cov12_base": [],
            "mean1_jit":  [], "mean2_jit":  [], "var1_jit":  [], "var2_jit":  [], "cov12_jit":  [],
            # effective log gain for this site (paired w0+w_site)
            "lg_eff_row": [],
            # store the 5 gains (same for whole rep)
            "g0": [], "g1": [], "g2": [], "g3": [], "g4": [],
        }

        w0_base = 1.0

        for row in pre_df.itertuples(index=False):
            # one set of trialwise gain latents, reused for baseline AND jittered
            g1_lat, g2_lat = draw_gain_latents(N_TRIALS, RHO_SHARED, row.sign, rng)

            # ----- Baseline weights (all gains = 1) -----
            w0_b   = 1.0
            wsite1_b = row.wsite1
            wsite2_b = row.wsite2
            r1_base_b = rate_from_frozen_weights(row.Lc1, row.Lpnp1, w0_b, wsite1_b, row.sig1, row.base1)
            r2_base_b = rate_from_frozen_weights(row.Lc2, row.Lpnp2, w0_b, wsite2_b, row.sig2, row.base2)
            r1_t_b = apply_gain_to_rate(r1_base_b, g1_lat, SIGMA_GAIN, KAPPA)
            r2_t_b = apply_gain_to_rate(r2_base_b, g2_lat, SIGMA_GAIN, KAPPA)
            lam1_b = np.clip(r1_t_b * DT_SECONDS, 0, None)
            lam2_b = np.clip(r2_t_b * DT_SECONDS, 0, None)
            k1_b = rng.poisson(lam1_b)
            k2_b = rng.poisson(lam2_b)
            mean1_b = float(np.mean(k1_b))
            mean2_b = float(np.mean(k2_b))
            var1_b = float(np.var(k1_b, ddof=1))
            var2_b = float(np.var(k2_b, ddof=1))
            cov12_b = float(np.cov(k1_b, k2_b, ddof=1)[0, 1])
            r_sc_b = fast_corr(k1_b, k2_b)
            if FORCE_SIGN_WITH_BC and np.isfinite(r_sc_b):
                r_sc_b = float(row.sign) * abs(r_sc_b)
            f1_b = weight_fraction_from_frozen(w0_b, wsite1_b, row.sig1)
            f2_b = weight_fraction_from_frozen(w0_b, wsite2_b, row.sig2)
            F_geo_b = float(np.sqrt(np.clip(f1_b * f2_b, 0.0, 1.0)))

            # ----- Jittered weights (apply global gamma gains) -----
            w0_j = w0_base * g0
            wsite1_j = row.wsite1 * [g1, g2, g3, g4][row.site]
            wsite2_j = row.wsite2 * [g1, g2, g3, g4][row.site]
            r1_base_j = rate_from_frozen_weights(row.Lc1, row.Lpnp1, w0_j, wsite1_j, row.sig1, row.base1)
            r2_base_j = rate_from_frozen_weights(row.Lc2, row.Lpnp2, w0_j, wsite2_j, row.sig2, row.base2)
            r1_t_j = apply_gain_to_rate(r1_base_j, g1_lat, SIGMA_GAIN, KAPPA)
            r2_t_j = apply_gain_to_rate(r2_base_j, g2_lat, SIGMA_GAIN, KAPPA)
            lam1_j = np.clip(r1_t_j * DT_SECONDS, 0, None)
            lam2_j = np.clip(r2_t_j * DT_SECONDS, 0, None)
            k1_j = rng.poisson(lam1_j)
            k2_j = rng.poisson(lam2_j)
            mean1_j = float(np.mean(k1_j))
            mean2_j = float(np.mean(k2_j))
            var1_j = float(np.var(k1_j, ddof=1))
            var2_j = float(np.var(k2_j, ddof=1))
            cov12_j = float(np.cov(k1_j, k2_j, ddof=1)[0, 1])
            r_sc_j = fast_corr(k1_j, k2_j)
            if FORCE_SIGN_WITH_BC and np.isfinite(r_sc_j):
                r_sc_j = float(row.sign) * abs(r_sc_j)
            f1_j = weight_fraction_from_frozen(w0_j, wsite1_j, row.sig1)
            f2_j = weight_fraction_from_frozen(w0_j, wsite2_j, row.sig2)
            F_geo_j = float(np.sqrt(np.clip(f1_j * f2_j, 0.0, 1.0)))

            # ----- deltas (pure effect of weight noise, given identical gain latents)
            delta_rsc = r_sc_j - r_sc_b
            delta_F   = F_geo_j - F_geo_b

            # effective log gain for this site's F (log g0 + log g_site)
            lg_site   = [lg1, lg2, lg3, lg4][row.site]
            lg_eff_row = lg0 + lg_site

            # append row
            out["unit_i"].append(row.unit_i); out["unit_j"].append(row.unit_j)
            out["cond_index"].append(row.cond_index); out["site"].append(row.site)
            out["center_is_pref"].append(bool(row.center_is_pref))
            out["periph_is_pref"].append(bool(row.periph_is_pref))
            out["bhattCoef"].append(float(row.bhattCoef)); out["sign"].append(float(row.sign))

            out["F_base"].append(F_geo_b)
            out["F_jitter"].append(F_geo_j)
            out["F_delta"].append(delta_F)
            out["corr_base"].append(r_sc_b)
            out["corr_jitter"].append(r_sc_j)
            out["corr_delta"].append(delta_rsc)

            out["mean1_base"].append(mean1_b)
            out["mean2_base"].append(mean2_b)
            out["var1_base"].append(var1_b)
            out["var2_base"].append(var2_b)
            out["cov12_base"].append(cov12_b)

            out["mean1_jit"].append(mean1_j)
            out["mean2_jit"].append(mean2_j)
            out["var1_jit"].append(var1_j)
            out["var2_jit"].append(var2_j)
            out["cov12_jit"].append(cov12_j)

            out["lg_eff_row"].append(float(lg_eff_row))
            out["g0"].append(g0)
            out["g1"].append(g1)
            out["g2"].append(g2)
            out["g3"].append(g3)
            out["g4"].append(g4)

        df_rep = pd.DataFrame(out)

        # realized spreads (diagnostics)
        F_std_rep   = float(np.nanstd(df_rep["F_jitter"].to_numpy()))
        Fd_std_rep  = float(np.nanstd(df_rep["F_delta"].to_numpy()))
        Cd_std_rep  = float(np.nanstd(df_rep["corr_delta"].to_numpy()))

        # per-loop noise summary
        loop_summary = dict(F_std_rep=F_std_rep, Fdelta_std_rep=Fd_std_rep, Cdelta_std_rep=Cd_std_rep,
                            g0=g0, g1=g1, g2=g2, g3=g3, g4=g4, **noise_idx)
        return df_rep, loop_summary


    # ---------- Per-rep slope diagnostics ----------
    def slope_diagnostics(df, bc_thresh=BC_THRESH):
        """
        Returns slopes for:
          - jittered:   correlation_jitter ~ F_jitter
          - delta:      corr_delta ~ F_delta
        computed separately for high-BC and low-BC subsets.
        """
        out = {}
        def slopes(x, y):
            x = np.asarray(x); y = np.asarray(y)
            mask = np.isfinite(x) & np.isfinite(y)
            if mask.sum() < 5 or np.std(x[mask]) == 0:
                return np.nan, np.nan
            rho, _ = spearmanr(x[mask], y[mask])
            beta = np.polyfit(x[mask], y[mask], 1)[0]
            return float(rho), float(beta)

        # High BC
        dhi = df[df["bhattCoef"] >= bc_thresh]
        rho_hi_j, beta_hi_j = slopes(dhi["F_jitter"],  dhi["corr_jitter"])
        rho_hi_d, beta_hi_d = slopes(dhi["F_delta"],   dhi["corr_delta"])

        # Low BC
        dlo = df[df["bhattCoef"] <  bc_thresh]
        rho_lo_j, beta_lo_j = slopes(dlo["F_jitter"],  dlo["corr_jitter"])
        rho_lo_d, beta_lo_d = slopes(dlo["F_delta"],   dlo["corr_delta"])

        out.update(dict(
            rho_highBC_jitter=rho_hi_j, beta_highBC_jitter=beta_hi_j,
            rho_lowBC_jitter=rho_lo_j,  beta_lowBC_jitter=beta_lo_j,
            rho_highBC_delta=rho_hi_d,  beta_highBC_delta=beta_hi_d,
            rho_lowBC_delta=rho_lo_d,   beta_lowBC_delta=beta_lo_d
        ))
        return out


    # ---------- Driver ----------
    def run_simulation_weight_gains_delta(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
                                          n_reps=N_REPS, rng_seed=RNG_SEED):
        pre_df, params_by_unit = precompute_pair_cond_table(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup)
        rng_master = np.random.default_rng(rng_seed)

        long_rows, summaries = [], []
        for rep in range(n_reps):
            print(rep)
            rng = np.random.default_rng(rng_master.integers(1, 2**31-1))
            df_rep, loop_idx = simulate_one_rep(pre_df, params_by_unit, rng)
            df_rep["rep"] = rep
            long_rows.append(df_rep)

            diag = slope_diagnostics(df_rep, bc_thresh=BC_THRESH)
            loop_idx.update(diag)
            loop_idx["rep"] = rep
            summaries.append(loop_idx)

        full_df = pd.concat(long_rows, ignore_index=True)
        summary_df = pd.DataFrame(summaries)
        return full_df, summary_df


    # ===================== RUN =====================
    full_df, summary_df = run_simulation_weight_gains_delta(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_reps=N_REPS, rng_seed=RNG_SEED
    )

    print("Full DF shape:", full_df.shape)
    print("Summary (cols):", list(summary_df.columns))
    print(summary_df[["rep", "noise_L2", "noise_eff_mean", "noise_eff_std",
                      "F_std_rep", "Fdelta_std_rep", "Cdelta_std_rep",
                      "beta_highBC_jitter", "beta_lowBC_jitter", "beta_highBC_delta",
                      "beta_lowBC_delta"]].head())

    # ===================== SAVE =====================
    full_df.to_csv("paired_weightgains_long_delta.csv", index=False)
    summary_df.to_csv("paired_weightgains_summary_delta.csv", index=False)

    # Save MAT (numeric columns only)
    full_numeric = {k: full_df[k].to_numpy() for k in full_df.columns if np.issubdtype(full_df[k].dtype, np.number)}
    sum_numeric  = {k: summary_df[k].to_numpy() for k in summary_df.columns if np.issubdtype(summary_df[k].dtype, np.number)}
    savemat("paired_weightgains_long_delta.mat", full_numeric)
    savemat("paired_weightgains_summary_delta.mat", sum_numeric)
    print("Saved CSV & MAT files.")

    # ===================== QUICK PLOTS =====================

    # (1) Δ correlation vs Δ F (one rep, sanity check)
    rep_to_plot = 0
    rep_df = full_df[full_df["rep"] == rep_to_plot]
    plt.figure(figsize=(6.8,5.0))
    plt.scatter(rep_df["F_delta"], rep_df["corr_delta"], c=rep_df["bhattCoef"],
                s=12, cmap="viridis")
    plt.axhline(0, color="k", lw=1); plt.axvline(0, color="k", lw=1)
    plt.xlabel("ΔF  (F_jitter − F_base)")
    plt.ylabel("Δ r_sc  (jitter − base)")
    plt.title(f"Within-loop Δ effect (rep {rep_to_plot})")
    plt.colorbar(label="BC")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    # (2) Per-loop slopes vs loop noise (use whichever metric you like)
    NOISE_COL = "noise_L2"  # alternatives: "noise_eff_std", "noise_eff_mean"
    plt.figure(figsize=(8,5))
    plt.scatter(summary_df[NOISE_COL], summary_df["beta_highBC_delta"], alpha=0.7, label=f"Δ slope (High-BC ≥ {BC_THRESH})")
    plt.scatter(summary_df[NOISE_COL], summary_df["beta_lowBC_delta"],  alpha=0.7, label=f"Δ slope (Low-BC < {BC_THRESH})")
    # trend lines
    for col in ["beta_highBC_delta","beta_lowBC_delta"]:
        m = summary_df[[NOISE_COL, col]].dropna()
        if len(m) > 2:
            x = m[NOISE_COL].to_numpy(); y = m[col].to_numpy()
            b, a = np.polyfit(x, y, 1)  # y ≈ a + b*x
            xs = np.linspace(x.min(), x.max(), 100)
            plt.plot(xs, a + b*xs, lw=2)
    plt.axhline(0, color="k", lw=1)
    plt.xlabel(f"Loop noise level ({NOISE_COL})")
    plt.ylabel("OLS slope:  Δ r_sc  ~  Δ F")
    plt.title("Per-loop Δ slope vs weight-noise")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    # (3) Jittered scatter F vs BC colored by r_sc (one rep)
    plt.figure(figsize=(6.8, 5.2))
    sc = plt.scatter(rep_df["F_jitter"], rep_df["bhattCoef"],
                     c=rep_df["corr_delta"], s=10, cmap="coolwarm", vmin=-0.3, vmax=0.3)
    plt.colorbar(sc, label="r_sc (jittered)")
    plt.xlabel("Geometric mean F (jittered)")
    plt.ylabel("Bhattacharyya coefficient (BC)")
    plt.title(f"Corr vs F × BC (jittered, rep {rep_to_plot})")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

#working
for filler in range(1):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.io import savemat


    # ---------- helpers ----------
    def decode_paired_condition(i):
        """Return center_is_pref (bool), periph_is_pref (bool), site_idx (0..3)."""
        if 10 <= i <= 17:
            site = (i - 10) // 2
            periph_is_pref = (i - 10) % 2 == 1
            center_is_pref = False
        elif 19 <= i <= 26:
            site = (i - 19) // 2
            periph_is_pref = (i - 19) % 2 == 1
            center_is_pref = True
        else:
            raise ValueError(f"Condition {i} is not paired (10..17 or 19..26).")
        return center_is_pref, periph_is_pref, site


    def sign_from_bc(bc, thresh=BC_THRESH):
        return +1.0 if (np.isfinite(bc) and bc >= thresh) else -1.0


    def rate_from_weights(Lc, Lpnp, w0, ws, sigma, base):
        num = w0 * Lc + ws * Lpnp
        den = w0 + ws + np.maximum(1e-9, sigma)
        return num / den + base


    def weight_fraction(w0, ws, sigma):
        den = w0 + ws + np.maximum(1e-9, sigma)
        return np.clip((w0 + ws) / den, 0.0, 1.0)


    # Online (vectorized) update for means, variances, covariance
    def online_update(x, y, n, mean1, mean2, M2_1, M2_2, C12):
        """
        x,y are vectors (one sample per row)
        Updates Welford accumulators.
        """
        n_new = n + 1.0
        dx = x - mean1
        mean1_new = mean1 + dx / n_new
        dy = y - mean2
        mean2_new = mean2 + dy / n_new
        M2_1_new = M2_1 + dx * (x - mean1_new)
        M2_2_new = M2_2 + dy * (y - mean2_new)
        C12_new  = C12  + dx * (y - mean2_new)
        return n_new, mean1_new, mean2_new, M2_1_new, M2_2_new, C12_new


    # ---------- Precompute static per (pair, cond) table ----------
    def precompute_table(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup):
        params = {unit: rfWeightNeuronFits[i] for i, unit in enumerate(masterGoodUnits)}
        paired_indices = list(range(10, 18)) + list(range(19, 27))
        rows = []
        for (u1, u2) in all_combs:
            if (u1 not in params) or (u2 not in params):
                continue
            p1 = params[u1]; p2 = params[u2]
            Lp1, Lnp1, w1_1, w2_1, w3_1, w4_1, sig1, base1 = p1
            Lp2, Lnp2, w1_2, w2_2, w3_2, w4_2, sig2, base2 = p2
            key = tuple(sorted((u1, u2)))
            bc  = float(bc_lookup.get(key, np.nan))
            sgn = sign_from_bc(bc)

            for cond_idx in paired_indices:
                center_is_pref, periph_is_pref, site = decode_paired_condition(cond_idx)
                ws1 = [w1_1, w2_1, w3_1, w4_1][site]
                ws2 = [w1_2, w2_2, w3_2, w4_2][site]
                Lc1   = Lp1 if center_is_pref  else Lnp1
                Lpnp1 = Lp1 if periph_is_pref  else Lnp1
                Lc2   = Lp2 if center_is_pref  else Lnp2
                Lpnp2 = Lp2 if periph_is_pref  else Lnp2

                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1, Lpnp1, Lc2, Lpnp2, ws1, ws2, sig1, base1, sig2, base2, bc, sgn))

        cols = ["unit_i", "unit_j", "cond_index", "site", "center_is_pref",
                "periph_is_pref", "Lc1", "Lpnp1", "Lc2", "Lpnp2", "wsite1",
                "wsite2", "sig1", "base1", "sig2", "base2", "bhattCoef", "sign"]
        pre = pd.DataFrame(rows, columns=cols)

        # Baseline F (no weight noise), per row
        w0 = 1.0
        f1_base = weight_fraction(w0, pre["wsite1"].to_numpy(), pre["sig1"].to_numpy())
        f2_base = weight_fraction(w0, pre["wsite2"].to_numpy(), pre["sig2"].to_numpy())
        F_geo_base = np.sqrt(f1_base * f2_base)
        pre["F_base"] = F_geo_base
        return pre


    # ---------- Main simulation (trial-by-trial weight gains) ----------
    def simulate_trialwise_weight_noise(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    ):
        rng = np.random.default_rng(rng_seed)
        pre = precompute_table(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup)

        # 1) Does F correlate with baseline rates in your fitted set?
        r1_nom = rate_from_weights(pre["Lc1"], pre["Lpnp1"], 1.0, pre["wsite1"], pre["sig1"], pre["base1"])
        r2_nom = rate_from_weights(pre["Lc2"], pre["Lpnp2"], 1.0, pre["wsite2"], pre["sig2"], pre["base2"])
        rate_mean = 0.5 * (r1_nom + r2_nom)
        print("corr(F_base, mean rate) =", np.corrcoef(pre["F_base"], rate_mean)[0, 1])

        # 2) Does r_sc grow with the mean rate (as the theory says)?
        x = full_df["r_sc"].to_numpy()
        y = 0.5 * (full_df["mean1"] + full_df["mean2"])
        m = np.isfinite(x) & np.isfinite(y)
        print("corr(r_sc, mean rate) =", np.corrcoef(x[m], y[m])[0, 1])

        # Also check the two blocks separately:
        mask_periph_pref = pre["periph_is_pref"].astype(bool)
        print("periph-pref block corr =", np.corrcoef(pre.loc[mask_periph_pref, "F_base"],
                                                      rate_mean[mask_periph_pref])[0, 1])
        mask_center_pref = pre["center_is_pref"].astype(bool)
        print("center-pref block corr =", np.corrcoef(pre.loc[mask_center_pref, "F_base"],
                                                      rate_mean[mask_center_pref])[0, 1])

        # Pull out as numpy arrays for speed
        Lc1   = pre["Lc1"].to_numpy()
        Lpnp1 = pre["Lpnp1"].to_numpy()
        Lc2   = pre["Lc2"].to_numpy()
        Lpnp2 = pre["Lpnp2"].to_numpy()
        ws1   = pre["wsite1"].to_numpy()
        ws2   = pre["wsite2"].to_numpy()
        sig1  = pre["sig1"].to_numpy()
        sig2  = pre["sig2"].to_numpy()
        b1    = pre["base1"].to_numpy()
        b2    = pre["base2"].to_numpy()
        site  = pre["site"].to_numpy().astype(int)      # 0..3
        sgn   = pre["sign"].to_numpy().astype(float)    # +1 / -1
        bc    = pre["bhattCoef"].to_numpy()

        n_rows = len(pre)
        w0_base = 1.0

        # ---- per-row shared fraction for the latent, scaled by BC magnitude ----
        if USE_BC_SCALED_RHO:
            rho_eff = RHO_MIN + (RHO_MAX - RHO_MIN) * np.clip(bc, 0.0, 1.0)
        else:
            rho_eff = np.full(n_rows, RHO_SHARED, dtype=float)

        # Online accumulators per row
        n     = np.zeros(n_rows, dtype=float)
        m1    = np.zeros(n_rows, dtype=float)
        m2    = np.zeros(n_rows, dtype=float)
        M2_1  = np.zeros(n_rows, dtype=float)
        M2_2  = np.zeros(n_rows, dtype=float)
        C12   = np.zeros(n_rows, dtype=float)

        # Trial loop
        for t in range(n_trials):
            # 5 global gains (shared by all rows/pairs) for this trial
            g0 = rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) if APPLY_W0_GAIN else 1.0
            g  = np.array([rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) for _ in range(4)], dtype=float)  # g1..g4

            # Select site gain per row
            g_site = g[site]  # shape (n_rows,)

            # Jittered weights for this trial
            w0_j   = w0_base * g0                     # same for all rows
            ws1_j  = ws1 * g_site
            ws2_j  = ws2 * g_site

            # Base rates (given this trial's weight gains)
            r1_base = rate_from_weights(Lc1, Lpnp1, w0_j, ws1_j, sig1, b1)
            r2_base = rate_from_weights(Lc2, Lpnp2, w0_j, ws2_j, sig2, b2)

            # Trialwise latent (BC-scaled shared fraction, BC-signed)
            eta1 = rng.normal(size=n_rows)
            eta2 = rng.normal(size=n_rows)
            xi   = rng.normal(size=n_rows)           # shared per row
            s1 = np.sqrt(1.0 - rho_eff)
            s2 = np.sqrt(rho_eff)

            lin1 = s1 * eta1 + s2 * xi
            lin2 = s1 * eta2 + sgn * s2 * xi

            gain1 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin1)
            gain2 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin2)

            # Trial rates with latent gain
            r1_t = r1_base * gain1
            r2_t = r2_base * gain2

            # Poisson counts for this trial
            lam1 = np.clip(r1_t * DT_SECONDS, 0, None)
            lam2 = np.clip(r2_t * DT_SECONDS, 0, None)
            k1 = rng.poisson(lam1)
            k2 = rng.poisson(lam2)

            # Online update
            n, m1, m2, M2_1, M2_2, C12 = online_update(k1, k2, n, m1, m2, M2_1, M2_2, C12)

        # Final stats
        dof = np.maximum(1.0, n - 1.0)
        var1 = M2_1 / dof
        var2 = M2_2 / dof
        cov12 = C12 / dof
        denom = np.sqrt(np.clip(var1 * var2, 1e-12, None))
        r_sc = np.where(denom > 0, cov12 / denom, np.nan)

        # Pack result DataFrame — one row per (pair × condition)
        out = pre.copy()
        out["mean1"] = m1
        out["mean2"] = m2
        out["var1"]  = var1
        out["var2"]  = var2
        out["cov12"] = cov12
        out["r_sc"]  = r_sc
        out["trials"] = n.astype(int)
        # Optionally store the per-row rho used
        out["rho_eff"] = rho_eff

        return out

    # ================= USER SETTINGS =================
    N_TRIALS      = 10000   # trials per (pair × condition)
    DT_SECONDS    = 0.25    # 250 ms counting window

    # ---- Per-trial weight gains (global across all neurons & pairs) ----
    WGAIN_SHAPE   = 10.0    # Gamma(shape, scale=1/shape); mean=1, var=1/shape
    APPLY_W0_GAIN = True    # if False, g0=1 always

    # ---- Trialwise latent gain (imposes r_sc & sign by BC) ----
    SIGMA_GAIN    = 0.4     # log-SD of per-trial rate gain
    RHO_SHARED    = 0.4     # used only if USE_BC_SCALED_RHO=False
    KAPPA         = 0.3     # 0 => gain fully modulates rate; 1 => cancels in denom

    # ---- Make the shared fraction depend on BC magnitude ----
    USE_BC_SCALED_RHO = True   # True => use BC to scale shared fraction
    RHO_MIN = 0.20             # shared fraction when BC = 0
    RHO_MAX = 0.80             # shared fraction when BC = 1

    # ---- BC sign gate ----
    BC_THRESH     = 0.65    # ≥thr => positive shared sign; <thr => negative

    # ---- RNG ----
    RNG_SEED      = 123

    # ========== PREREQS you already have in memory ==========
    # all_combs: list of (unit_i, unit_j)
    # masterGoodUnits: np.array of unit IDs (strings), aligned with rfWeightNeuronFits
    # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
    # bc_lookup[tuple(sorted((u1,u2)))] = BC float

    # ===================== RUN =====================
    full_df = simulate_trialwise_weight_noise(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    )
    print(full_df.head())

    # Save numeric outputs to .mat for MATLAB mesh plotting (F_base × BC → r_sc)
    mat_path = "trialwise_weightnoise_results.mat"
    numeric_cols = {k: full_df[k].to_numpy() for k in full_df.columns if np.issubdtype(full_df[k].dtype, np.number)}
    savemat(mat_path, numeric_cols)

    # Quick sanity plot
    vlim=0.1
    x = full_df["F_base"].to_numpy()
    y = full_df["bhattCoef"].to_numpy()
    z = full_df["r_sc"].to_numpy()
    msk = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    plt.figure(figsize=(6.6, 5))
    sc = plt.scatter(x[msk], y[msk], c=z[msk], s=10, cmap="coolwarm", vmin=-vlim, vmax=+vlim)
    plt.colorbar(sc, label="r_sc")
    plt.xlabel("Baseline F (geom. mean weight fraction)")
    plt.ylabel("Bhattacharyya coefficient (BC)")
    plt.title("r_sc over F × BC (trialwise weight noise, BC-scaled shared fraction)")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

# working version for 10000 trials draw
for filler in range(1):
    import numpy as np
    import pandas as pd
    from scipy.io import savemat
    import matplotlib.pyplot as plt


    # ---------- helpers ----------
    def decode_paired_condition(i):
        """Return center_is_pref (bool), periph_is_pref (bool), site_idx (0..3)."""
        if 10 <= i <= 17:
            site = (i - 10) // 2
            periph_is_pref = (i - 10) % 2 == 1
            center_is_pref = False
        elif 19 <= i <= 26:
            site = (i - 19) // 2
            periph_is_pref = (i - 19) % 2 == 1
            center_is_pref = True
        else:
            raise ValueError(f"Condition {i} is not paired (10..17 or 19..26).")
        return center_is_pref, periph_is_pref, site


    def sign_from_bc(bc, thresh=BC_THRESH):
        return +1.0 if (np.isfinite(bc) and bc >= thresh) else -1.0


    def rate_from_weights(Lc, Lpnp, w0, ws, sigma, base):
        num = w0 * Lc + ws * Lpnp
        den = w0 + ws + np.maximum(1e-9, sigma)
        return num / den + base


    def weight_fraction(w0, ws, sigma):
        den = w0 + ws + np.maximum(1e-9, sigma)
        return np.clip((w0 + ws) / den, 0.0, 1.0)


    def fast_corr(x, y):
        sx = x.std(); sy = y.std()
        if sx == 0 or sy == 0:
            return np.nan
        return float(np.corrcoef(x, y)[0, 1])


    # Online (vectorized) update for means, variances, covariance
    def online_update(x, y, n, mean1, mean2, M2_1, M2_2, C12):
        """
        x,y are vectors (one sample per row)
        Updates Welford accumulators.
        """
        n_new = n + 1.0
        dx = x - mean1
        mean1_new = mean1 + dx / n_new
        dy = y - mean2
        mean2_new = mean2 + dy / n_new
        M2_1_new = M2_1 + dx * (x - mean1_new)
        M2_2_new = M2_2 + dy * (y - mean2_new)
        C12_new  = C12  + dx * (y - mean2_new)
        return n_new, mean1_new, mean2_new, M2_1_new, M2_2_new, C12_new


    # ---------- Precompute table (paired + derived singles) ----------
    def precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup):
        """
        Build a long design matrix with 3 rows per 'paired' condition:
          stim_kind = 0: paired (center + that periphery site shown together)
          stim_kind = 1: single_center (center shown, periphery blank -> Lpnp=0)
          stim_kind = 2: single_periph (that periphery site shown, center blank -> Lc=0)
        """
        params = {unit: rfWeightNeuronFits[i] for i, unit in enumerate(masterGoodUnits)}
        paired_indices = list(range(10, 18)) + list(range(19, 27))
        rows = []

        for (u1, u2) in all_combs:
            if (u1 not in params) or (u2 not in params):
                continue
            p1 = params[u1]; p2 = params[u2]
            Lp1, Lnp1, w1_1, w2_1, w3_1, w4_1, sig1, base1 = p1
            Lp2, Lnp2, w1_2, w2_2, w3_2, w4_2, sig2, base2 = p2
            key = tuple(sorted((u1, u2)))
            bc  = float(bc_lookup.get(key, np.nan))
            sgn = sign_from_bc(bc)

            for cond_idx in paired_indices:
                center_is_pref, periph_is_pref, site = decode_paired_condition(cond_idx)
                ws1 = [w1_1, w2_1, w3_1, w4_1][site]
                ws2 = [w1_2, w2_2, w3_2, w4_2][site]
                # L's for the paired condition
                Lc1_p   = Lp1 if center_is_pref  else Lnp1
                Lpnp1_p = Lp1 if periph_is_pref  else Lnp1
                Lc2_p   = Lp2 if center_is_pref  else Lnp2
                Lpnp2_p = Lp2 if periph_is_pref  else Lnp2

                # ---- (0) paired row ----
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, Lpnp1_p, Lc2_p, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 0))

                # ---- (1) single_center row (periphery blank -> Lpnp=0) ----
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, 0.0, Lc2_p, 0.0, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 1))

                # ---- (2) single_periph row (center blank -> Lc=0) ----
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             0.0, Lpnp1_p, 0.0, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 2))

        cols = ["unit_i","unit_j","cond_index","site","center_is_pref","periph_is_pref",
                "Lc1","Lpnp1","Lc2","Lpnp2","wsite1","wsite2","sig1","base1","sig2","base2",
                "bhattCoef","sign","stim_kind"]
        pre = pd.DataFrame(rows, columns=cols)

        # Baseline F (no weight noise), per row (same across the 3 rows of a trio)
        w0 = 1.0
        f1_base = weight_fraction(w0, pre["wsite1"].to_numpy(), pre["sig1"].to_numpy())
        f2_base = weight_fraction(w0, pre["wsite2"].to_numpy(), pre["sig2"].to_numpy())
        F_geo_base = np.sqrt(f1_base * f2_base)
        pre["F_base"] = F_geo_base

        # Optional: an "effective F" (normalized drive fraction) that tracks rate (not required)
        num1 = w0*pre["Lc1"].to_numpy() + pre["wsite1"].to_numpy()*pre["Lpnp1"].to_numpy()
        num2 = w0*pre["Lc2"].to_numpy() + pre["wsite2"].to_numpy()*pre["Lpnp2"].to_numpy()
        den1 = w0 + pre["wsite1"].to_numpy() + pre["sig1"].to_numpy()
        den2 = w0 + pre["wsite2"].to_numpy() + pre["sig2"].to_numpy()
        F_eff1 = np.divide(num1, den1, out=np.zeros_like(num1), where=den1>0)
        F_eff2 = np.divide(num2, den2, out=np.zeros_like(num2), where=den2>0)
        pre["F_base_eff"] = np.sqrt(np.clip(F_eff1*F_eff2, 0.0, 1.0))

        # Handy boolean flags
        pre["is_paired"]        = (pre["stim_kind"] == 0).astype(int)
        pre["is_single_center"] = (pre["stim_kind"] == 1).astype(int)
        pre["is_single_periph"] = (pre["stim_kind"] == 2).astype(int)

        return pre


    # ---------- Main simulation (trial-by-trial weight gains, with singles) ----------
    def simulate_trialwise_weight_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    ):
        rng = np.random.default_rng(rng_seed)
        pre = precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup)

        # Pull out as numpy arrays for speed
        Lc1   = pre["Lc1"].to_numpy()
        Lpnp1 = pre["Lpnp1"].to_numpy()
        Lc2   = pre["Lc2"].to_numpy()
        Lpnp2 = pre["Lpnp2"].to_numpy()
        ws1   = pre["wsite1"].to_numpy()
        ws2   = pre["wsite2"].to_numpy()
        sig1  = pre["sig1"].to_numpy()
        sig2  = pre["sig2"].to_numpy()
        b1    = pre["base1"].to_numpy()
        b2    = pre["base2"].to_numpy()
        site  = pre["site"].to_numpy().astype(int)      # 0..3
        sgn   = pre["sign"].to_numpy().astype(float)    # +1 / -1
        bc    = pre["bhattCoef"].to_numpy()

        n_rows = len(pre)
        w0_base = 1.0

        # per-row shared fraction for the latent, possibly scaled by BC magnitude
        if USE_BC_SCALED_RHO:
            rho_eff = RHO_MIN + (RHO_MAX - RHO_MIN) * np.clip(bc, 0.0, 1.0)
        else:
            rho_eff = np.full(n_rows, RHO_SHARED, dtype=float)

        # Online accumulators per row
        n     = np.zeros(n_rows, dtype=float)
        m1    = np.zeros(n_rows, dtype=float)
        m2    = np.zeros(n_rows, dtype=float)
        M2_1  = np.zeros(n_rows, dtype=float)
        M2_2  = np.zeros(n_rows, dtype=float)
        C12   = np.zeros(n_rows, dtype=float)

        # Trial loop
        for _ in range(n_trials):
            # 5 global gains (shared by all rows/pairs) for this trial
            g0 = rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) if APPLY_W0_GAIN else 1.0
            g  = np.array([rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) for _ in range(4)], dtype=float)  # g1..g4

            # Select site gain per row
            g_site = g[site]  # shape (n_rows,)

            # Jittered weights for this trial
            w0_j   = w0_base * g0                     # same for all rows
            ws1_j  = ws1 * g_site
            ws2_j  = ws2 * g_site

            # Base rates (given this trial's weight gains), row-specific L's
            r1_base = rate_from_weights(Lc1, Lpnp1, w0_j, ws1_j, sig1, b1)
            r2_base = rate_from_weights(Lc2, Lpnp2, w0_j, ws2_j, sig2, b2)

            # Trialwise latent (BC-scaled shared fraction, BC-signed)
            eta1 = rng.normal(size=n_rows)
            eta2 = rng.normal(size=n_rows)
            xi   = rng.normal(size=n_rows)           # shared per row
            s1 = np.sqrt(1.0 - rho_eff)
            s2 = np.sqrt(rho_eff)

            lin1 = s1 * eta1 + s2 * xi
            lin2 = s1 * eta2 + sgn * s2 * xi

            gain1 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin1)
            gain2 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin2)

            # Trial rates with latent gain
            r1_t = r1_base * gain1
            r2_t = r2_base * gain2

            # Poisson counts for this trial
            lam1 = np.clip(r1_t * DT_SECONDS, 0, None)
            lam2 = np.clip(r2_t * DT_SECONDS, 0, None)
            k1 = rng.poisson(lam1)
            k2 = rng.poisson(lam2)

            # Online update
            n, m1, m2, M2_1, M2_2, C12 = online_update(k1, k2, n, m1, m2, M2_1, M2_2, C12)

        # Final stats
        dof  = np.maximum(1.0, n - 1.0)
        var1 = M2_1 / dof
        var2 = M2_2 / dof
        cov12 = C12 / dof
        denom = np.sqrt(np.clip(var1 * var2, 1e-12, None))
        r_sc  = np.where(denom > 0, cov12 / denom, np.nan)

        # Pack result DataFrame — one row per (pair × condition × stim_kind)
        out = pre.copy()
        out["mean1"] = m1
        out["mean2"] = m2
        out["var1"]  = var1
        out["var2"]  = var2
        out["cov12"] = cov12
        out["r_sc"]  = r_sc
        out["trials"] = n.astype(int)
        out["rho_eff"] = rho_eff

        return out

    # ================= USER SETTINGS =================
    N_TRIALS      = 10000   # trials per (pair × condition or single)
    DT_SECONDS    = 0.25    # 250 ms counting window

    # ---- Per-trial weight gains (global across all neurons & pairs) ----
    WGAIN_SHAPE   = 10.0    # Gamma(shape, scale=1/shape); mean=1, var=1/shape
    APPLY_W0_GAIN = True    # if False, g0=1 always

    # ---- Trialwise latent gain (imposes r_sc & sign by BC) ----
    SIGMA_GAIN    = 0.4     # log-SD of per-trial rate gain
    RHO_SHARED    = 0.4     # used if USE_BC_SCALED_RHO=False
    KAPPA         = 0.3     # 0 => gain fully modulates rate; 1 => cancels in denom

    # ---- Make the shared fraction depend on BC magnitude (optional) ----
    USE_BC_SCALED_RHO = True
    RHO_MIN = 0.20          # shared fraction when BC = 0
    RHO_MAX = 0.80          # shared fraction when BC = 1

    # ---- BC sign gate ----
    BC_THRESH     = 0.65    # ≥thr => positive shared sign; <thr => negative

    # ---- RNG ----
    RNG_SEED      = 123

    # ========== PREREQS (must already exist in your session) ==========
    # all_combs: list of (unit_i, unit_j)
    # masterGoodUnits: np.array of unit IDs (strings), aligned with rfWeightNeuronFits
    # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
    # bc_lookup[tuple(sorted((u1,u2)))] = BC float

    # ===================== RUN =====================
    full_df = simulate_trialwise_weight_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    )

    print(full_df[["stim_kind", "F_base", "bhattCoef", "r_sc", "mean1", "mean2"]].head())

    # ===================== SAVE (MAT for MATLAB meshes) =====================
    mat_path = "trialwise_weightnoise_with_singles.mat"
    numeric_cols = {
        k: full_df[k].to_numpy()
        for k in full_df.columns
        if np.issubdtype(full_df[k].dtype, np.number)
    }
    # Also save stim_kind mapping as comments here (MAT can load the numeric column):
    # stim_kind: 0=paired, 1=single_center, 2=single_periph
    savemat(mat_path, numeric_cols)
    print(f"Saved {mat_path}")

    # ===================== QUICK VISUAL (optional) =====================
    # Example: show r_sc over F×BC for paired rows only
    try:
        m = (full_df["stim_kind"] == 0)
        x = full_df.loc[m, "F_base"].to_numpy()
        y = full_df.loc[m, "bhattCoef"].to_numpy()
        z = full_df.loc[m, "r_sc"].to_numpy()
        ms = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        vlim = 0.15
        plt.figure(figsize=(6.6, 5))
        sc = plt.scatter(x[ms], y[ms], c=z[ms], s=10, cmap="coolwarm", vmin=-vlim, vmax=+vlim)
        plt.colorbar(sc, label="r_sc")
        plt.xlabel("F_base (geom. mean weight fraction)")
        plt.ylabel("Bhattacharyya coefficient (BC)")
        plt.title("Paired only: r_sc over F × BC")
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show()
    except Exception as e:
        print("Quick plot skipped:", e)

# new version
for filler in range(1):
    import numpy as np
    import pandas as pd
    from scipy.io import savemat
    import matplotlib.pyplot as plt

    # ================= USER SETTINGS =================
    N_TRIALS      = 10000   # trials per (pair × condition or single)
    DT_SECONDS    = 0.25    # 250 ms counting window

    # ---- Per-trial weight gains (global across all neurons & pairs) ----
    WGAIN_SHAPE   = 10.0    # Gamma(shape, scale=1/shape); mean=1, var=1/shape
    APPLY_W0_GAIN = True    # if False, g0=1 always

    # ---- Trialwise latent gain (imposes r_sc & sign by BC) ----
    SIGMA_GAIN    = 0.4     # log-SD of per-trial rate gain
    RHO_SHARED    = 0.4     # used if USE_BC_SCALED_RHO=False
    KAPPA         = 0.3     # 0 => gain fully modulates rate; 1 => cancels in denom

    # ---- Make the shared fraction depend on BC × drive (using F_base_eff) ----
    USE_BC_SCALED_RHO = True
    RHO_MIN = 0.20  # floor when BC=0 or drive is tiny
    RHO_MAX = 0.80  # ceiling when BC=1 and drive is large
    RHO_DRIVE_EXP = 1.5  # f(drive) = drive_norm ** RHO_DRIVE_EXP
    # Notes:
    # - We’ll use pre["F_base_eff"] as the drive proxy (already computed in your precompute_table)
    # - We robustly normalize it to [0,1] with 5–95th percentiles

    # ---- BC sign gate ----
    BC_THRESH     = 0.65    # ≥thr => positive shared sign; <thr => negative

    # ---- RNG ----
    RNG_SEED      = 123

    # ---------- helpers ----------
    def decode_paired_condition(i):
        """Return center_is_pref (bool), periph_is_pref (bool), site_idx (0..3)."""
        if 10 <= i <= 17:
            site = (i - 10) // 2
            periph_is_pref = (i - 10) % 2 == 1
            center_is_pref = False
        elif 19 <= i <= 26:
            site = (i - 19) // 2
            periph_is_pref = (i - 19) % 2 == 1
            center_is_pref = True
        else:
            raise ValueError(f"Condition {i} is not paired (10..17 or 19..26).")
        return center_is_pref, periph_is_pref, site


    def sign_from_bc(bc, thresh=BC_THRESH):
        return +1.0 if (np.isfinite(bc) and bc >= thresh) else -1.0


    def rate_from_weights(Lc, Lpnp, w0, ws, sigma, base):
        num = w0 * Lc + ws * Lpnp
        den = w0 + ws + np.maximum(1e-9, sigma)
        return num / den + base


    def weight_fraction(w0, ws, sigma):
        den = w0 + ws + np.maximum(1e-9, sigma)
        return np.clip((w0 + ws) / den, 0.0, 1.0)


    def fast_corr(x, y):
        sx = x.std(); sy = y.std()
        if sx == 0 or sy == 0:
            return np.nan
        return float(np.corrcoef(x, y)[0, 1])

    # Online (vectorized) update for means, variances, covariance
    def online_update(x, y, n, mean1, mean2, M2_1, M2_2, C12):
        """
        x,y are vectors (one sample per row)
        Updates Welford accumulators.
        """
        n_new = n + 1.0
        dx = x - mean1
        mean1_new = mean1 + dx / n_new
        dy = y - mean2
        mean2_new = mean2 + dy / n_new
        M2_1_new = M2_1 + dx * (x - mean1_new)
        M2_2_new = M2_2 + dy * (y - mean2_new)
        C12_new  = C12  + dx * (y - mean2_new)
        return n_new, mean1_new, mean2_new, M2_1_new, M2_2_new, C12_new

    # ---------- Precompute table (paired + derived singles) ----------
    def precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup):
        """
        Build a long design matrix with 3 rows per 'paired' condition:
          stim_kind = 0: paired (center + that periphery site shown together)
          stim_kind = 1: single_center (center shown, periphery blank -> Lpnp=0)
          stim_kind = 2: single_periph (that periphery site shown, center blank -> Lc=0)
        """
        params = {unit: rfWeightNeuronFits[i] for i, unit in enumerate(masterGoodUnits)}
        paired_indices = list(range(10, 18)) + list(range(19, 27))
        rows = []

        for (u1, u2) in all_combs:
            if (u1 not in params) or (u2 not in params):
                continue
            p1 = params[u1]; p2 = params[u2]
            Lp1, Lnp1, w1_1, w2_1, w3_1, w4_1, sig1, base1 = p1
            Lp2, Lnp2, w1_2, w2_2, w3_2, w4_2, sig2, base2 = p2
            key = tuple(sorted((u1, u2)))
            bc  = float(bc_lookup.get(key, np.nan))
            sgn = sign_from_bc(bc)

            for cond_idx in paired_indices:
                center_is_pref, periph_is_pref, site = decode_paired_condition(cond_idx)
                ws1 = [w1_1, w2_1, w3_1, w4_1][site]
                ws2 = [w1_2, w2_2, w3_2, w4_2][site]
                # L's for the paired condition
                Lc1_p   = Lp1 if center_is_pref  else Lnp1
                Lpnp1_p = Lp1 if periph_is_pref  else Lnp1
                Lc2_p   = Lp2 if center_is_pref  else Lnp2
                Lpnp2_p = Lp2 if periph_is_pref  else Lnp2

                # ---- (0) paired row ----
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, Lpnp1_p, Lc2_p, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 0))

                # ---- (1) single_center row (periphery blank -> Lpnp=0) ----
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, 0.0, Lc2_p, 0.0, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 1))

                # ---- (2) single_periph row (center blank -> Lc=0) ----
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             0.0, Lpnp1_p, 0.0, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 2))

        cols = ["unit_i","unit_j","cond_index","site","center_is_pref","periph_is_pref",
                "Lc1","Lpnp1","Lc2","Lpnp2","wsite1","wsite2","sig1","base1","sig2","base2",
                "bhattCoef","sign","stim_kind"]
        pre = pd.DataFrame(rows, columns=cols)

        # Baseline F (no weight noise), per row (same across the 3 rows of a trio)
        w0 = 1.0
        f1_base = weight_fraction(w0, pre["wsite1"].to_numpy(), pre["sig1"].to_numpy())
        f2_base = weight_fraction(w0, pre["wsite2"].to_numpy(), pre["sig2"].to_numpy())
        F_geo_base = np.sqrt(f1_base * f2_base)
        pre["F_base"] = F_geo_base

        # Optional: an "effective F" (normalized drive fraction) that tracks rate (not required)
        num1 = w0*pre["Lc1"].to_numpy() + pre["wsite1"].to_numpy()*pre["Lpnp1"].to_numpy()
        num2 = w0*pre["Lc2"].to_numpy() + pre["wsite2"].to_numpy()*pre["Lpnp2"].to_numpy()
        den1 = w0 + pre["wsite1"].to_numpy() + pre["sig1"].to_numpy()
        den2 = w0 + pre["wsite2"].to_numpy() + pre["sig2"].to_numpy()
        F_eff1 = np.divide(num1, den1, out=np.zeros_like(num1), where=den1>0)
        F_eff2 = np.divide(num2, den2, out=np.zeros_like(num2), where=den2>0)
        pre["F_base_eff"] = np.sqrt(np.clip(F_eff1*F_eff2, 0.0, 1.0))

        # Handy boolean flags
        pre["is_paired"]        = (pre["stim_kind"] == 0).astype(int)
        pre["is_single_center"] = (pre["stim_kind"] == 1).astype(int)
        pre["is_single_periph"] = (pre["stim_kind"] == 2).astype(int)

        return pre

    # ---------- Main simulation (trial-by-trial weight gains, with singles) ----------
    def simulate_trialwise_weight_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    ):
        rng = np.random.default_rng(rng_seed)
        pre = precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup)

        # Pull out as numpy arrays for speed
        Lc1   = pre["Lc1"].to_numpy()
        Lpnp1 = pre["Lpnp1"].to_numpy()
        Lc2   = pre["Lc2"].to_numpy()
        Lpnp2 = pre["Lpnp2"].to_numpy()
        ws1   = pre["wsite1"].to_numpy()
        ws2   = pre["wsite2"].to_numpy()
        sig1  = pre["sig1"].to_numpy()
        sig2  = pre["sig2"].to_numpy()
        b1    = pre["base1"].to_numpy()
        b2    = pre["base2"].to_numpy()
        site  = pre["site"].to_numpy().astype(int)      # 0..3
        sgn   = pre["sign"].to_numpy().astype(float)    # +1 / -1
        bc    = pre["bhattCoef"].to_numpy()

        n_rows = len(pre)
        w0_base = 1.0

        # per-row shared fraction for the latent, now scaled by BC × drive (F_base_eff)
        if USE_BC_SCALED_RHO:
            # drive proxy (already computed in precompute_table_with_singles)
            drive = pre["F_base_eff"].to_numpy()

            # robustly normalize drive to [0,1]
            lo, hi = np.nanpercentile(drive, [5, 95])
            span = max(1e-9, hi - lo)
            drive_norm = np.clip((drive - lo) / span, 0.0, 1.0) ** RHO_DRIVE_EXP

            # BC magnitude in [0,1]
            bc01 = np.clip(bc, 0.0, 1.0)

            # final per-row shared fraction
            rho_eff = RHO_MIN + (RHO_MAX - RHO_MIN) * (bc01 * drive_norm)
        else:
            rho_eff = np.full(n_rows, RHO_SHARED, dtype=float)

        # Online accumulators per row
        n     = np.zeros(n_rows, dtype=float)
        m1    = np.zeros(n_rows, dtype=float)
        m2    = np.zeros(n_rows, dtype=float)
        M2_1  = np.zeros(n_rows, dtype=float)
        M2_2  = np.zeros(n_rows, dtype=float)
        C12   = np.zeros(n_rows, dtype=float)

        # Trial loop
        for _ in range(n_trials):
            # 5 global gains (shared by all rows/pairs) for this trial
            g0 = rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) if APPLY_W0_GAIN else 1.0
            g  = np.array([rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) for _ in range(4)], dtype=float)  # g1..g4

            # Select site gain per row
            g_site = g[site]  # shape (n_rows,)

            # Jittered weights for this trial
            w0_j   = w0_base * g0                     # same for all rows
            ws1_j  = ws1 * g_site
            ws2_j  = ws2 * g_site

            # Base rates (given this trial's weight gains), row-specific L's
            r1_base = rate_from_weights(Lc1, Lpnp1, w0_j, ws1_j, sig1, b1)
            r2_base = rate_from_weights(Lc2, Lpnp2, w0_j, ws2_j, sig2, b2)

            # Trialwise latent (BC-scaled shared fraction, BC-signed)
            eta1 = rng.normal(size=n_rows)
            eta2 = rng.normal(size=n_rows)
            xi   = rng.normal(size=n_rows)           # shared per row
            s1 = np.sqrt(1.0 - rho_eff)
            s2 = np.sqrt(rho_eff)

            lin1 = s1 * eta1 + s2 * xi
            lin2 = s1 * eta2 + sgn * s2 * xi

            gain1 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin1)
            gain2 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin2)

            # Trial rates with latent gain
            r1_t = r1_base * gain1
            r2_t = r2_base * gain2

            # Poisson counts for this trial
            lam1 = np.clip(r1_t * DT_SECONDS, 0, None)
            lam2 = np.clip(r2_t * DT_SECONDS, 0, None)
            k1 = rng.poisson(lam1)
            k2 = rng.poisson(lam2)

            # Online update
            n, m1, m2, M2_1, M2_2, C12 = online_update(k1, k2, n, m1, m2, M2_1, M2_2, C12)

        # Final stats
        dof  = np.maximum(1.0, n - 1.0)
        var1 = M2_1 / dof
        var2 = M2_2 / dof
        cov12 = C12 / dof
        denom = np.sqrt(np.clip(var1 * var2, 1e-12, None))
        r_sc  = np.where(denom > 0, cov12 / denom, np.nan)

        # Pack result DataFrame — one row per (pair × condition × stim_kind)
        out = pre.copy()
        out["mean1"] = m1
        out["mean2"] = m2
        out["var1"]  = var1
        out["var2"]  = var2
        out["cov12"] = cov12
        out["r_sc"]  = r_sc
        out["trials"] = n.astype(int)
        out["rho_eff"] = rho_eff

        return out

    # ========== PREREQS (must already exist in your session) ==========
    # all_combs: list of (unit_i, unit_j)
    # masterGoodUnits: np.array of unit IDs (strings), aligned with rfWeightNeuronFits
    # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
    # bc_lookup[tuple(sorted((u1,u2)))] = BC float

    # ===================== RUN =====================
    full_df = simulate_trialwise_weight_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    )

    print(full_df[["stim_kind", "F_base", "bhattCoef", "r_sc", "mean1", "mean2"]].head())

    # ===================== SAVE (MAT for MATLAB meshes) =====================
    mat_path = "trialwise_weightnoise_with_singles.mat"
    numeric_cols = {
        k: full_df[k].to_numpy()
        for k in full_df.columns
        if np.issubdtype(full_df[k].dtype, np.number)
    }
    # Also save stim_kind mapping as comments here (MAT can load the numeric column):
    # stim_kind: 0=paired, 1=single_center, 2=single_periph
    savemat(mat_path, numeric_cols)
    print(f"Saved {mat_path}")

    # ===================== QUICK VISUAL (optional) =====================
    # Example: show r_sc over F×BC for paired rows only
    try:
        m = (full_df["stim_kind"] == 0)
        x = full_df.loc[m, "F_base"].to_numpy()
        y = full_df.loc[m, "bhattCoef"].to_numpy()
        z = full_df.loc[m, "r_sc"].to_numpy()
        ms = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        vlim = 0.15
        plt.figure(figsize=(6.6, 5))
        sc = plt.scatter(x[ms], y[ms], c=z[ms], s=10, cmap="coolwarm", vmin=-vlim, vmax=+vlim)
        plt.colorbar(sc, label="r_sc")
        plt.xlabel("F_base (geom. mean weight fraction)")
        plt.ylabel("Bhattacharyya coefficient (BC)")
        plt.title("Paired only: r_sc over F × BC")
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show()
    except Exception as e:
        print("Quick plot skipped:", e)


# with spon noise floor
# new version (with spontaneous/baseline noise floor)
for filler in range(1):
    import numpy as np
    import pandas as pd
    from scipy.io import savemat
    import matplotlib.pyplot as plt

    # ================= USER SETTINGS =================
    N_TRIALS      = 10000   # trials per (pair × condition or single)
    DT_SECONDS    = 0.25    # 250 ms counting window

    # ---- Per-trial weight gains (global across all neurons & pairs) ----
    WGAIN_SHAPE   = 3.0    # Gamma(shape, scale=1/shape); mean=1, var=1/shape
    APPLY_W0_GAIN = True    # if False, g0=1 always

    # ---- Trialwise latent gain (imposes r_sc & sign by BC) ----
    SIGMA_GAIN    = 0.0     # log-SD of per-trial rate gain
    RHO_SHARED    = 0.4     # used if USE_BC_SCALED_RHO=False
    KAPPA         = 0.3     # 0 => gain fully modulates rate; 1 => cancels in denom

    # ---- Make the shared fraction depend on BC × drive (using F_base_eff) ----
    USE_BC_SCALED_RHO = True
    RHO_MIN = 0.20  # floor when BC=0 or drive is tiny
    RHO_MAX = 0.80  # ceiling when BC=1 and drive is large
    RHO_DRIVE_EXP = 1.5  # f(drive) = drive_norm ** RHO_DRIVE_EXP

    # ---- BC sign gate ----
    BC_THRESH     = 0.65    # ≥thr => positive shared sign; <thr => negative

    # ---- RNG ----
    RNG_SEED      = 123

    # ---- Spontaneous / baseline noise floor ----
    # Two modes:
    #   "lognormal_mult": base_t = base * lognormal_factor   (most conservative)
    #   "lognormal_add" : base_t = base + max(0, base*(factor-1))  (pure floor-up noise)
    BASE_NOISE_MODE   = "lognormal_mult"
    SIGMA_BASE        = 0.0   # log-SD for baseline noise (0 => off)
    RHO_BASE_SHARED   = 0.30   # shared fraction for baseline noise (0..1)
    BASE_SIGNED_BY_BC = False  # if True, shared baseline uses BC sign gate (often keep False)
    BASE_NOISE_CLIP   = 0.0    # ensure baseline stays >= this

    # ---------- helpers ----------
    def decode_paired_condition(i):
        """Return center_is_pref (bool), periph_is_pref (bool), site_idx (0..3)."""
        if 10 <= i <= 17:
            site = (i - 10) // 2
            periph_is_pref = (i - 10) % 2 == 1
            center_is_pref = False
        elif 19 <= i <= 26:
            site = (i - 19) // 2
            periph_is_pref = (i - 19) % 2 == 1
            center_is_pref = True
        else:
            raise ValueError(f"Condition {i} is not paired (10..17 or 19..26).")
        return center_is_pref, periph_is_pref, site


    def sign_from_bc(bc, thresh=BC_THRESH):
        return +1.0 if (np.isfinite(bc) and bc >= thresh) else -1.0


    def rate_from_weights(Lc, Lpnp, w0, ws, sigma, base):
        num = w0 * Lc + ws * Lpnp
        den = w0 + ws + np.maximum(1e-9, sigma)
        return num / den + base


    def weight_fraction(w0, ws, sigma):
        den = w0 + ws + np.maximum(1e-9, sigma)
        return np.clip((w0 + ws) / den, 0.0, 1.0)


    # Online (vectorized) update for means, variances, covariance
    def online_update(x, y, n, mean1, mean2, M2_1, M2_2, C12):
        """
        x,y are vectors (one sample per row)
        Updates Welford accumulators.
        """
        n_new = n + 1.0
        dx = x - mean1
        mean1_new = mean1 + dx / n_new
        dy = y - mean2
        mean2_new = mean2 + dy / n_new
        M2_1_new = M2_1 + dx * (x - mean1_new)
        M2_2_new = M2_2 + dy * (y - mean2_new)
        C12_new  = C12  + dx * (y - mean2_new)
        return n_new, mean1_new, mean2_new, M2_1_new, M2_2_new, C12_new


    def add_spont_noise_to_base(b, rng, n_rows, rho_base, sigma_base,
                               sgn=None, signed=False, mode="lognormal_mult",
                               base_clip=0.0):
        """
        Trial-specific baseline noise model.

        Inputs
        ------
        b : array (n_rows,)
            Baseline (spontaneous rate) from fits.
        rho_base : float in [0,1]
            Shared fraction for baseline noise between neuron 1 and neuron 2.
        sigma_base : float >= 0
            Log-SD of the baseline noise.
        signed : bool
            If True, apply BC sign-gate to the shared term (rarely needed; often False).
        mode : "lognormal_mult" or "lognormal_add"
            - mult: base_t = b * f
            - add : base_t = b + max(0, b*(f-1))  (floor-up only)

        Returns
        -------
        b1_t, b2_t : arrays (n_rows,)
            Trial-specific baselines for neuron 1 and 2.
        """
        if sigma_base <= 0:
            return b, b

        # correlated standard normals
        e1 = rng.normal(size=n_rows)
        e2 = rng.normal(size=n_rows)
        z  = rng.normal(size=n_rows)  # shared

        s1 = np.sqrt(1.0 - rho_base)
        s2 = np.sqrt(rho_base)

        lin1 = s1 * e1 + s2 * z
        if signed and (sgn is not None):
            lin2 = s1 * e2 + (sgn * s2 * z)
        else:
            lin2 = s1 * e2 + s2 * z

        # lognormal factor with mean ~1
        f1 = np.exp(sigma_base * lin1 - 0.5 * sigma_base**2)
        f2 = np.exp(sigma_base * lin2 - 0.5 * sigma_base**2)

        if mode == "lognormal_mult":
            b1_t = np.clip(b * f1, base_clip, None)
            b2_t = np.clip(b * f2, base_clip, None)
        elif mode == "lognormal_add":
            add1 = np.clip(b * (f1 - 1.0), 0.0, None)
            add2 = np.clip(b * (f2 - 1.0), 0.0, None)
            b1_t = np.clip(b + add1, base_clip, None)
            b2_t = np.clip(b + add2, base_clip, None)
        else:
            raise ValueError("BASE_NOISE_MODE must be 'lognormal_mult' or 'lognormal_add'")

        return b1_t, b2_t


    # ---------- Precompute table (paired + derived singles) ----------
    def precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup):
        """
        Build a long design matrix with 3 rows per 'paired' condition:
          stim_kind = 0: paired (center + that periphery site shown together)
          stim_kind = 1: single_center (center shown, periphery blank -> Lpnp=0)
          stim_kind = 2: single_periph (that periphery site shown, center blank -> Lc=0)
        """
        params = {unit: rfWeightNeuronFits[i] for i, unit in enumerate(masterGoodUnits)}
        paired_indices = list(range(10, 18)) + list(range(19, 27))
        rows = []

        for (u1, u2) in all_combs:
            if (u1 not in params) or (u2 not in params):
                continue
            p1 = params[u1]; p2 = params[u2]
            Lp1, Lnp1, w1_1, w2_1, w3_1, w4_1, sig1, base1 = p1
            Lp2, Lnp2, w1_2, w2_2, w3_2, w4_2, sig2, base2 = p2
            key = tuple(sorted((u1, u2)))
            bc  = float(bc_lookup.get(key, np.nan))
            sgn = sign_from_bc(bc)

            for cond_idx in paired_indices:
                center_is_pref, periph_is_pref, site = decode_paired_condition(cond_idx)
                ws1 = [w1_1, w2_1, w3_1, w4_1][site]
                ws2 = [w1_2, w2_2, w3_2, w4_2][site]

                # L's for the paired condition
                Lc1_p   = Lp1 if center_is_pref  else Lnp1
                Lpnp1_p = Lp1 if periph_is_pref  else Lnp1
                Lc2_p   = Lp2 if center_is_pref  else Lnp2
                Lpnp2_p = Lp2 if periph_is_pref  else Lnp2

                # (0) paired
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, Lpnp1_p, Lc2_p, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 0))

                # (1) single_center (periphery blank)
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, 0.0, Lc2_p, 0.0, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 1))

                # (2) single_periph (center blank)
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             0.0, Lpnp1_p, 0.0, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 2))

        cols = ["unit_i","unit_j","cond_index","site","center_is_pref","periph_is_pref",
                "Lc1","Lpnp1","Lc2","Lpnp2","wsite1","wsite2","sig1","base1","sig2","base2",
                "bhattCoef","sign","stim_kind"]
        pre = pd.DataFrame(rows, columns=cols)

        # Baseline F (no weight noise)
        w0 = 1.0
        f1_base = weight_fraction(w0, pre["wsite1"].to_numpy(), pre["sig1"].to_numpy())
        f2_base = weight_fraction(w0, pre["wsite2"].to_numpy(), pre["sig2"].to_numpy())
        pre["F_base"] = np.sqrt(f1_base * f2_base)

        # Drive proxy (more rate-like)
        num1 = w0*pre["Lc1"].to_numpy() + pre["wsite1"].to_numpy()*pre["Lpnp1"].to_numpy()
        num2 = w0*pre["Lc2"].to_numpy() + pre["wsite2"].to_numpy()*pre["Lpnp2"].to_numpy()
        den1 = w0 + pre["wsite1"].to_numpy() + pre["sig1"].to_numpy()
        den2 = w0 + pre["wsite2"].to_numpy() + pre["sig2"].to_numpy()
        F_eff1 = np.divide(num1, den1, out=np.zeros_like(num1), where=den1>0)
        F_eff2 = np.divide(num2, den2, out=np.zeros_like(num2), where=den2>0)
        pre["F_base_eff"] = np.sqrt(np.clip(F_eff1*F_eff2, 0.0, 1.0))

        pre["is_paired"]        = (pre["stim_kind"] == 0).astype(int)
        pre["is_single_center"] = (pre["stim_kind"] == 1).astype(int)
        pre["is_single_periph"] = (pre["stim_kind"] == 2).astype(int)

        return pre


    # ---------- Main simulation ----------
    def simulate_trialwise_weight_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    ):
        rng = np.random.default_rng(rng_seed)
        pre = precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup)

        # numpy arrays for speed
        Lc1   = pre["Lc1"].to_numpy()
        Lpnp1 = pre["Lpnp1"].to_numpy()
        Lc2   = pre["Lc2"].to_numpy()
        Lpnp2 = pre["Lpnp2"].to_numpy()
        ws1   = pre["wsite1"].to_numpy()
        ws2   = pre["wsite2"].to_numpy()
        sig1  = pre["sig1"].to_numpy()
        sig2  = pre["sig2"].to_numpy()
        b1    = pre["base1"].to_numpy()
        b2    = pre["base2"].to_numpy()
        site  = pre["site"].to_numpy().astype(int)
        sgn   = pre["sign"].to_numpy().astype(float)
        bc    = pre["bhattCoef"].to_numpy()

        n_rows = len(pre)
        w0_base = 1.0

        # per-row shared fraction for latent gain, scaled by BC × drive
        if USE_BC_SCALED_RHO:
            drive = pre["F_base_eff"].to_numpy()
            lo, hi = np.nanpercentile(drive, [5, 95])
            span = max(1e-9, hi - lo)
            drive_norm = np.clip((drive - lo) / span, 0.0, 1.0) ** RHO_DRIVE_EXP
            bc01 = np.clip(bc, 0.0, 1.0)
            rho_eff = RHO_MIN + (RHO_MAX - RHO_MIN) * (bc01 * drive_norm)
        else:
            rho_eff = np.full(n_rows, RHO_SHARED, dtype=float)

        # Online accumulators per row
        n     = np.zeros(n_rows, dtype=float)
        m1    = np.zeros(n_rows, dtype=float)
        m2    = np.zeros(n_rows, dtype=float)
        M2_1  = np.zeros(n_rows, dtype=float)
        M2_2  = np.zeros(n_rows, dtype=float)
        C12   = np.zeros(n_rows, dtype=float)

        for _ in range(n_trials):
            # trialwise global weight gains
            g0 = rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) if APPLY_W0_GAIN else 1.0
            g  = np.array([rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) for _ in range(4)], dtype=float)
            g_site = g[site]

            # jittered weights
            w0_j  = w0_base * g0
            ws1_j = ws1 * g_site
            ws2_j = ws2 * g_site

            # ---- NEW: trial-specific spontaneous/baseline noise floor ----
            # Apply baseline noise separately to each neuron; can have shared component
            b1_t, b2_t = add_spont_noise_to_base(
                b=b1, rng=rng, n_rows=n_rows,
                rho_base=RHO_BASE_SHARED,
                sigma_base=SIGMA_BASE,
                sgn=sgn, signed=BASE_SIGNED_BY_BC,
                mode=BASE_NOISE_MODE,
                base_clip=BASE_NOISE_CLIP
            )

            # base rates including noisy baseline
            r1_base = rate_from_weights(Lc1, Lpnp1, w0_j, ws1_j, sig1, b1_t)
            r2_base = rate_from_weights(Lc2, Lpnp2, w0_j, ws2_j, sig2, b2_t)

            # latent gain (creates r_sc, with BC-scaled rho and BC-gated sign)
            eta1 = rng.normal(size=n_rows)
            eta2 = rng.normal(size=n_rows)
            xi   = rng.normal(size=n_rows)

            s1 = np.sqrt(1.0 - rho_eff)
            s2 = np.sqrt(rho_eff)

            lin1 = s1 * eta1 + s2 * xi
            lin2 = s1 * eta2 + sgn * s2 * xi

            gain1 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin1)
            gain2 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin2)

            r1_t = r1_base * gain1
            r2_t = r2_base * gain2

            lam1 = np.clip(r1_t * DT_SECONDS, 0, None)
            lam2 = np.clip(r2_t * DT_SECONDS, 0, None)
            k1 = rng.poisson(lam1)
            k2 = rng.poisson(lam2)

            n, m1, m2, M2_1, M2_2, C12 = online_update(k1, k2, n, m1, m2, M2_1, M2_2, C12)

        # final stats
        dof   = np.maximum(1.0, n - 1.0)
        var1  = M2_1 / dof
        var2  = M2_2 / dof
        cov12 = C12 / dof
        denom = np.sqrt(np.clip(var1 * var2, 1e-12, None))
        r_sc  = np.where(denom > 0, cov12 / denom, np.nan)

        out = pre.copy()
        out["mean1"]   = m1
        out["mean2"]   = m2
        out["var1"]    = var1
        out["var2"]    = var2
        out["cov12"]   = cov12
        out["r_sc"]    = r_sc
        out["trials"]  = n.astype(int)
        out["rho_eff"] = rho_eff

        # convenient: store baseline noise settings too
        out["sigma_base"]      = SIGMA_BASE
        out["rho_base_shared"] = RHO_BASE_SHARED

        return out


    # ========== PREREQS (must already exist in your session) ==========
    # all_combs: list of (unit_i, unit_j)
    # masterGoodUnits: np.array of unit IDs (strings), aligned with rfWeightNeuronFits
    # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
    # bc_lookup[tuple(sorted((u1,u2)))] = BC float

    # ===================== RUN =====================
    full_df = simulate_trialwise_weight_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    )

    print(full_df[["stim_kind", "F_base", "bhattCoef", "r_sc", "mean1", "mean2"]].head())

    # ===================== SAVE (MAT for MATLAB meshes) =====================
    mat_path = "trialwise_weightnoise_with_singles.mat"
    numeric_cols = {
        k: full_df[k].to_numpy()
        for k in full_df.columns
        if np.issubdtype(full_df[k].dtype, np.number)
    }
    # stim_kind: 0=paired, 1=single_center, 2=single_periph
    savemat(mat_path, numeric_cols)
    print(f"Saved {mat_path}")

    # ===================== QUICK VISUAL (optional) =====================
    try:
        m = (full_df["stim_kind"] == 0)
        x = full_df.loc[m, "F_base"].to_numpy()
        y = full_df.loc[m, "bhattCoef"].to_numpy()
        z = full_df.loc[m, "r_sc"].to_numpy()
        ms = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        vlim = 0.15
        plt.figure(figsize=(6.6, 5))
        sc = plt.scatter(x[ms], y[ms], c=z[ms], s=10, cmap="coolwarm", vmin=-vlim, vmax=+vlim)
        plt.colorbar(sc, label="r_sc")
        plt.xlabel("F_base (geom. mean weight fraction)")
        plt.ylabel("Bhattacharyya coefficient (BC)")
        plt.title("Paired only: r_sc over F × BC")
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show()
    except Exception as e:
        print("Quick plot skipped:", e)


# this version tries to replicate my results by adding a noisy normalization pool to the denom
for filler in range(1):
    import numpy as np
    import pandas as pd
    from scipy.io import savemat
    import matplotlib.pyplot as plt

    # ================= USER SETTINGS =================
    N_TRIALS      = 10000   # trials per (pair × condition or single)
    DT_SECONDS    = 0.25    # 250 ms counting window

    # ---- Per-trial weight gains (global across all neurons & pairs) ----
    WGAIN_SHAPE   = 10.0    # Gamma(shape, scale=1/shape); mean=1, var=1/shape
    APPLY_W0_GAIN = True    # if False, g0=1 always

    # ---- Trialwise latent gain (imposes r_sc & sign by BC) ----
    SIGMA_GAIN    = 0.0     # log-SD of per-trial rate gain
    RHO_SHARED    = 0.4     # used if USE_BC_SCALED_RHO=False
    KAPPA         = 0.3     # 0 => gain fully modulates rate; 1 => cancels in denom

    # ---- Make the shared fraction depend on BC × drive (using F_base_eff) ----
    USE_BC_SCALED_RHO = True
    RHO_MIN = 0.20  # floor when BC=0 or drive is tiny
    RHO_MAX = 0.80  # ceiling when BC=1 and drive is large
    RHO_DRIVE_EXP = 1.5  # f(drive) = drive_norm ** RHO_DRIVE_EXP

    # ---- BC sign gate ----
    BC_THRESH     = 0.65    # ≥thr => positive shared sign; <thr => negative

    # ---- RNG ----
    RNG_SEED      = 123

    # ---- Spontaneous / baseline noise floor ----
    BASE_NOISE_MODE   = "lognormal_mult"
    SIGMA_BASE        = 0.0   # log-SD for baseline noise (0 => off)
    RHO_BASE_SHARED   = 0.30   # shared fraction for baseline noise (0..1)
    BASE_SIGNED_BY_BC = False  # if True, shared baseline uses BC sign gate (often keep False)
    BASE_NOISE_CLIP   = 0.0    # ensure baseline stays >= this

    # ============================================================
    # NEW: Noisy normalization pool driven by (w0 + ws)
    #   r = E / (w0 + ws + sigma + alpha*(w0+ws)*N(t)) + b
    # ============================================================
    USE_NOISY_POOL   = True
    ALPHA_POOL       = 0.35   # coupling strength to pool noise (try 0.1–0.6)
    POOL_NOISE_SD    = 1.0    # N(t) ~ Normal(0, POOL_NOISE_SD)
    POOL_DEN_CLIP    = 1e-6   # keep denominator > 0

    # ---------- helpers ----------
    def decode_paired_condition(i):
        """Return center_is_pref (bool), periph_is_pref (bool), site_idx (0..3)."""
        if 10 <= i <= 17:
            site = (i - 10) // 2
            periph_is_pref = (i - 10) % 2 == 1
            center_is_pref = False
        elif 19 <= i <= 26:
            site = (i - 19) // 2
            periph_is_pref = (i - 19) % 2 == 1
            center_is_pref = True
        else:
            raise ValueError(f"Condition {i} is not paired (10..17 or 19..26).")
        return center_is_pref, periph_is_pref, site


    def sign_from_bc(bc, thresh=BC_THRESH):
        return +1.0 if (np.isfinite(bc) and bc >= thresh) else -1.0


    def rate_from_weights(Lc, Lpnp, w0, ws, sigma, base):
        """Original deterministic normalization (no pool noise)."""
        num = w0 * Lc + ws * Lpnp
        den = w0 + ws + np.maximum(1e-9, sigma)
        return num / den + base


    def rate_from_weights_noisy_pool(Lc, Lpnp, w0, ws, sigma, base, N_shared, alpha=ALPHA_POOL):
        """
        Noisy pool denominator:
            den = (w0 + ws + sigma) + alpha*(w0+ws)*N_shared
        N_shared is trialwise and shared between neuron 1 and 2 for each row.
        """
        num = w0 * Lc + ws * Lpnp
        wsum = (w0 + ws)
        den = wsum + sigma + alpha * wsum * N_shared
        den = np.clip(den, POOL_DEN_CLIP, None)
        return num / den + base


    def weight_fraction(w0, ws, sigma):
        den = w0 + ws + np.maximum(1e-9, sigma)
        return np.clip((w0 + ws) / den, 0.0, 1.0)


    # Online (vectorized) update for means, variances, covariance
    def online_update(x, y, n, mean1, mean2, M2_1, M2_2, C12):
        """
        x,y are vectors (one sample per row)
        Updates Welford accumulators.
        """
        n_new = n + 1.0
        dx = x - mean1
        mean1_new = mean1 + dx / n_new
        dy = y - mean2
        mean2_new = mean2 + dy / n_new
        M2_1_new = M2_1 + dx * (x - mean1_new)
        M2_2_new = M2_2 + dy * (y - mean2_new)
        C12_new  = C12  + dx * (y - mean2_new)
        return n_new, mean1_new, mean2_new, M2_1_new, M2_2_new, C12_new


    def add_spont_noise_to_base(b, rng, n_rows, rho_base, sigma_base,
                               sgn=None, signed=False, mode="lognormal_mult",
                               base_clip=0.0):
        """
        Trial-specific baseline noise model.

        Returns b1_t, b2_t : arrays (n_rows,)
        """
        if sigma_base <= 0:
            return b, b

        # correlated standard normals
        e1 = rng.normal(size=n_rows)
        e2 = rng.normal(size=n_rows)
        z  = rng.normal(size=n_rows)  # shared

        s1 = np.sqrt(1.0 - rho_base)
        s2 = np.sqrt(rho_base)

        lin1 = s1 * e1 + s2 * z
        if signed and (sgn is not None):
            lin2 = s1 * e2 + (sgn * s2 * z)
        else:
            lin2 = s1 * e2 + s2 * z

        # lognormal factor with mean ~1
        f1 = np.exp(sigma_base * lin1 - 0.5 * sigma_base**2)
        f2 = np.exp(sigma_base * lin2 - 0.5 * sigma_base**2)

        if mode == "lognormal_mult":
            b1_t = np.clip(b * f1, base_clip, None)
            b2_t = np.clip(b * f2, base_clip, None)
        elif mode == "lognormal_add":
            add1 = np.clip(b * (f1 - 1.0), 0.0, None)
            add2 = np.clip(b * (f2 - 1.0), 0.0, None)
            b1_t = np.clip(b + add1, base_clip, None)
            b2_t = np.clip(b + add2, base_clip, None)
        else:
            raise ValueError("BASE_NOISE_MODE must be 'lognormal_mult' or 'lognormal_add'")

        return b1_t, b2_t


    # ---------- Precompute table (paired + derived singles) ----------
    def precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup):
        """
        Build a long design matrix with 3 rows per 'paired' condition:
          stim_kind = 0: paired (center + that periphery site shown together)
          stim_kind = 1: single_center (center shown, periphery blank -> Lpnp=0)
          stim_kind = 2: single_periph (that periphery site shown, center blank -> Lc=0)
        """
        params = {unit: rfWeightNeuronFits[i] for i, unit in enumerate(masterGoodUnits)}
        paired_indices = list(range(10, 18)) + list(range(19, 27))
        rows = []

        for (u1, u2) in all_combs:
            if (u1 not in params) or (u2 not in params):
                continue
            p1 = params[u1]; p2 = params[u2]
            Lp1, Lnp1, w1_1, w2_1, w3_1, w4_1, sig1, base1 = p1
            Lp2, Lnp2, w1_2, w2_2, w3_2, w4_2, sig2, base2 = p2
            key = tuple(sorted((u1, u2)))
            bc  = float(bc_lookup.get(key, np.nan))
            sgn = sign_from_bc(bc)

            for cond_idx in paired_indices:
                center_is_pref, periph_is_pref, site = decode_paired_condition(cond_idx)
                ws1 = [w1_1, w2_1, w3_1, w4_1][site]
                ws2 = [w1_2, w2_2, w3_2, w4_2][site]

                # L's for the paired condition
                Lc1_p   = Lp1 if center_is_pref  else Lnp1
                Lpnp1_p = Lp1 if periph_is_pref  else Lnp1
                Lc2_p   = Lp2 if center_is_pref  else Lnp2
                Lpnp2_p = Lp2 if periph_is_pref  else Lnp2

                # (0) paired
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, Lpnp1_p, Lc2_p, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 0))

                # (1) single_center (periphery blank)
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, 0.0, Lc2_p, 0.0, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 1))

                # (2) single_periph (center blank)
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             0.0, Lpnp1_p, 0.0, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 2))

        cols = ["unit_i","unit_j","cond_index","site","center_is_pref","periph_is_pref",
                "Lc1","Lpnp1","Lc2","Lpnp2","wsite1","wsite2","sig1","base1","sig2","base2",
                "bhattCoef","sign","stim_kind"]
        pre = pd.DataFrame(rows, columns=cols)

        # Baseline F (no weight noise)
        w0 = 1.0
        f1_base = weight_fraction(w0, pre["wsite1"].to_numpy(), pre["sig1"].to_numpy())
        f2_base = weight_fraction(w0, pre["wsite2"].to_numpy(), pre["sig2"].to_numpy())
        pre["F_base"] = np.sqrt(f1_base * f2_base)

        # Drive proxy (more rate-like)
        num1 = w0*pre["Lc1"].to_numpy() + pre["wsite1"].to_numpy()*pre["Lpnp1"].to_numpy()
        num2 = w0*pre["Lc2"].to_numpy() + pre["wsite2"].to_numpy()*pre["Lpnp2"].to_numpy()
        den1 = w0 + pre["wsite1"].to_numpy() + pre["sig1"].to_numpy()
        den2 = w0 + pre["wsite2"].to_numpy() + pre["sig2"].to_numpy()
        F_eff1 = np.divide(num1, den1, out=np.zeros_like(num1), where=den1>0)
        F_eff2 = np.divide(num2, den2, out=np.zeros_like(num2), where=den2>0)
        pre["F_base_eff"] = np.sqrt(np.clip(F_eff1*F_eff2, 0.0, 1.0))

        pre["is_paired"]        = (pre["stim_kind"] == 0).astype(int)
        pre["is_single_center"] = (pre["stim_kind"] == 1).astype(int)
        pre["is_single_periph"] = (pre["stim_kind"] == 2).astype(int)

        return pre


    # ---------- Main simulation ----------
    def simulate_trialwise_weight_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    ):
        rng = np.random.default_rng(rng_seed)
        pre = precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup)

        # numpy arrays for speed
        Lc1   = pre["Lc1"].to_numpy()
        Lpnp1 = pre["Lpnp1"].to_numpy()
        Lc2   = pre["Lc2"].to_numpy()
        Lpnp2 = pre["Lpnp2"].to_numpy()
        ws1   = pre["wsite1"].to_numpy()
        ws2   = pre["wsite2"].to_numpy()
        sig1  = pre["sig1"].to_numpy()
        sig2  = pre["sig2"].to_numpy()
        b1    = pre["base1"].to_numpy()
        b2    = pre["base2"].to_numpy()
        site  = pre["site"].to_numpy().astype(int)
        sgn   = pre["sign"].to_numpy().astype(float)
        bc    = pre["bhattCoef"].to_numpy()

        n_rows = len(pre)
        w0_base = 1.0

        # per-row shared fraction for latent gain, scaled by BC × drive
        if USE_BC_SCALED_RHO:
            drive = pre["F_base_eff"].to_numpy()
            lo, hi = np.nanpercentile(drive, [5, 95])
            span = max(1e-9, hi - lo)
            drive_norm = np.clip((drive - lo) / span, 0.0, 1.0) ** RHO_DRIVE_EXP
            bc01 = np.clip(bc, 0.0, 1.0)
            rho_eff = RHO_MIN + (RHO_MAX - RHO_MIN) * (bc01 * drive_norm)
        else:
            rho_eff = np.full(n_rows, RHO_SHARED, dtype=float)

        # Online accumulators per row
        n     = np.zeros(n_rows, dtype=float)
        m1    = np.zeros(n_rows, dtype=float)
        m2    = np.zeros(n_rows, dtype=float)
        M2_1  = np.zeros(n_rows, dtype=float)
        M2_2  = np.zeros(n_rows, dtype=float)
        C12   = np.zeros(n_rows, dtype=float)

        for _ in range(n_trials):
            # trialwise global weight gains
            g0 = rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) if APPLY_W0_GAIN else 1.0
            g  = np.array([rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) for _ in range(4)], dtype=float)
            g_site = g[site]

            # jittered weights
            w0_j  = w0_base * g0
            ws1_j = ws1 * g_site
            ws2_j = ws2 * g_site

            # ---- trial-specific spontaneous/baseline noise floor ----
            b1_t, b2_t = add_spont_noise_to_base(
                b=b1, rng=rng, n_rows=n_rows,
                rho_base=RHO_BASE_SHARED,
                sigma_base=SIGMA_BASE,
                sgn=sgn, signed=BASE_SIGNED_BY_BC,
                mode=BASE_NOISE_MODE,
                base_clip=BASE_NOISE_CLIP
            )

            # ============================================================
            # NEW: trialwise noisy pool term N(t), shared between neurons
            # ============================================================
            # if USE_NOISY_POOL:
            #     N_shared = rng.normal(loc=0.0, scale=POOL_NOISE_SD, size=n_rows)
            #
            #     r1_base = rate_from_weights_noisy_pool(Lc1, Lpnp1, w0_j, ws1_j, sig1, b1_t, N_shared, alpha=ALPHA_POOL)
            #     r2_base = rate_from_weights_noisy_pool(Lc2, Lpnp2, w0_j, ws2_j, sig2, b2_t, N_shared, alpha=ALPHA_POOL)
            # else:
            #     r1_base = rate_from_weights(Lc1, Lpnp1, w0_j, ws1_j, sig1, b1_t)
            #     r2_base = rate_from_weights(Lc2, Lpnp2, w0_j, ws2_j, sig2, b2_t)

            if USE_NOISY_POOL:
                # trialwise pool noise, shared across the pair but BC-signed for neuron 2
                N_shared = rng.normal(loc=0.0, scale=POOL_NOISE_SD, size=n_rows)

                # apply BC sign-gate: high BC -> same sign; low BC -> opposite sign
                N1 = N_shared
                N2 = sgn * N_shared  # <-- THIS is the sign flip

                r1_base = rate_from_weights_noisy_pool(
                    Lc1, Lpnp1, w0_j, ws1_j, sig1, b1_t, N1, alpha=ALPHA_POOL
                )
                r2_base = rate_from_weights_noisy_pool(
                    Lc2, Lpnp2, w0_j, ws2_j, sig2, b2_t, N2, alpha=ALPHA_POOL
                )
            else:
                r1_base = rate_from_weights(Lc1, Lpnp1, w0_j, ws1_j, sig1, b1_t)
                r2_base = rate_from_weights(Lc2, Lpnp2, w0_j, ws2_j, sig2, b2_t)

            # latent gain (creates r_sc, with BC-scaled rho and BC-gated sign)
            eta1 = rng.normal(size=n_rows)
            eta2 = rng.normal(size=n_rows)
            xi   = rng.normal(size=n_rows)

            s1 = np.sqrt(1.0 - rho_eff)
            s2 = np.sqrt(rho_eff)

            lin1 = s1 * eta1 + s2 * xi
            lin2 = s1 * eta2 + sgn * s2 * xi

            gain1 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin1)
            gain2 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin2)

            r1_t = r1_base * gain1
            r2_t = r2_base * gain2

            lam1 = np.clip(r1_t * DT_SECONDS, 0, None)
            lam2 = np.clip(r2_t * DT_SECONDS, 0, None)
            k1 = rng.poisson(lam1)
            k2 = rng.poisson(lam2)

            n, m1, m2, M2_1, M2_2, C12 = online_update(k1, k2, n, m1, m2, M2_1, M2_2, C12)

        # final stats
        dof   = np.maximum(1.0, n - 1.0)
        var1  = M2_1 / dof
        var2  = M2_2 / dof
        cov12 = C12 / dof
        denom = np.sqrt(np.clip(var1 * var2, 1e-12, None))
        r_sc  = np.where(denom > 0, cov12 / denom, np.nan)

        out = pre.copy()
        out["mean1"]   = m1
        out["mean2"]   = m2
        out["var1"]    = var1
        out["var2"]    = var2
        out["cov12"]   = cov12
        out["r_sc"]    = r_sc
        out["trials"]  = n.astype(int)
        out["rho_eff"] = rho_eff

        # convenient: store baseline noise settings too
        out["sigma_base"]      = SIGMA_BASE
        out["rho_base_shared"] = RHO_BASE_SHARED

        # store noisy-pool settings
        out["use_noisy_pool"]  = int(USE_NOISY_POOL)
        out["alpha_pool"]      = float(ALPHA_POOL)
        out["pool_noise_sd"]   = float(POOL_NOISE_SD)

        return out


    # ========== PREREQS (must already exist in your session) ==========
    # all_combs: list of (unit_i, unit_j)
    # masterGoodUnits: np.array of unit IDs (strings), aligned with rfWeightNeuronFits
    # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
    # bc_lookup[tuple(sorted((u1,u2)))] = BC float

    # ===================== RUN =====================
    full_df = simulate_trialwise_weight_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    )

    print(full_df[["stim_kind", "F_base", "bhattCoef", "r_sc", "mean1", "mean2"]].head())

    # ===================== SAVE (MAT for MATLAB meshes) =====================
    mat_path = "trialwise_weightnoise_with_singles_version2.mat"
    numeric_cols = {
        k: full_df[k].to_numpy()
        for k in full_df.columns
        if np.issubdtype(full_df[k].dtype, np.number)
    }
    # stim_kind: 0=paired, 1=single_center, 2=single_periph
    savemat(mat_path, numeric_cols)
    print(f"Saved {mat_path}")


# this version adds a noisy normalization pool to the denom that does not scale with wo+ws
# ============================================================
# Trialwise weight-noise + "fair" noisy normalization pool (NO (w0+ws) multiplier)
#
# Key change vs version2:
#   den = (w0 + ws + sigma) + alpha_pool * N_shared
#   (i.e., pool noise amplitude does NOT scale with (w0+ws))
#
# Still supports:
# - optional BC sign-gate for the shared pool term (keep if you want; set to False to disable)
# - optional latent gain block (default off; SIGMA_GAIN=0)
# - optional baseline noise floor (default off; SIGMA_BASE=0)
#
# Saves: trialwise_weightnoise_with_singles_version3.mat
# ============================================================
for filler in range(1):
    import numpy as np
    import pandas as pd
    from scipy.io import savemat
    import matplotlib.pyplot as plt

    # ================= USER SETTINGS =================
    N_TRIALS      = 10000   # trials per (pair × condition or single)
    DT_SECONDS    = 0.25    # 250 ms counting window

    # ---- Per-trial weight gains (global across all neurons & pairs) ----
    WGAIN_SHAPE   = 10.0    # Gamma(shape, scale=1/shape); mean=1, var=1/shape
    APPLY_W0_GAIN = False    # if False, g0=1 always

    APPLY_SITE_WEIGHT_NOISE = False

    # ---- Trialwise latent gain (imposes r_sc & sign by BC) ----
    SIGMA_GAIN    = 0.0     # log-SD of per-trial rate gain (0 => OFF)
    RHO_SHARED    = 0.4     # used if USE_BC_SCALED_RHO=False
    KAPPA         = 0.3     # 0 => gain fully modulates rate; 1 => cancels in denom (kept for compatibility)

    # ---- Make the shared fraction depend on BC × drive (using F_base_eff) ----
    USE_BC_SCALED_RHO = True
    RHO_MIN = 0.20
    RHO_MAX = 0.80
    RHO_DRIVE_EXP = 1.5

    # ---- BC sign gate ----
    BC_THRESH     = 0.65    # ≥thr => positive shared sign; <thr => negative

    # ---- RNG ----
    RNG_SEED      = 123

    # ---- Spontaneous / baseline noise floor ----
    BASE_NOISE_MODE   = "lognormal_mult"
    SIGMA_BASE        = 0.0    # log-SD for baseline noise (0 => off)
    RHO_BASE_SHARED   = 0.30
    BASE_SIGNED_BY_BC = False
    BASE_NOISE_CLIP   = 0.0

    # ============================================================
    # FAIR noisy normalization pool (NO (w0+ws) multiplier)
    #
    #   r = E / (w0 + ws + sigma + alpha_pool*N(t)) + b
    #
    # Notes:
    # - alpha_pool has units of "denominator units"
    # - N(t) ~ Normal(0, POOL_NOISE_SD)
    # ============================================================
    USE_NOISY_POOL     = True
    ALPHA_POOL         = 0.35   # absolute coupling strength (not scaled by w0+ws)
    POOL_NOISE_SD      = 1.0
    POOL_DEN_CLIP      = 1e-6   # keep denominator > 0

    # Optional: keep / disable BC sign-gating for pool noise.
    # If False: both neurons receive same N_shared (no sign flip).
    # If True : neuron2 receives sgn*N_shared (high BC => +, low BC => -).
    POOL_SIGNED_BY_BC  = True

    # ---------- helpers ----------
    def decode_paired_condition(i):
        """Return center_is_pref (bool), periph_is_pref (bool), site_idx (0..3)."""
        if 10 <= i <= 17:
            site = (i - 10) // 2
            periph_is_pref = (i - 10) % 2 == 1
            center_is_pref = False
        elif 19 <= i <= 26:
            site = (i - 19) // 2
            periph_is_pref = (i - 19) % 2 == 1
            center_is_pref = True
        else:
            raise ValueError(f"Condition {i} is not paired (10..17 or 19..26).")
        return center_is_pref, periph_is_pref, site

    def sign_from_bc(bc, thresh=BC_THRESH):
        return +1.0 if (np.isfinite(bc) and bc >= thresh) else -1.0

    def rate_from_weights(Lc, Lpnp, w0, ws, sigma, base):
        """Original deterministic normalization (no pool noise)."""
        num = w0 * Lc + ws * Lpnp
        den = w0 + ws + np.maximum(1e-9, sigma)
        return num / den + base

    def rate_from_weights_noisy_pool_fair(Lc, Lpnp, w0, ws, sigma, base, N_shared, alpha=ALPHA_POOL):
        """
        FAIR noisy pool denominator:
            den = (w0 + ws + sigma) + alpha*N_shared
        (NO multiplier by (w0+ws))
        """
        num = w0 * Lc + ws * Lpnp
        den = (w0 + ws + sigma) + alpha * N_shared
        den = np.clip(den, POOL_DEN_CLIP, None)
        return num / den + base

    def weight_fraction(w0, ws, sigma):
        den = w0 + ws + np.maximum(1e-9, sigma)
        return np.clip((w0 + ws) / den, 0.0, 1.0)

    # Online (vectorized) update for means, variances, covariance
    def online_update(x, y, n, mean1, mean2, M2_1, M2_2, C12):
        n_new = n + 1.0
        dx = x - mean1
        mean1_new = mean1 + dx / n_new
        dy = y - mean2
        mean2_new = mean2 + dy / n_new
        M2_1_new = M2_1 + dx * (x - mean1_new)
        M2_2_new = M2_2 + dy * (y - mean2_new)
        C12_new  = C12  + dx * (y - mean2_new)
        return n_new, mean1_new, mean2_new, M2_1_new, M2_2_new, C12_new

    def add_spont_noise_to_base(b, rng, n_rows, rho_base, sigma_base,
                               sgn=None, signed=False, mode="lognormal_mult",
                               base_clip=0.0):
        """Trial-specific baseline noise model. Returns b1_t, b2_t."""
        if sigma_base <= 0:
            return b, b

        e1 = rng.normal(size=n_rows)
        e2 = rng.normal(size=n_rows)
        z  = rng.normal(size=n_rows)

        s1 = np.sqrt(1.0 - rho_base)
        s2 = np.sqrt(rho_base)

        lin1 = s1 * e1 + s2 * z
        if signed and (sgn is not None):
            lin2 = s1 * e2 + (sgn * s2 * z)
        else:
            lin2 = s1 * e2 + s2 * z

        f1 = np.exp(sigma_base * lin1 - 0.5 * sigma_base**2)
        f2 = np.exp(sigma_base * lin2 - 0.5 * sigma_base**2)

        if mode == "lognormal_mult":
            b1_t = np.clip(b * f1, base_clip, None)
            b2_t = np.clip(b * f2, base_clip, None)
        elif mode == "lognormal_add":
            add1 = np.clip(b * (f1 - 1.0), 0.0, None)
            add2 = np.clip(b * (f2 - 1.0), 0.0, None)
            b1_t = np.clip(b + add1, base_clip, None)
            b2_t = np.clip(b + add2, base_clip, None)
        else:
            raise ValueError("BASE_NOISE_MODE must be 'lognormal_mult' or 'lognormal_add'")
        return b1_t, b2_t

    # ---------- Precompute table (paired + derived singles) ----------
    def precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup):
        """
        Build a long design matrix with 3 rows per 'paired' condition:
          stim_kind = 0: paired
          stim_kind = 1: single_center
          stim_kind = 2: single_periph
        """
        params = {unit: rfWeightNeuronFits[i] for i, unit in enumerate(masterGoodUnits)}
        paired_indices = list(range(10, 18)) + list(range(19, 27))
        rows = []

        for (u1, u2) in all_combs:
            if (u1 not in params) or (u2 not in params):
                continue
            p1 = params[u1]; p2 = params[u2]
            Lp1, Lnp1, w1_1, w2_1, w3_1, w4_1, sig1, base1 = p1
            Lp2, Lnp2, w1_2, w2_2, w3_2, w4_2, sig2, base2 = p2

            key = tuple(sorted((u1, u2)))
            bc  = float(bc_lookup.get(key, np.nan))
            sgn = sign_from_bc(bc)

            for cond_idx in paired_indices:
                center_is_pref, periph_is_pref, site = decode_paired_condition(cond_idx)
                ws1 = [w1_1, w2_1, w3_1, w4_1][site]
                ws2 = [w1_2, w2_2, w3_2, w4_2][site]

                Lc1_p   = Lp1 if center_is_pref  else Lnp1
                Lpnp1_p = Lp1 if periph_is_pref  else Lnp1
                Lc2_p   = Lp2 if center_is_pref  else Lnp2
                Lpnp2_p = Lp2 if periph_is_pref  else Lnp2

                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, Lpnp1_p, Lc2_p, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 0))

                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, 0.0, Lc2_p, 0.0, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 1))

                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             0.0, Lpnp1_p, 0.0, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 2))

        cols = ["unit_i","unit_j","cond_index","site","center_is_pref","periph_is_pref",
                "Lc1","Lpnp1","Lc2","Lpnp2","wsite1","wsite2","sig1","base1","sig2","base2",
                "bhattCoef","sign","stim_kind"]
        pre = pd.DataFrame(rows, columns=cols)

        # Baseline F (no weight noise)
        w0 = 1.0
        f1_base = weight_fraction(w0, pre["wsite1"].to_numpy(), pre["sig1"].to_numpy())
        f2_base = weight_fraction(w0, pre["wsite2"].to_numpy(), pre["sig2"].to_numpy())
        pre["F_base"] = np.sqrt(f1_base * f2_base)

        # Drive proxy (more rate-like)
        num1 = w0*pre["Lc1"].to_numpy() + pre["wsite1"].to_numpy()*pre["Lpnp1"].to_numpy()
        num2 = w0*pre["Lc2"].to_numpy() + pre["wsite2"].to_numpy()*pre["Lpnp2"].to_numpy()
        den1 = w0 + pre["wsite1"].to_numpy() + pre["sig1"].to_numpy()
        den2 = w0 + pre["wsite2"].to_numpy() + pre["sig2"].to_numpy()
        F_eff1 = np.divide(num1, den1, out=np.zeros_like(num1), where=den1>0)
        F_eff2 = np.divide(num2, den2, out=np.zeros_like(num2), where=den2>0)
        pre["F_base_eff"] = np.sqrt(np.clip(F_eff1*F_eff2, 0.0, 1.0))

        pre["is_paired"]        = (pre["stim_kind"] == 0).astype(int)
        pre["is_single_center"] = (pre["stim_kind"] == 1).astype(int)
        pre["is_single_periph"] = (pre["stim_kind"] == 2).astype(int)

        # MATLAB compatibility fallback
        pre["nonPrefSupp"] = pre["F_base"]

        return pre

    # ---------- Main simulation ----------
    def simulate_trialwise_weight_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    ):
        rng = np.random.default_rng(rng_seed)
        pre = precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup)

        Lc1   = pre["Lc1"].to_numpy()
        Lpnp1 = pre["Lpnp1"].to_numpy()
        Lc2   = pre["Lc2"].to_numpy()
        Lpnp2 = pre["Lpnp2"].to_numpy()
        ws1   = pre["wsite1"].to_numpy()
        ws2   = pre["wsite2"].to_numpy()
        sig1  = pre["sig1"].to_numpy()
        sig2  = pre["sig2"].to_numpy()
        b1    = pre["base1"].to_numpy()
        b2    = pre["base2"].to_numpy()
        site  = pre["site"].to_numpy().astype(int)
        sgn   = pre["sign"].to_numpy().astype(float)
        bc    = pre["bhattCoef"].to_numpy()

        n_rows = len(pre)
        w0_base = 1.0

        # per-row shared fraction for latent gain, scaled by BC × drive
        if USE_BC_SCALED_RHO:
            drive = pre["F_base_eff"].to_numpy()
            lo, hi = np.nanpercentile(drive, [5, 95])
            span = max(1e-9, hi - lo)
            drive_norm = np.clip((drive - lo) / span, 0.0, 1.0) ** RHO_DRIVE_EXP
            bc01 = np.clip(bc, 0.0, 1.0)
            rho_eff = RHO_MIN + (RHO_MAX - RHO_MIN) * (bc01 * drive_norm)
        else:
            rho_eff = np.full(n_rows, RHO_SHARED, dtype=float)

        # Online accumulators per row
        n     = np.zeros(n_rows, dtype=float)
        m1    = np.zeros(n_rows, dtype=float)
        m2    = np.zeros(n_rows, dtype=float)
        M2_1  = np.zeros(n_rows, dtype=float)
        M2_2  = np.zeros(n_rows, dtype=float)
        C12   = np.zeros(n_rows, dtype=float)

        for _ in range(n_trials):
            # trialwise global weight gains
            if APPLY_SITE_WEIGHT_NOISE:
                g0 = rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) if APPLY_W0_GAIN else 1.0
                g  = np.array([rng.gamma(WGAIN_SHAPE, 1.0/WGAIN_SHAPE) for _ in range(4)], dtype=float)
                g_site = g[site]
            else:
                g0 = 1.0
                g_site = 1.0

            # jittered weights
            w0_j  = w0_base * g0
            ws1_j = ws1 * g_site
            ws2_j = ws2 * g_site

            # ---- trial-specific spontaneous/baseline noise floor ----
            b1_t, b2_t = add_spont_noise_to_base(
                b=b1, rng=rng, n_rows=n_rows,
                rho_base=RHO_BASE_SHARED,
                sigma_base=SIGMA_BASE,
                sgn=sgn, signed=BASE_SIGNED_BY_BC,
                mode=BASE_NOISE_MODE,
                base_clip=BASE_NOISE_CLIP
            )

            # ============================================================
            # FAIR noisy pool: shared N(t) enters denominator additively
            # ============================================================
            if USE_NOISY_POOL:
                N_shared = rng.normal(loc=0.0, scale=POOL_NOISE_SD, size=n_rows)
                if POOL_SIGNED_BY_BC:
                    N1 = N_shared
                    N2 = sgn * N_shared
                else:
                    N1 = N_shared
                    N2 = N_shared

                r1_base = rate_from_weights_noisy_pool_fair(
                    Lc1, Lpnp1, w0_j, ws1_j, sig1, b1_t, N1, alpha=ALPHA_POOL
                )
                r2_base = rate_from_weights_noisy_pool_fair(
                    Lc2, Lpnp2, w0_j, ws2_j, sig2, b2_t, N2, alpha=ALPHA_POOL
                )
            else:
                r1_base = rate_from_weights(Lc1, Lpnp1, w0_j, ws1_j, sig1, b1_t)
                r2_base = rate_from_weights(Lc2, Lpnp2, w0_j, ws2_j, sig2, b2_t)

            # latent gain (optional; default off when SIGMA_GAIN=0)
            if SIGMA_GAIN > 0:
                eta1 = rng.normal(size=n_rows)
                eta2 = rng.normal(size=n_rows)
                xi   = rng.normal(size=n_rows)

                s1 = np.sqrt(1.0 - rho_eff)
                s2 = np.sqrt(rho_eff)

                lin1 = s1 * eta1 + s2 * xi
                lin2 = s1 * eta2 + sgn * s2 * xi

                gain1 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin1)
                gain2 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin2)

                r1_t = r1_base * gain1
                r2_t = r2_base * gain2
            else:
                r1_t = r1_base
                r2_t = r2_base

            lam1 = np.clip(r1_t * DT_SECONDS, 0, None)
            lam2 = np.clip(r2_t * DT_SECONDS, 0, None)
            k1 = rng.poisson(lam1)
            k2 = rng.poisson(lam2)

            n, m1, m2, M2_1, M2_2, C12 = online_update(k1, k2, n, m1, m2, M2_1, M2_2, C12)

        # final stats
        dof   = np.maximum(1.0, n - 1.0)
        var1  = M2_1 / dof
        var2  = M2_2 / dof
        cov12 = C12 / dof
        denom = np.sqrt(np.clip(var1 * var2, 1e-12, None))
        r_sc  = np.where(denom > 0, cov12 / denom, np.nan)

        out = pre.copy()
        out["mean1"]   = m1
        out["mean2"]   = m2
        out["var1"]    = var1
        out["var2"]    = var2
        out["cov12"]   = cov12
        out["r_sc"]    = r_sc
        out["trials"]  = n.astype(int)
        out["rho_eff"] = rho_eff

        # store baseline noise settings too
        out["sigma_base"]      = SIGMA_BASE
        out["rho_base_shared"] = RHO_BASE_SHARED

        # store noisy-pool settings
        out["use_noisy_pool"]      = int(USE_NOISY_POOL)
        out["alpha_pool"]          = float(ALPHA_POOL)
        out["pool_noise_sd"]       = float(POOL_NOISE_SD)
        out["pool_signed_by_bc"]   = int(POOL_SIGNED_BY_BC)
        out["pool_den_clip"]       = float(POOL_DEN_CLIP)

        # store weight-gain settings
        out["wgain_shape"]     = float(WGAIN_SHAPE)
        out["apply_w0_gain"]   = int(APPLY_W0_GAIN)

        return out

    # ========== PREREQS (must already exist in your session) ==========
    # all_combs: list of (unit_i, unit_j)
    # masterGoodUnits: np.array of unit IDs (strings), aligned with rfWeightNeuronFits
    # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
    # bc_lookup[tuple(sorted((u1,u2)))] = BC float

    # ===================== RUN =====================
    full_df = simulate_trialwise_weight_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    )

    print(full_df[["stim_kind", "F_base", "bhattCoef", "r_sc", "mean1", "mean2"]].head())

    # ===================== SAVE (MAT for MATLAB meshes) =====================
    mat_path = "trialwise_weightnoise_with_singles_version3.mat"
    numeric_cols = {
        k: full_df[k].to_numpy()
        for k in full_df.columns
        if np.issubdtype(full_df[k].dtype, np.number)
    }
    savemat(mat_path, numeric_cols)
    print(f"Saved {mat_path}")

    # ===================== QUICK VISUAL (optional) =====================
    try:
        m = (full_df["stim_kind"] == 0)
        x = full_df.loc[m, "F_base"].to_numpy()
        y = full_df.loc[m, "bhattCoef"].to_numpy()
        z = full_df.loc[m, "r_sc"].to_numpy()
        ms = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        vlim = 0.15
        plt.figure(figsize=(6.6, 5))
        sc = plt.scatter(x[ms], y[ms], c=z[ms], s=10, cmap="coolwarm", vmin=-vlim, vmax=+vlim)
        plt.colorbar(sc, label="r_sc")
        plt.xlabel("F_base (geom. mean weight fraction)")
        plt.ylabel("Bhattacharyya coefficient (BC)")
        plt.title("Paired only: r_sc over F × BC (FAIR pool noise)")
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show()
    except Exception as e:
        print("Quick plot skipped:", e)


# # this version only applied multiplicative gain to the denominator including sigma + weights
# ============================================================
# Trialwise DENOM-only intensity gain (weights + sigma/ib ONLY in denominator),
# with singles, NO numerator gain.
#
# Response equation used:
#   R = (w0*Lc + ws*Lpnp) / ( (w0*g_den) + (ws*g_den) + (sigma*g_den) ) + base
#
# i.e. the gain ONLY affects the denominator "intensity" terms.
#
# Output MAT: trialwise_weightnoiseDenom_with_singles.mat
# ============================================================
for filler in range(1):
    import numpy as np
    import pandas as pd
    from scipy.io import savemat
    import matplotlib.pyplot as plt

    # ================= USER SETTINGS =================
    N_TRIALS      = 10000   # trials per (pair × condition or single)
    DT_SECONDS    = 0.25    # 250 ms counting window

    # ---- GLOBAL per-trial denom "intensity" gain ----
    # Gamma(shape, scale=1/shape): mean=1, var=1/shape
    DENOM_GAIN_SHAPE = 1.5    # smaller => stronger fluctuations => stronger correlations
    APPLY_DENOM_GAIN = True

    # ---- Trialwise latent gain (OFF by default; keep 0 to isolate denom gain) ----
    SIGMA_GAIN    = 0.0
    RHO_SHARED    = 0.4
    KAPPA         = 0.3

    # ---- Make the shared fraction depend on BC × drive (using F_base_eff) ----
    USE_BC_SCALED_RHO = True
    RHO_MIN = 0.20
    RHO_MAX = 0.80
    RHO_DRIVE_EXP = 1.5

    # ---- BC sign gate ----
    BC_THRESH     = 0.65

    # ---- RNG ----
    RNG_SEED      = 123

    # ---- Spontaneous / baseline noise floor ----
    BASE_NOISE_MODE   = "lognormal_mult"
    SIGMA_BASE        = 0.0
    RHO_BASE_SHARED   = 0.30
    BASE_SIGNED_BY_BC = False
    BASE_NOISE_CLIP   = 0.0

    # ---- MATLAB compatibility ----
    WRITE_NONPREFSUPP = True

    # ---------- helpers ----------
    def decode_paired_condition(i):
        """Return center_is_pref (bool), periph_is_pref (bool), site_idx (0..3)."""
        if 10 <= i <= 17:
            site = (i - 10) // 2
            periph_is_pref = (i - 10) % 2 == 1
            center_is_pref = False
        elif 19 <= i <= 26:
            site = (i - 19) // 2
            periph_is_pref = (i - 19) % 2 == 1
            center_is_pref = True
        else:
            raise ValueError(f"Condition {i} is not paired (10..17 or 19..26).")
        return center_is_pref, periph_is_pref, site


    def sign_from_bc(bc, thresh=BC_THRESH):
        return +1.0 if (np.isfinite(bc) and bc >= thresh) else -1.0


    def rate_from_weights_denom_only_gain(Lc, Lpnp, w0, ws, sigma, base, g_den):
        """
        Denominator-only gain:
          R = (w0*Lc + ws*Lpnp) / (g_den*(w0 + ws + sigma)) + base

        IMPORTANT: numerator uses UN-gained weights.
        """
        sig = np.maximum(1e-9, sigma)
        den = g_den * (w0 + ws + sig)
        num = (w0 * Lc) + (ws * Lpnp)
        return np.divide(num, den, out=np.zeros_like(num, dtype=float), where=den > 0) + base


    def weight_fraction(w0, ws, sigma):
        den = w0 + ws + np.maximum(1e-9, sigma)
        return np.clip((w0 + ws) / den, 0.0, 1.0)


    def online_update(x, y, n, mean1, mean2, M2_1, M2_2, C12):
        """Vectorized Welford update."""
        n_new = n + 1.0
        dx = x - mean1
        mean1_new = mean1 + dx / n_new
        dy = y - mean2
        mean2_new = mean2 + dy / n_new
        M2_1_new = M2_1 + dx * (x - mean1_new)
        M2_2_new = M2_2 + dy * (y - mean2_new)
        C12_new  = C12  + dx * (y - mean2_new)
        return n_new, mean1_new, mean2_new, M2_1_new, M2_2_new, C12_new


    def add_spont_noise_to_base(b, rng, n_rows, rho_base, sigma_base,
                               sgn=None, signed=False, mode="lognormal_mult",
                               base_clip=0.0):
        """Trial-specific baseline noise model."""
        if sigma_base <= 0:
            return b, b

        e1 = rng.normal(size=n_rows)
        e2 = rng.normal(size=n_rows)
        z  = rng.normal(size=n_rows)  # shared

        s1 = np.sqrt(1.0 - rho_base)
        s2 = np.sqrt(rho_base)

        lin1 = s1 * e1 + s2 * z
        if signed and (sgn is not None):
            lin2 = s1 * e2 + (sgn * s2 * z)
        else:
            lin2 = s1 * e2 + s2 * z

        f1 = np.exp(sigma_base * lin1 - 0.5 * sigma_base**2)
        f2 = np.exp(sigma_base * lin2 - 0.5 * sigma_base**2)

        if mode == "lognormal_mult":
            b1_t = np.clip(b * f1, base_clip, None)
            b2_t = np.clip(b * f2, base_clip, None)
        elif mode == "lognormal_add":
            add1 = np.clip(b * (f1 - 1.0), 0.0, None)
            add2 = np.clip(b * (f2 - 1.0), 0.0, None)
            b1_t = np.clip(b + add1, base_clip, None)
            b2_t = np.clip(b + add2, base_clip, None)
        else:
            raise ValueError("BASE_NOISE_MODE must be 'lognormal_mult' or 'lognormal_add'")

        return b1_t, b2_t


    # ---------- Precompute table (paired + derived singles) ----------
    def precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup):
        params = {unit: rfWeightNeuronFits[i] for i, unit in enumerate(masterGoodUnits)}
        paired_indices = list(range(10, 18)) + list(range(19, 27))
        rows = []

        for (u1, u2) in all_combs:
            if (u1 not in params) or (u2 not in params):
                continue
            p1 = params[u1]; p2 = params[u2]
            # [Lp, Lnp, w1, w2, w3, w4, sigma, base]
            Lp1, Lnp1, w1_1, w2_1, w3_1, w4_1, sig1, base1 = p1
            Lp2, Lnp2, w1_2, w2_2, w3_2, w4_2, sig2, base2 = p2

            key = tuple(sorted((u1, u2)))
            bc  = float(bc_lookup.get(key, np.nan))
            sgn = sign_from_bc(bc)

            for cond_idx in paired_indices:
                center_is_pref, periph_is_pref, site = decode_paired_condition(cond_idx)
                ws1 = [w1_1, w2_1, w3_1, w4_1][site]
                ws2 = [w1_2, w2_2, w3_2, w4_2][site]

                Lc1_p   = Lp1 if center_is_pref  else Lnp1
                Lpnp1_p = Lp1 if periph_is_pref  else Lnp1
                Lc2_p   = Lp2 if center_is_pref  else Lnp2
                Lpnp2_p = Lp2 if periph_is_pref  else Lnp2

                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, Lpnp1_p, Lc2_p, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 0))

                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, 0.0, Lc2_p, 0.0, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 1))

                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             0.0, Lpnp1_p, 0.0, Lpnp2_p, ws1, ws2,
                             sig1, base1, sig2, base2, bc, sgn, 2))

        cols = ["unit_i","unit_j","cond_index","site","center_is_pref","periph_is_pref",
                "Lc1","Lpnp1","Lc2","Lpnp2","wsite1","wsite2","sig1","base1","sig2","base2",
                "bhattCoef","sign","stim_kind"]
        pre = pd.DataFrame(rows, columns=cols)

        w0 = 1.0
        f1_base = weight_fraction(w0, pre["wsite1"].to_numpy(), pre["sig1"].to_numpy())
        f2_base = weight_fraction(w0, pre["wsite2"].to_numpy(), pre["sig2"].to_numpy())
        pre["F_base"] = np.sqrt(f1_base * f2_base)
        if WRITE_NONPREFSUPP:
            pre["nonPrefSupp"] = pre["F_base"]

        # Drive proxy (no denom gain; just baseline)
        num1 = w0*pre["Lc1"].to_numpy() + pre["wsite1"].to_numpy()*pre["Lpnp1"].to_numpy()
        num2 = w0*pre["Lc2"].to_numpy() + pre["wsite2"].to_numpy()*pre["Lpnp2"].to_numpy()
        den1 = w0 + pre["wsite1"].to_numpy() + pre["sig1"].to_numpy()
        den2 = w0 + pre["wsite2"].to_numpy() + pre["sig2"].to_numpy()
        F_eff1 = np.divide(num1, den1, out=np.zeros_like(num1), where=den1>0)
        F_eff2 = np.divide(num2, den2, out=np.zeros_like(num2), where=den2>0)
        pre["F_base_eff"] = np.sqrt(np.clip(F_eff1*F_eff2, 0.0, 1.0))

        pre["is_paired"]        = (pre["stim_kind"] == 0).astype(int)
        pre["is_single_center"] = (pre["stim_kind"] == 1).astype(int)
        pre["is_single_periph"] = (pre["stim_kind"] == 2).astype(int)

        return pre


    # ---------- Main simulation ----------
    def simulate_trialwise_weightnoise_denom_only_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    ):
        rng = np.random.default_rng(rng_seed)
        pre = precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup)

        Lc1   = pre["Lc1"].to_numpy()
        Lpnp1 = pre["Lpnp1"].to_numpy()
        Lc2   = pre["Lc2"].to_numpy()
        Lpnp2 = pre["Lpnp2"].to_numpy()
        ws1   = pre["wsite1"].to_numpy()
        ws2   = pre["wsite2"].to_numpy()
        sig1  = pre["sig1"].to_numpy()
        sig2  = pre["sig2"].to_numpy()
        b1    = pre["base1"].to_numpy()
        b2    = pre["base2"].to_numpy()
        sgn   = pre["sign"].to_numpy().astype(float)
        bc    = pre["bhattCoef"].to_numpy()

        n_rows = len(pre)
        w0 = 1.0

        if USE_BC_SCALED_RHO:
            drive = pre["F_base_eff"].to_numpy()
            lo, hi = np.nanpercentile(drive, [5, 95])
            span = max(1e-9, hi - lo)
            drive_norm = np.clip((drive - lo) / span, 0.0, 1.0) ** RHO_DRIVE_EXP
            bc01 = np.clip(bc, 0.0, 1.0)
            rho_eff = RHO_MIN + (RHO_MAX - RHO_MIN) * (bc01 * drive_norm)
        else:
            rho_eff = np.full(n_rows, RHO_SHARED, dtype=float)

        n     = np.zeros(n_rows, dtype=float)
        m1    = np.zeros(n_rows, dtype=float)
        m2    = np.zeros(n_rows, dtype=float)
        M2_1  = np.zeros(n_rows, dtype=float)
        M2_2  = np.zeros(n_rows, dtype=float)
        C12   = np.zeros(n_rows, dtype=float)

        for _ in range(n_trials):
            # ---- GLOBAL denom gain (shared across all rows) ----
            g_den = rng.gamma(DENOM_GAIN_SHAPE, 1.0 / DENOM_GAIN_SHAPE) if APPLY_DENOM_GAIN else 1.0

            # baseline noise (optional)
            b1_t, b2_t = add_spont_noise_to_base(
                b=b1, rng=rng, n_rows=n_rows,
                rho_base=RHO_BASE_SHARED,
                sigma_base=SIGMA_BASE,
                sgn=sgn, signed=BASE_SIGNED_BY_BC,
                mode=BASE_NOISE_MODE,
                base_clip=BASE_NOISE_CLIP
            )

            # rates: denom-only gain (NO numerator gain)
            r1_base = rate_from_weights_denom_only_gain(Lc1, Lpnp1, w0, ws1, sig1, b1_t, g_den=g_den)
            r2_base = rate_from_weights_denom_only_gain(Lc2, Lpnp2, w0, ws2, sig2, b2_t, g_den=g_den)

            # optional latent gain (usually OFF)
            if SIGMA_GAIN > 0:
                eta1 = rng.normal(size=n_rows)
                eta2 = rng.normal(size=n_rows)
                xi   = rng.normal(size=n_rows)

                s1 = np.sqrt(1.0 - rho_eff)
                s2 = np.sqrt(rho_eff)

                lin1 = s1 * eta1 + s2 * xi
                lin2 = s1 * eta2 + sgn * s2 * xi

                gain1 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin1)
                gain2 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin2)

                r1_t = r1_base * gain1
                r2_t = r2_base * gain2
            else:
                r1_t = r1_base
                r2_t = r2_base

            lam1 = np.clip(r1_t * DT_SECONDS, 0, None)
            lam2 = np.clip(r2_t * DT_SECONDS, 0, None)
            k1 = rng.poisson(lam1)
            k2 = rng.poisson(lam2)

            n, m1, m2, M2_1, M2_2, C12 = online_update(k1, k2, n, m1, m2, M2_1, M2_2, C12)

        dof   = np.maximum(1.0, n - 1.0)
        var1  = M2_1 / dof
        var2  = M2_2 / dof
        cov12 = C12 / dof
        denom = np.sqrt(np.clip(var1 * var2, 1e-12, None))
        r_sc  = np.where(denom > 0, cov12 / denom, np.nan)

        out = pre.copy()
        out["mean1"]   = m1
        out["mean2"]   = m2
        out["var1"]    = var1
        out["var2"]    = var2
        out["cov12"]   = cov12
        out["r_sc"]    = r_sc
        out["trials"]  = n.astype(int)
        out["rho_eff"] = rho_eff

        # bookkeeping
        out["denom_gain_shape"] = DENOM_GAIN_SHAPE
        out["apply_denom_gain"] = int(APPLY_DENOM_GAIN)
        out["sigma_gain"]       = SIGMA_GAIN
        out["sigma_base"]       = SIGMA_BASE
        out["rho_base_shared"]  = RHO_BASE_SHARED

        # MATLAB fallback field
        if WRITE_NONPREFSUPP and ("nonPrefSupp" not in out.columns):
            out["nonPrefSupp"] = out["F_base"]

        return out


    # ===================== RUN =====================
    full_df = simulate_trialwise_weightnoise_denom_only_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    )

    print(full_df[["stim_kind", "F_base", "bhattCoef", "r_sc", "mean1", "mean2"]].head())

    # ===================== SAVE (MAT) =====================
    mat_path = "trialwise_weightnoiseDenom_with_singles.mat"
    numeric_cols = {
        k: full_df[k].to_numpy()
        for k in full_df.columns
        if np.issubdtype(full_df[k].dtype, np.number)
    }
    savemat(mat_path, numeric_cols)
    print(f"Saved {mat_path}")

    # ===================== QUICK VISUAL (optional) =====================
    try:
        m = (full_df["stim_kind"] == 0)  # paired only
        x = full_df.loc[m, "F_base"].to_numpy()
        y = full_df.loc[m, "bhattCoef"].to_numpy()
        z = full_df.loc[m, "r_sc"].to_numpy()
        ms = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)

        vlim = 0.15
        plt.figure(figsize=(6.6, 5))
        sc = plt.scatter(x[ms], y[ms], c=z[ms], s=10, cmap="coolwarm", vmin=-vlim, vmax=+vlim)
        plt.colorbar(sc, label="r_sc")
        plt.xlabel("F_base (geom. mean weight fraction)")
        plt.ylabel("Bhattacharyya coefficient (BC)")
        plt.title("Paired only: r_sc over F × BC (DENOM-only gain; NO numerator gain)")
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show()
    except Exception as e:
        print("Quick plot skipped:", e)


# this version adds sigma noise and has that fluctuate globally and nothing else
# ============================================================
# Trialwise SIGMA/ib modulation (NO weight modulation)
# R = (w0*Lc + ws*Lpnp + Lb*ib) / (w0 + ws + ib)
# - One global ib gain per trial, applied to every neuron/pair/row
# - Baseline Lb set near 0 by default (override), and enters via Lb*ib
# - Optional: keep latent gain OFF to isolate sigma/ib-driven correlations
# ============================================================
for filler in range(1):
    import numpy as np
    import pandas as pd
    from scipy.io import savemat
    import matplotlib.pyplot as plt

    # ================= USER SETTINGS =================
    N_TRIALS      = 10000   # trials per (pair × condition or single)
    DT_SECONDS    = 0.25    # 250 ms counting window
    RNG_SEED      = 123

    # ---- Trialwise GLOBAL ib (sigma) gain: one scalar per trial ----
    # Gamma(shape, scale=1/shape): mean=1, var=1/shape
    IBGAIN_SHAPE   = 3.0     # smaller => stronger fluctuations => stronger correlations
    APPLY_IB_GAIN  = True

    # ---- Baseline in the NEW equation (multiplies ib in numerator) ----
    # If you want baseline ~0, keep this tiny (or 0).
    LB_OVERRIDE    = 1e-6
    USE_FIT_BASE   = False   # if True, use fitted "base" as Lb (usually keep False)

    # ---- OPTIONAL latent gain (OFF by default) ----
    # Turn ON only if you want additional shared gain structure beyond ib modulation.
    SIGMA_GAIN     = 0.0     # log-SD of per-trial rate gain; keep 0.0 to isolate ib effects
    KAPPA          = 0.0     # if SIGMA_GAIN>0: 0 => gain modulates rate; 1 => cancels in denom
    USE_BC_SCALED_RHO = False
    RHO_SHARED     = 0.0
    RHO_MIN        = 0.20
    RHO_MAX        = 0.80
    RHO_DRIVE_EXP  = 1.5
    BC_THRESH      = 0.65    # ≥thr => positive shared sign; <thr => negative

    # ---- Plot settings ----
    VCOLOR_LIM     = 0.15

    # ---------- helpers ----------
    def decode_paired_condition(i):
        """Return center_is_pref (bool), periph_is_pref (bool), site_idx (0..3)."""
        if 10 <= i <= 17:
            site = (i - 10) // 2
            periph_is_pref = (i - 10) % 2 == 1
            center_is_pref = False
        elif 19 <= i <= 26:
            site = (i - 19) // 2
            periph_is_pref = (i - 19) % 2 == 1
            center_is_pref = True
        else:
            raise ValueError(f"Condition {i} is not paired (10..17 or 19..26).")
        return center_is_pref, periph_is_pref, site


    def sign_from_bc(bc, thresh=BC_THRESH):
        return +1.0 if (np.isfinite(bc) and bc >= thresh) else -1.0


    def rate_from_weights_ib(Lc, Lpnp, w0, ws, ib, Lb):
        """
        NEW response equation:
          R = (w0*Lc + ws*Lpnp + Lb*ib) / (w0 + ws + ib)
        All inputs can be vectors.
        """
        ib_safe = np.maximum(1e-9, ib)
        den = (w0 + ws + ib_safe)
        num = (w0 * Lc) + (ws * Lpnp) + (Lb * ib_safe)
        return np.divide(num, den, out=np.zeros_like(num, dtype=float), where=den > 0)


    def weight_fraction_ib(w0, ws, ib):
        """
        Weight fraction with ib in the denominator:
          F = (w0 + ws) / (w0 + ws + ib)
        """
        den = w0 + ws + np.maximum(1e-9, ib)
        return np.clip((w0 + ws) / den, 0.0, 1.0)


    # Online (vectorized) update for means, variances, covariance
    def online_update(x, y, n, mean1, mean2, M2_1, M2_2, C12):
        """
        x,y are vectors (one sample per row)
        Updates Welford accumulators.
        """
        n_new = n + 1.0
        dx = x - mean1
        mean1_new = mean1 + dx / n_new
        dy = y - mean2
        mean2_new = mean2 + dy / n_new
        M2_1_new = M2_1 + dx * (x - mean1_new)
        M2_2_new = M2_2 + dy * (y - mean2_new)
        C12_new  = C12  + dx * (y - mean2_new)
        return n_new, mean1_new, mean2_new, M2_1_new, M2_2_new, C12_new


    # ---------- Precompute table (paired + derived singles) ----------
    def precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup):
        """
        Build a long design matrix with 3 rows per 'paired' condition:
          stim_kind = 0: paired (center + that periphery site shown together)
          stim_kind = 1: single_center (center shown, periphery blank -> Lpnp=0)
          stim_kind = 2: single_periph (that periphery site shown, center blank -> Lc=0)
        """
        params = {unit: rfWeightNeuronFits[i] for i, unit in enumerate(masterGoodUnits)}
        paired_indices = list(range(10, 18)) + list(range(19, 27))
        rows = []

        for (u1, u2) in all_combs:
            if (u1 not in params) or (u2 not in params):
                continue
            p1 = params[u1]; p2 = params[u2]
            # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
            Lp1, Lnp1, w1_1, w2_1, w3_1, w4_1, ib1, base1 = p1
            Lp2, Lnp2, w1_2, w2_2, w3_2, w4_2, ib2, base2 = p2

            key = tuple(sorted((u1, u2)))
            bc  = float(bc_lookup.get(key, np.nan))
            sgn = sign_from_bc(bc)

            for cond_idx in paired_indices:
                center_is_pref, periph_is_pref, site = decode_paired_condition(cond_idx)
                ws1 = [w1_1, w2_1, w3_1, w4_1][site]
                ws2 = [w1_2, w2_2, w3_2, w4_2][site]

                # L's for the paired condition
                Lc1_p   = Lp1 if center_is_pref  else Lnp1
                Lpnp1_p = Lp1 if periph_is_pref  else Lnp1
                Lc2_p   = Lp2 if center_is_pref  else Lnp2
                Lpnp2_p = Lp2 if periph_is_pref  else Lnp2

                # (0) paired
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, Lpnp1_p, Lc2_p, Lpnp2_p, ws1, ws2,
                             ib1, base1, ib2, base2, bc, sgn, 0))

                # (1) single_center
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, 0.0, Lc2_p, 0.0, ws1, ws2,
                             ib1, base1, ib2, base2, bc, sgn, 1))

                # (2) single_periph
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             0.0, Lpnp1_p, 0.0, Lpnp2_p, ws1, ws2,
                             ib1, base1, ib2, base2, bc, sgn, 2))

        cols = ["unit_i","unit_j","cond_index","site","center_is_pref","periph_is_pref",
                "Lc1","Lpnp1","Lc2","Lpnp2","wsite1","wsite2","ib1","base1","ib2","base2",
                "bhattCoef","sign","stim_kind"]
        pre = pd.DataFrame(rows, columns=cols)

        # Baseline F metric using ib in denom (no noise)
        w0 = 1.0
        f1_base = weight_fraction_ib(w0, pre["wsite1"].to_numpy(), pre["ib1"].to_numpy())
        f2_base = weight_fraction_ib(w0, pre["wsite2"].to_numpy(), pre["ib2"].to_numpy())
        pre["F_base_ib"] = np.sqrt(f1_base * f2_base)
        pre["F_base"] = pre["F_base_ib"]
        pre["nonPrefSupp"] = pre["F_base_ib"]

        # Drive proxy (rate-ish) for optional rho scaling
        # Using NEW equation with Lb ~ 0 (drive dominated by stimulus terms)
        Lb_tmp = np.full(len(pre), LB_OVERRIDE, dtype=float)
        r1_drive = rate_from_weights_ib(pre["Lc1"].to_numpy(), pre["Lpnp1"].to_numpy(),
                                        w0, pre["wsite1"].to_numpy(),
                                        pre["ib1"].to_numpy(), Lb_tmp)
        r2_drive = rate_from_weights_ib(pre["Lc2"].to_numpy(), pre["Lpnp2"].to_numpy(),
                                        w0, pre["wsite2"].to_numpy(),
                                        pre["ib2"].to_numpy(), Lb_tmp)
        pre["drive_proxy"] = np.sqrt(np.clip(r1_drive * r2_drive, 0.0, None))

        pre["is_paired"]        = (pre["stim_kind"] == 0).astype(int)
        pre["is_single_center"] = (pre["stim_kind"] == 1).astype(int)
        pre["is_single_periph"] = (pre["stim_kind"] == 2).astype(int)

        return pre


    # ---------- Main simulation (ib-only modulation) ----------
    def simulate_trialwise_ib_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    ):
        rng = np.random.default_rng(rng_seed)
        pre = precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup)

        # numpy arrays for speed
        Lc1   = pre["Lc1"].to_numpy()
        Lpnp1 = pre["Lpnp1"].to_numpy()
        Lc2   = pre["Lc2"].to_numpy()
        Lpnp2 = pre["Lpnp2"].to_numpy()
        ws1   = pre["wsite1"].to_numpy()
        ws2   = pre["wsite2"].to_numpy()
        ib1   = pre["ib1"].to_numpy()
        ib2   = pre["ib2"].to_numpy()
        base1 = pre["base1"].to_numpy()
        base2 = pre["base2"].to_numpy()
        bc    = pre["bhattCoef"].to_numpy()
        sgn   = pre["sign"].to_numpy().astype(float)

        n_rows = len(pre)
        w0 = 1.0  # fixed, no modulation

        # Choose Lb for each neuron (baseline multiplier on ib in numerator)
        if USE_FIT_BASE:
            Lb1 = np.clip(base1, 0.0, None)
            Lb2 = np.clip(base2, 0.0, None)
        else:
            Lb1 = np.full(n_rows, LB_OVERRIDE, dtype=float)
            Lb2 = np.full(n_rows, LB_OVERRIDE, dtype=float)

        # Optional rho scaling if latent gain enabled
        if USE_BC_SCALED_RHO:
            drive = pre["drive_proxy"].to_numpy()
            lo, hi = np.nanpercentile(drive, [5, 95])
            span = max(1e-9, hi - lo)
            drive_norm = np.clip((drive - lo) / span, 0.0, 1.0) ** RHO_DRIVE_EXP
            bc01 = np.clip(bc, 0.0, 1.0)
            rho_eff = RHO_MIN + (RHO_MAX - RHO_MIN) * (bc01 * drive_norm)
        else:
            rho_eff = np.full(n_rows, RHO_SHARED, dtype=float)

        # Online accumulators per row
        n     = np.zeros(n_rows, dtype=float)
        m1    = np.zeros(n_rows, dtype=float)
        m2    = np.zeros(n_rows, dtype=float)
        M2_1  = np.zeros(n_rows, dtype=float)
        M2_2  = np.zeros(n_rows, dtype=float)
        C12   = np.zeros(n_rows, dtype=float)

        for _ in range(n_trials):
            # ---- GLOBAL per-trial ib gain (shared across all neurons/rows) ----
            g_ib = rng.gamma(IBGAIN_SHAPE, 1.0 / IBGAIN_SHAPE) if APPLY_IB_GAIN else 1.0

            ib1_t = np.maximum(1e-9, ib1 * g_ib)
            ib2_t = np.maximum(1e-9, ib2 * g_ib)

            # Base rates from NEW equation
            r1_base = rate_from_weights_ib(Lc1, Lpnp1, w0, ws1, ib1_t, Lb1)
            r2_base = rate_from_weights_ib(Lc2, Lpnp2, w0, ws2, ib2_t, Lb2)

            # ---- OPTIONAL latent gain ----
            if SIGMA_GAIN > 0:
                eta1 = rng.normal(size=n_rows)
                eta2 = rng.normal(size=n_rows)
                xi   = rng.normal(size=n_rows)

                s1 = np.sqrt(1.0 - rho_eff)
                s2 = np.sqrt(rho_eff)

                lin1 = s1 * eta1 + s2 * xi
                lin2 = s1 * eta2 + sgn * s2 * xi

                gain1 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin1)
                gain2 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin2)

                r1_t = r1_base * gain1
                r2_t = r2_base * gain2
            else:
                r1_t = r1_base
                r2_t = r2_base

            lam1 = np.clip(r1_t * DT_SECONDS, 0, None)
            lam2 = np.clip(r2_t * DT_SECONDS, 0, None)
            k1 = rng.poisson(lam1)
            k2 = rng.poisson(lam2)

            n, m1, m2, M2_1, M2_2, C12 = online_update(k1, k2, n, m1, m2, M2_1, M2_2, C12)

        # final stats
        dof   = np.maximum(1.0, n - 1.0)
        var1  = M2_1 / dof
        var2  = M2_2 / dof
        cov12 = C12 / dof
        denom = np.sqrt(np.clip(var1 * var2, 1e-12, None))
        r_sc  = np.where(denom > 0, cov12 / denom, np.nan)

        out = pre.copy()
        out["mean1"]   = m1
        out["mean2"]   = m2
        out["var1"]    = var1
        out["var2"]    = var2
        out["cov12"]   = cov12
        out["r_sc"]    = r_sc
        out["trials"]  = n.astype(int)
        out["rho_eff"] = rho_eff

        # store settings
        out["ibgain_shape"] = IBGAIN_SHAPE
        out["Lb_override"]  = LB_OVERRIDE
        out["used_fit_base"] = int(USE_FIT_BASE)
        out["sigma_gain"] = SIGMA_GAIN
        out["F_base"] = out["F_base_ib"]
        out["nonPrefSupp"] = out["F_base_ib"]

        return out


    # ========== PREREQS (must already exist in your session) ==========
    # all_combs: list of (unit_i, unit_j)
    # masterGoodUnits: np.array of unit IDs (strings), aligned with rfWeightNeuronFits
    # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
    # bc_lookup[tuple(sorted((u1,u2)))] = BC float

    # ===================== RUN =====================
    full_df = simulate_trialwise_ib_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    )

    print(full_df[["stim_kind", "F_base_ib", "bhattCoef", "r_sc", "mean1", "mean2"]].head())

    # ===================== SAVE (MAT for MATLAB meshes) =====================
    mat_path = "trialwise_ibnoise_with_singles.mat"
    numeric_cols = {
        k: full_df[k].to_numpy()
        for k in full_df.columns
        if np.issubdtype(full_df[k].dtype, np.number)
    }
    savemat(mat_path, numeric_cols)
    print(f"Saved {mat_path}")

    # ===================== QUICK VISUAL (optional) =====================
    try:
        m = (full_df["stim_kind"] == 0)  # paired only
        x = full_df.loc[m, "F_base_ib"].to_numpy()
        y = full_df.loc[m, "bhattCoef"].to_numpy()
        z = full_df.loc[m, "r_sc"].to_numpy()
        ms = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)

        plt.figure(figsize=(6.6, 5))
        sc = plt.scatter(x[ms], y[ms], c=z[ms], s=10, cmap="coolwarm",
                         vmin=-VCOLOR_LIM, vmax=+VCOLOR_LIM)
        plt.colorbar(sc, label="r_sc")
        plt.xlabel("F_base_ib (geom. mean weight fraction with ib)")
        plt.ylabel("Bhattacharyya coefficient (BC)")
        plt.title("Paired only: r_sc over F × BC (trialwise global ib gain)")
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show()
    except Exception as e:
        print("Quick plot skipped:", e)


# this version adds noise to the L terms that fluctuate globally.
# ============================================================
# Trialwise L-noise (NO sigma modulation)
# R = (w0*(Lc*gL) + ws*(Lpnp*gL) + Lb*ib) / (w0 + ws + ib)
# - One global L gain per trial, shared across ALL neurons/rows
# - Sigma/ib stays fixed at fitted value (no modulation)
# - Saves: trialwise_Lnoise_with_singles.mat
# ============================================================
for filler in range(1):
    import numpy as np
    import pandas as pd
    from scipy.io import savemat
    import matplotlib.pyplot as plt

    # ================= USER SETTINGS =================
    N_TRIALS      = 10000   # trials per (pair × condition or single)
    DT_SECONDS    = 0.25    # 250 ms counting window
    RNG_SEED      = 123

    # ---- Trialwise GLOBAL gain on L terms (numerator) ----
    # Gamma(shape, scale=1/shape): mean=1, var=1/shape
    LGAIN_SHAPE   = 3.0     # smaller => stronger fluctuations => stronger correlations
    APPLY_L_GAIN  = True

    # ---- Baseline in the equation (multiplies ib in numerator) ----
    # If you want baseline ~0, keep this tiny (or 0).
    LB_OVERRIDE    = 1e-6
    USE_FIT_BASE   = False   # if True, use fitted "base" as Lb (usually keep False)

    # ---- OPTIONAL latent gain (OFF by default) ----
    # Independent additional shared gain, if you want it (usually keep OFF)
    SIGMA_GAIN     = 0.0
    KAPPA          = 0.0
    USE_BC_SCALED_RHO = False
    RHO_SHARED     = 0.0
    RHO_MIN        = 0.20
    RHO_MAX        = 0.80
    RHO_DRIVE_EXP  = 1.5
    BC_THRESH      = 0.65    # ≥thr => positive shared sign; <thr => negative

    # ---- Plot settings ----
    VCOLOR_LIM     = 0.15

    # ---------- helpers ----------
    def decode_paired_condition(i):
        """Return center_is_pref (bool), periph_is_pref (bool), site_idx (0..3)."""
        if 10 <= i <= 17:
            site = (i - 10) // 2
            periph_is_pref = (i - 10) % 2 == 1
            center_is_pref = False
        elif 19 <= i <= 26:
            site = (i - 19) // 2
            periph_is_pref = (i - 19) % 2 == 1
            center_is_pref = True
        else:
            raise ValueError(f"Condition {i} is not paired (10..17 or 19..26).")
        return center_is_pref, periph_is_pref, site


    def sign_from_bc(bc, thresh=BC_THRESH):
        return +1.0 if (np.isfinite(bc) and bc >= thresh) else -1.0


    def rate_from_weights_Lnoise(Lc, Lpnp, w0, ws, ib, Lb, gL):
        """
        Equation with global L-gain on numerator stimulus terms only:
          R = (w0*(Lc*gL) + ws*(Lpnp*gL) + Lb*ib) / (w0 + ws + ib)
        """
        ib_safe = np.maximum(1e-9, ib)
        den = (w0 + ws + ib_safe)
        num = (w0 * (Lc * gL)) + (ws * (Lpnp * gL)) + (Lb * ib_safe)
        return np.divide(num, den, out=np.zeros_like(num, dtype=float), where=den > 0)


    def weight_fraction_ib(w0, ws, ib):
        """
        Weight fraction with ib in the denominator:
          F = (w0 + ws) / (w0 + ws + ib)
        """
        den = w0 + ws + np.maximum(1e-9, ib)
        return np.clip((w0 + ws) / den, 0.0, 1.0)


    def online_update(x, y, n, mean1, mean2, M2_1, M2_2, C12):
        """Vectorized Welford update."""
        n_new = n + 1.0
        dx = x - mean1
        mean1_new = mean1 + dx / n_new
        dy = y - mean2
        mean2_new = mean2 + dy / n_new
        M2_1_new = M2_1 + dx * (x - mean1_new)
        M2_2_new = M2_2 + dy * (y - mean2_new)
        C12_new  = C12  + dx * (y - mean2_new)
        return n_new, mean1_new, mean2_new, M2_1_new, M2_2_new, C12_new


    # ---------- Precompute table (paired + derived singles) ----------
    def precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup):
        """
        Build a long design matrix with 3 rows per 'paired' condition:
          stim_kind = 0: paired (center + that periphery site shown together)
          stim_kind = 1: single_center (center shown, periphery blank -> Lpnp=0)
          stim_kind = 2: single_periph (that periphery site shown, center blank -> Lc=0)
        """
        params = {unit: rfWeightNeuronFits[i] for i, unit in enumerate(masterGoodUnits)}
        paired_indices = list(range(10, 18)) + list(range(19, 27))
        rows = []

        for (u1, u2) in all_combs:
            if (u1 not in params) or (u2 not in params):
                continue
            p1 = params[u1]; p2 = params[u2]
            # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
            Lp1, Lnp1, w1_1, w2_1, w3_1, w4_1, ib1, base1 = p1
            Lp2, Lnp2, w1_2, w2_2, w3_2, w4_2, ib2, base2 = p2

            key = tuple(sorted((u1, u2)))
            bc  = float(bc_lookup.get(key, np.nan))
            sgn = sign_from_bc(bc)

            for cond_idx in paired_indices:
                center_is_pref, periph_is_pref, site = decode_paired_condition(cond_idx)
                ws1 = [w1_1, w2_1, w3_1, w4_1][site]
                ws2 = [w1_2, w2_2, w3_2, w4_2][site]

                # L's for the paired condition
                Lc1_p   = Lp1 if center_is_pref  else Lnp1
                Lpnp1_p = Lp1 if periph_is_pref  else Lnp1
                Lc2_p   = Lp2 if center_is_pref  else Lnp2
                Lpnp2_p = Lp2 if periph_is_pref  else Lnp2

                # (0) paired
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, Lpnp1_p, Lc2_p, Lpnp2_p, ws1, ws2,
                             ib1, base1, ib2, base2, bc, sgn, 0))

                # (1) single_center
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             Lc1_p, 0.0, Lc2_p, 0.0, ws1, ws2,
                             ib1, base1, ib2, base2, bc, sgn, 1))

                # (2) single_periph
                rows.append((u1, u2, cond_idx, site, center_is_pref, periph_is_pref,
                             0.0, Lpnp1_p, 0.0, Lpnp2_p, ws1, ws2,
                             ib1, base1, ib2, base2, bc, sgn, 2))

        cols = ["unit_i","unit_j","cond_index","site","center_is_pref","periph_is_pref",
                "Lc1","Lpnp1","Lc2","Lpnp2","wsite1","wsite2","ib1","base1","ib2","base2",
                "bhattCoef","sign","stim_kind"]
        pre = pd.DataFrame(rows, columns=cols)

        # Baseline F metric using ib in denom (no noise)
        w0 = 1.0
        f1_base = weight_fraction_ib(w0, pre["wsite1"].to_numpy(), pre["ib1"].to_numpy())
        f2_base = weight_fraction_ib(w0, pre["wsite2"].to_numpy(), pre["ib2"].to_numpy())
        pre["F_base_ib"] = np.sqrt(f1_base * f2_base)
        pre["F_base"] = pre["F_base_ib"]          # MATLAB expects F_base
        pre["nonPrefSupp"] = pre["F_base_ib"]     # MATLAB fallback

        # Drive proxy for optional rho scaling (use gL=1 here)
        if USE_FIT_BASE:
            Lb1 = np.clip(pre["base1"].to_numpy(), 0.0, None)
            Lb2 = np.clip(pre["base2"].to_numpy(), 0.0, None)
        else:
            Lb1 = np.full(len(pre), LB_OVERRIDE, dtype=float)
            Lb2 = np.full(len(pre), LB_OVERRIDE, dtype=float)

        r1_drive = rate_from_weights_Lnoise(pre["Lc1"].to_numpy(), pre["Lpnp1"].to_numpy(),
                                            w0, pre["wsite1"].to_numpy(),
                                            pre["ib1"].to_numpy(), Lb1, gL=1.0)
        r2_drive = rate_from_weights_Lnoise(pre["Lc2"].to_numpy(), pre["Lpnp2"].to_numpy(),
                                            w0, pre["wsite2"].to_numpy(),
                                            pre["ib2"].to_numpy(), Lb2, gL=1.0)
        pre["drive_proxy"] = np.sqrt(np.clip(r1_drive * r2_drive, 0.0, None))

        pre["is_paired"]        = (pre["stim_kind"] == 0).astype(int)
        pre["is_single_center"] = (pre["stim_kind"] == 1).astype(int)
        pre["is_single_periph"] = (pre["stim_kind"] == 2).astype(int)

        return pre


    # ---------- Main simulation (L-only modulation) ----------
    def simulate_trialwise_L_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    ):
        rng = np.random.default_rng(rng_seed)
        pre = precompute_table_with_singles(all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup)

        # numpy arrays for speed
        Lc1   = pre["Lc1"].to_numpy()
        Lpnp1 = pre["Lpnp1"].to_numpy()
        Lc2   = pre["Lc2"].to_numpy()
        Lpnp2 = pre["Lpnp2"].to_numpy()
        ws1   = pre["wsite1"].to_numpy()
        ws2   = pre["wsite2"].to_numpy()
        ib1   = pre["ib1"].to_numpy()    # FIXED sigma/ib (no modulation)
        ib2   = pre["ib2"].to_numpy()
        base1 = pre["base1"].to_numpy()
        base2 = pre["base2"].to_numpy()
        bc    = pre["bhattCoef"].to_numpy()
        sgn   = pre["sign"].to_numpy().astype(float)

        n_rows = len(pre)
        w0 = 1.0  # fixed

        # Choose Lb for each neuron (baseline multiplier on ib in numerator)
        if USE_FIT_BASE:
            Lb1 = np.clip(base1, 0.0, None)
            Lb2 = np.clip(base2, 0.0, None)
        else:
            Lb1 = np.full(n_rows, LB_OVERRIDE, dtype=float)
            Lb2 = np.full(n_rows, LB_OVERRIDE, dtype=float)

        # Optional rho scaling if latent gain enabled
        if USE_BC_SCALED_RHO:
            drive = pre["drive_proxy"].to_numpy()
            lo, hi = np.nanpercentile(drive, [5, 95])
            span = max(1e-9, hi - lo)
            drive_norm = np.clip((drive - lo) / span, 0.0, 1.0) ** RHO_DRIVE_EXP
            bc01 = np.clip(bc, 0.0, 1.0)
            rho_eff = RHO_MIN + (RHO_MAX - RHO_MIN) * (bc01 * drive_norm)
        else:
            rho_eff = np.full(n_rows, RHO_SHARED, dtype=float)

        # Online accumulators per row
        n     = np.zeros(n_rows, dtype=float)
        m1    = np.zeros(n_rows, dtype=float)
        m2    = np.zeros(n_rows, dtype=float)
        M2_1  = np.zeros(n_rows, dtype=float)
        M2_2  = np.zeros(n_rows, dtype=float)
        C12   = np.zeros(n_rows, dtype=float)

        for _ in range(n_trials):
            # ---- GLOBAL per-trial L gain (shared across all neurons/rows) ----
            gL = rng.gamma(LGAIN_SHAPE, 1.0 / LGAIN_SHAPE) if APPLY_L_GAIN else 1.0

            # Base rates from equation with L-noise
            r1_base = rate_from_weights_Lnoise(Lc1, Lpnp1, w0, ws1, ib1, Lb1, gL=gL)
            r2_base = rate_from_weights_Lnoise(Lc2, Lpnp2, w0, ws2, ib2, Lb2, gL=gL)

            # ---- OPTIONAL latent gain ----
            if SIGMA_GAIN > 0:
                eta1 = rng.normal(size=n_rows)
                eta2 = rng.normal(size=n_rows)
                xi   = rng.normal(size=n_rows)

                s1 = np.sqrt(1.0 - rho_eff)
                s2 = np.sqrt(rho_eff)

                lin1 = s1 * eta1 + s2 * xi
                lin2 = s1 * eta2 + sgn * s2 * xi

                gain1 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin1)
                gain2 = np.exp((1.0 - KAPPA) * SIGMA_GAIN * lin2)

                r1_t = r1_base * gain1
                r2_t = r2_base * gain2
            else:
                r1_t = r1_base
                r2_t = r2_base

            lam1 = np.clip(r1_t * DT_SECONDS, 0, None)
            lam2 = np.clip(r2_t * DT_SECONDS, 0, None)
            k1 = rng.poisson(lam1)
            k2 = rng.poisson(lam2)

            n, m1, m2, M2_1, M2_2, C12 = online_update(k1, k2, n, m1, m2, M2_1, M2_2, C12)

        # final stats
        dof   = np.maximum(1.0, n - 1.0)
        var1  = M2_1 / dof
        var2  = M2_2 / dof
        cov12 = C12 / dof
        denom = np.sqrt(np.clip(var1 * var2, 1e-12, None))
        r_sc  = np.where(denom > 0, cov12 / denom, np.nan)

        out = pre.copy()
        out["mean1"]   = m1
        out["mean2"]   = m2
        out["var1"]    = var1
        out["var2"]    = var2
        out["cov12"]   = cov12
        out["r_sc"]    = r_sc
        out["trials"]  = n.astype(int)
        out["rho_eff"] = rho_eff

        # store settings
        out["lgain_shape"] = LGAIN_SHAPE
        out["apply_l_gain"] = int(APPLY_L_GAIN)
        out["Lb_override"]  = LB_OVERRIDE
        out["used_fit_base"] = int(USE_FIT_BASE)
        out["sigma_gain"] = SIGMA_GAIN

        # MATLAB compatibility fields
        out["F_base"] = out["F_base_ib"]
        out["nonPrefSupp"] = out["F_base_ib"]

        return out


    # ========== PREREQS (must already exist in your session) ==========
    # all_combs: list of (unit_i, unit_j)
    # masterGoodUnits: np.array of unit IDs (strings), aligned with rfWeightNeuronFits
    # rfWeightNeuronFits[i] = [Lp, Lnp, w1, w2, w3, w4, sigma, base]
    # bc_lookup[tuple(sorted((u1,u2)))] = BC float

    # ===================== RUN =====================
    full_df = simulate_trialwise_L_noise_with_singles(
        all_combs, masterGoodUnits, rfWeightNeuronFits, bc_lookup,
        n_trials=N_TRIALS, rng_seed=RNG_SEED
    )

    print(full_df[["stim_kind", "F_base", "bhattCoef", "r_sc", "mean1", "mean2"]].head())

    # ===================== SAVE (MAT for MATLAB meshes) =====================
    mat_path = "trialwise_Lnoise_with_singles.mat"
    numeric_cols = {
        k: full_df[k].to_numpy()
        for k in full_df.columns
        if np.issubdtype(full_df[k].dtype, np.number)
    }
    savemat(mat_path, numeric_cols)
    print(f"Saved {mat_path}")

    # ===================== QUICK VISUAL (optional) =====================
    try:
        m = (full_df["stim_kind"] == 0)  # paired only
        x = full_df.loc[m, "F_base"].to_numpy()
        y = full_df.loc[m, "bhattCoef"].to_numpy()
        z = full_df.loc[m, "r_sc"].to_numpy()
        ms = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)

        plt.figure(figsize=(6.6, 5))
        sc = plt.scatter(x[ms], y[ms], c=z[ms], s=10, cmap="coolwarm",
                         vmin=-VCOLOR_LIM, vmax=+VCOLOR_LIM)
        plt.colorbar(sc, label="r_sc")
        plt.xlabel("F_base (geom. mean weight fraction)")
        plt.ylabel("Bhattacharyya coefficient (BC)")
        plt.title("Paired only: r_sc over F × BC (trialwise global L gain)")
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show()
    except Exception as e:
        print("Quick plot skipped:", e)



