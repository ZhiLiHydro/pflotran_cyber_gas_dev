#!/usr/bin/env python
"""Independent numpy oracle for the ML fracture-surrogate regression tests (no torch, no PFLOTRAN code).

For a single-cell (0-D) deck the surrogate's inputs are constant, so the porosity (and, for node_min_v1,
the four mineral states) can be recomputed from the weight file alone with the documented contract:

  stage 0 (sinh rate):  rate = 1e-5*sinh(trunk(x)*y_scale + y_mean); rate <- min(rate,0), >= -3;
                        phi <- max(phi + rate*dt, 0.02); cells with phi <= 0.0201 are skipped
  stage 1 (node_min_v1): rates from the bounded output stage; explicit Euler on the four states with the
                        pore-volume budget and the Ca/Si clamp; phi = phi0 - vf_csh - vf_m4ac - vf_other
  permeability = k0 * (phi/phi0)^3

Usage:
  python ml_surrogate_oracle.py <deck.in>              # after the deck has been run: checks <deck>.regression
                                                       # and the per-step diagnostic dump against the oracle
  python ml_surrogate_oracle.py <deck.in> -e pflotran  # runs pflotran first (~0.3 s), then checks
  python ml_surrogate_oracle.py --gold <deck.in>       # checks the committed <deck>.regression.gold instead
Exit status 0 = PASS. Tolerances: 1e-9 relative with a 1e-10 absolute floor (documented in the handoff).
"""
import argparse, csv, os, re, subprocess, sys
import numpy as np, h5py

# ------------------------------------------------------------------ weight file + forward pass (both stages)
def load(wf):
    with h5py.File(wf, "r") as f:
        W = [f[f"/layers/layer_{i:02d}/W"][:] for i in range(4)]; B = [f[f"/layers/layer_{i:02d}/b"][:] for i in range(4)]
        S = dict(Xm=f["/scalers/X_mean"][:], Xs=f["/scalers/X_scale"][:], logm=f["/transforms/log10_mask"][:].astype(bool),
                 lo=f["/transforms/clip_lo"][:], hi=f["/transforms/clip_hi"][:], attrs={k: v for k, v in f.attrs.items()})
        S["stage"] = int(f.attrs.get("output_stage_id", 0))
        if S["stage"] == 0: S["ym"] = f["/scalers/y_mean"][0]; S["ys"] = f["/scalers/y_scale"][0]
    return W, B, S

def trunk(W, B, S, raw):
    x = np.clip(raw, S["lo"], S["hi"]); x = np.where(S["logm"], np.log10(np.maximum(x, 1e-300)), x)
    h = (x - S["Xm"]) / S["Xs"]
    for i in range(3): h = np.maximum(W[i] @ h + B[i], 0.0)
    return W[3] @ h + B[3]

def softplus(v): return np.log1p(np.exp(-abs(v))) + max(v, 0.0)
def sigmoid(v): return 1.0 / (1.0 + np.exp(-v)) if v >= 0 else np.exp(v) / (1.0 + np.exp(v))

def rates_node_min(W, B, S, raw):
    a = S["attrs"]; o = trunk(W, B, S, raw); phi, vf, meancs = raw[0], raw[7], raw[8]; s_state = vf * meancs
    ss = min(max((phi - a["phi_res"]) / a["ramp_width"], 0.0), 1.0)
    ramp = ss if int(a.get("ramp_form_id", 0)) == 1 else (3 * ss**2 - 2 * ss**3) ** 2
    R = a["rate_scale"]; r_csh, r_m4, r_oth = R * softplus(o[0]) * ramp, R * softplus(o[2]) * ramp, R * softplus(o[3]) * ramp
    cs_new = a["cs_min"] + (a["cs_max"] - a["cs_min"]) * sigmoid(o[1]); tot = r_csh + r_m4 + r_oth
    if tot > a["rate_max"]: g = a["rate_max"] / tot; r_csh, r_m4, r_oth = r_csh * g, r_m4 * g, r_oth * g
    u = np.tanh(o[4]); r_x = a["k_exch"] * (max(u, 0) * (a["cs_max"] * vf - s_state) - max(-u, 0) * (s_state - a["cs_min"] * vf))
    return np.array([r_csh, cs_new * r_csh + r_x, r_m4, r_oth])

# ------------------------------------------------------------------ host update rules
PHI0, PHI_FLOOR, PHI_SEALED, RATE_MAX = 0.99, 0.02, 0.0201, 3.0

def step_phi(W, B, S, phi, ap, tot, dt):
    if phi <= PHI_SEALED: return phi
    z = trunk(W, B, S, np.array([phi, ap] + list(tot)))[0] * S["ys"] + S["ym"]
    r = max(min(1e-5 * np.sinh(z), 0.0), -RATE_MAX)
    return max(phi + r * dt, PHI_FLOOR)

def step_min(W, B, S, x, ap, tot, dt):
    a = S["attrs"]; vf, s, m4, oth = x; phi = a["phi0"] - vf - m4 - oth
    if phi <= PHI_SEALED: return x
    meancs = s / vf if vf > 1e-6 else 0.0
    d = rates_node_min(W, B, S, np.array([phi, ap] + list(tot) + [vf, meancs, m4])) * dt
    room = max(phi - a["phi_res"], 0.0); grow = d[0] + d[2] + d[3]
    sc = 0.0 if room <= 1e-9 else (room / grow if grow > room else 1.0)
    x = x + d * sc
    x[1] = min(max(x[1], a["cs_min"] * x[0]), a["cs_max"] * x[0])
    return x

# ------------------------------------------------------------------ file readers
def read_regression(path):
    out = {}; cur = None
    for ln in open(path):
        m = re.match(r"-- (\w+): (.*) --", ln.strip())
        if m: cur = m.group(2); out[cur] = {}; continue
        m = re.match(r"\s*(\S+):\s+([-+0-9.Ee]+)", ln)
        if m and cur: out[cur][m.group(1)] = float(m.group(2))
    return out

def read_diag(path):
    return list(csv.DictReader(open(path)))

def deck_params(deck):
    s = open(deck).read()
    w = re.search(r"FILENAME\s+\./?(\S+)", s).group(1); k0 = float(re.search(r"PERM_ISO\s+([-\d.dDeE+]+)", s).group(1).replace("d", "e"))
    diag = re.search(r"DIAGNOSTIC_FILE\s+(\S+)", s).group(1)
    return w, k0, diag

# ------------------------------------------------------------------ the check
def check(deck, gold=False, verbose=True):
    d = os.path.dirname(os.path.abspath(deck)); base = os.path.splitext(os.path.basename(deck))[0]
    w, k0, diag = deck_params(deck); W, B, S = load(os.path.join(d, w)); ap = np.sqrt(12.0 * k0)
    rows = read_diag(os.path.join(d, diag + ".0.csv")); reg = read_regression(os.path.join(d, base + (".regression.gold" if gold else ".regression")))
    keys = ["TotalH", "Ca", "Si", "Al", "Mg"]
    tots = np.array([[float(r[k]) for k in keys] for r in rows])
    times = np.array([float(r["time_y"]) for r in rows]); dts = np.diff(np.concatenate([[0.0], times]))
    worst = 0.0; msgs = []
    def rel(a, b): return abs(a - b) / max(abs(b), 1e-10)
    # inputs the surrogate saw must be constant (closed cell, no reactions)
    drift = np.abs(tots - tots[0]).max(0) / np.maximum(np.abs(tots[0]), 1e-30)
    msgs.append(f"totals drift over the run (relative): {drift.max():.2e}")
    if S["stage"] == 0:
        phi = PHI0
        for i, r in enumerate(rows):
            worst = max(worst, rel(float(r["porosity"]), phi))          # the input phi of step i is the oracle state
            phi = step_phi(W, B, S, phi, ap, tots[i], dts[i])
        phi_final = phi
    else:
        x = np.zeros(4)
        for i, r in enumerate(rows):
            a = S["attrs"]; phi_in = a["phi0"] - x[0] - x[2] - x[3]
            meancs = x[1] / x[0] if x[0] > 1e-6 else 0.0
            for got, exp in [(float(r["porosity"]), phi_in), (float(r["vf_CSH_total"]), x[0]), (float(r["CSH_meanCS"]), meancs), (float(r["vf_M4Ac"]), x[2]), (float(r["vf_other"]), x[3])]:
                worst = max(worst, rel(got, exp))
            x = step_min(W, B, S, x, ap, tots[i], dts[i])
        phi_final = S["attrs"]["phi0"] - x[0] - x[2] - x[3]
    por = reg["Porosity"]["1"]; perm = reg["Permeability X"]["1"] if "Permeability X" in reg else [v for k, v in reg.items() if k.startswith("Permeability")][0]["1"]
    worst = max(worst, rel(por, phi_final), rel(perm, k0 * (phi_final / PHI0) ** 3))
    msgs.append(f"{base}: stage {S['stage']}, {len(rows)} surrogate steps, final phi oracle {phi_final:.10f} vs {'gold' if gold else 'run'} {por:.10f}; worst rel err {worst:.2e}")
    if verbose: print("\n".join(msgs))
    return worst

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("decks", nargs="+"); ap.add_argument("-e", "--executable"); ap.add_argument("--gold", action="store_true")
    ap.add_argument("--tol", type=float, default=1e-9); cfg = ap.parse_args(); fail = 0
    for deck in cfg.decks:
        d = os.path.dirname(os.path.abspath(deck))
        if cfg.executable:
            p = subprocess.run([cfg.executable, "-pflotranin", os.path.basename(deck)], cwd=d, capture_output=True, text=True)
            if p.returncode != 0: print(f"{deck}: pflotran failed rc={p.returncode}"); fail += 1; continue
        try:
            worst = check(deck, gold=cfg.gold)
        except Exception as e:
            print(f"{deck}: oracle error: {e}"); fail += 1; continue
        if worst >= cfg.tol: print(f"  FAIL ({worst:.2e} >= {cfg.tol:.0e})"); fail += 1
        else: print("  PASS")
    sys.exit(1 if fail else 0)
