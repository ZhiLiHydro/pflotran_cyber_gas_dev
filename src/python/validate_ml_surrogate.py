#!/usr/bin/env python3
"""Python reference forward pass for the NUMO fracture-clogging MLP.

Mirrors MLPForward / MLPGelu in ml_surrogate_aux.F90 using numpy.  Used
to validate the Fortran implementation: |phi_fortran - phi_python| should
be below 1e-10 in float64.

Author: Heeho Park
Date: 2026-05-18
"""

import argparse
import sys

import h5py
import numpy as np


PHI_MIN = 0.02
PHI_MAX = 0.99


def gelu(x: np.ndarray) -> np.ndarray:
    c1 = np.sqrt(2.0 / np.pi)
    c2 = 0.044715
    return 0.5 * x * (1.0 + np.tanh(c1 * (x + c2 * x ** 3)))


def load_weights(path):
    with h5py.File(path, "r") as f:
        X_mean = f["/scalers/X_mean"][:]
        X_scale = f["/scalers/X_scale"][:]
        y_mean = f["/scalers/y_mean"][:][0]
        y_scale = f["/scalers/y_scale"][:][0]
        layers = []
        for i in range(5):
            W = f[f"/layers/layer_0{i}/W"][:]
            b = f[f"/layers/layer_0{i}/b"][:]
            layers.append((W, b))
    return X_mean, X_scale, y_mean, y_scale, layers


def forward(x_in, X_mean, X_scale, y_mean, y_scale, layers):
    h = (np.asarray(x_in, dtype=np.float64) - X_mean) / X_scale
    for i, (W, b) in enumerate(layers):
        h = W @ h + b
        if i < len(layers) - 1:
            h = gelu(h)
    phi_norm = float(h[0])
    phi = phi_norm * y_scale + y_mean
    return float(np.clip(phi, PHI_MIN, PHI_MAX))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("weights", help="Path to model_mlp_weights.h5")
    args = parser.parse_args()

    weights = load_weights(args.weights)

    cases = [
        (-9.0, -3.0, 500.0, 0.3),
        (-6.0, -2.0, 100.0, 0.1),
        (-12.0, -5.0, 900.0, 0.9),
    ]
    for x in cases:
        phi = forward(x, *weights)
        print(
            "input {:12.4f} {:12.4f} {:12.4f} {:12.4f}  -> phi = {:25.17e}".format(
                *x, phi
            )
        )


if __name__ == "__main__":
    sys.exit(main())
