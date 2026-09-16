# ML fracture-surrogate regression tests

Four 300-y katdiss decks (run.1 of the `lhs_1024_chem_new_db_cement_katdiss` sweep) with the
`FRACTURE_ML_SURROGATE` block active, one per shipped weight file:

| test | weight file | output stage |
|---|---|---|
| ml_bthuber | bthuber_weights.h5 | 0 (sinh rate; the deployed July 2026 MLP) |
| ml_node_phi | node_phi_weights.h5 | 0 (neural-ODE phi, drop-in) |
| ml_node_min_lin | node_min_lin_weights.h5 | 1 (node_min_v1, linear ramp; recommended) |
| ml_node_min | node_min_weights.h5 | 1 (node_min_v1, smoothstep2 ramp) |

The REGRESSION block records porosity and permeability of fracture cells 91..116. Gold files were
generated with the branch `heeho/ml-node-min` binary.

Run:  `make ml_surrogate` (from regression_tests/), or `make test` (the suite is part of `standard`).
Forward-pass golden vectors (no PFLOTRAN run): `make test_ml_surrogate` from src/pflotran/.

The weight files and test vectors are copies of `node_development/export/` and `chuber_handoff/`
in the ML repository; regenerate them there (export_node.py) and copy, then `--update` the gold files
only when a changed result is intended.
