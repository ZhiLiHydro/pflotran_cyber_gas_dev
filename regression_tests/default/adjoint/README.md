# Discrete adjoint of MODE GIRT_DAE

A deck's `ADJOINT` block (costs, observations and parameters) runs PETSc's
`TSAdjointSolve` after the forward run of `MODE GIRT_DAE`, alone or with
`MODE RICHARDS_DAE` flow as one model, and prints

```
=== ADJOINT SUMMARY ===
Cost J_1 = TerminalConc_CO2(aq)
  J_1(T)              = ...
  dJ_1/dp_1 [MICROBIAL_RATE_CONSTANT[1]] = ...
  max |dJ_1/du_phys(0)| = ...     initial primary concentrations
  max |dJ_1/dK(0)| = ...          initial kinetic inventories (if any)
  max |dJ_1/dP(0)| = ...          initial pressures (with flow)
=== END ADJOINT ===
```

The initial concentrations and kinetic inventories are independent: dJ/du0
holds the inventories fixed, also where the input initializes them from
the concentrations. Both pull dJ/dY(0) back through the initial state;
verify_girt_dae.py Taylor-tests dJ/dY(0), kinetic and multirate entries
included (a deck's initial volume fraction also sets the reference of the
surface area, so finite differences in it are not dJ/dK(0)).
`OUTPUT FULL_FIELD` in the block also writes
`<prefix>-adjoint-dJ<i>-du0.bin` (PETSc binary Vec, entry
`(natural cell id-1)*ncomp + species-1`) and, with flow,
`<prefix>-adjoint-dJ<i>-dP0.bin` (entry `natural cell id-1`).

## Observations

An `OBSERVATIONS` sub-block of `ADJOINT` holds data at any time of the run
and makes one more cost, their misfit
phi = 1/2 sum_i ((H_i - d_i)/sigma_i)^2, with H_i the simulated value:

```
ADJOINT
  OBSERVATIONS
    TIME_SERIES
      CELL_ID 1                          # or COORDINATE x y z
      OBSERVED_VARIABLE FREE_ION CO2(aq) # or IMMOBILE <species>
      STANDARD_DEVIATION 1.d-6           # or ERROR ABSOLUTE a RELATIVE r
      FILE co2.dat                       # TIME_UNITS, DATA_UNITS, then
    /                                    #   rows: time value
    MEASUREMENT
      TIME 3 h
      CELL_ID 1
      OBSERVED_VARIABLE FREE_ION O2(aq)
      VALUE 1.2d-3 mol/kg
      STANDARD_DEVIATION 1.d-5
    /
  /
  ...
END
```

Nothing has a default: every `TIME_SERIES` and `MEASUREMENT` gives its
cell, variable, error and times with units; sigma = a + r|d| is in the
data's units (mol/kg aqueous, mol/m^3 bulk immobile). An observation is
at one cell (a `REGION` is rejected: it would need an operator such as a
volume average); a `COORDINATE` on a face between cells of two processes
(unstructured grids) is rejected too. `TOTAL` concentrations are not
supported yet; for a species without aqueous complexes the free ion is the
total. Times are after the start and at most the final time.

Observation times are waypoints (they shorten the run's steps), so the
run has a step state at each;
the forward run records the simulated values there, and the backward sweep
runs in segments between those steps, adding (H_i - d_i)/sigma_i^2 dH_i/dY
to the objective's lambda at each: one sweep for all times. The summary
prints the objective as the last cost J_k, with chi^2/n = 2 J_k/n;
`<prefix>-observations.dat` lists every observation (time, cell, observed,
simulated, sigma, residual) and `<prefix>-adjoint-gradient.dat` the
objective and its gradient (parameter name, value, derivative).
`-adjoint_debug_split` runs the sweep in segments at every waypoint
interval (a test: the result is the same bit for bit).

## Optimization

An `OPTIMIZATION` sub-block of `ADJOINT` minimizes the misfit of the
`OBSERVATIONS` (required) over the parameters with PETSc TAO, from the
parameters' deck values:

```
ADJOINT
  OBSERVATIONS ... /
  PARAMETER RATE_CONSTANT ... /
  OPTIMIZATION
    TYPE BQNLS                       # default; a first-order TAO type
    MAXIMUM_ITERATIONS 50            # > 0; TAO's default otherwise
    GRADIENT_TOLERANCE 1.d-8 1.d-8   # absolute, relative >= 0 (TAO's defaults)
    PARAMETER 1                      # deck order of the PARAMETERs
      TRANSFORM LOG10                # default but for SMOOTHSTEP_INTERVAL
                                     # (NONE), already a log10 width
      BOUNDS 1.d-3 1.d0              # in the parameter's units
    /
  /
END
```

All parameters must be positive: a point where one is not (or where
10**q underflows or overflows) is not run and returns +Inf. The forward run
is quiet and is TAO's first evaluation. Each further
evaluation reruns the run in place from its recorded initial state with the
recorded forcing replayed (no input is read again), and its adjoint sweep;
a last rerun at the solution writes the output, the adjoint summary and a
`-- GENERIC: Optimization --` section of the regression file (iterations,
converged reason, objective, parameters). An evaluation without a gradient
(a failed forward run, a mineral regime transition, a parameter Jacobian
that cannot be evaluated) returns +Inf, a step out of the domain; at the
initial parameters it stops the run (run the deck without OPTIMIZATION to
see why). The line search is Armijo's unless `-adj_tao_ls_type` sets one:
it backtracks from +Inf, while TAO's More-Thuente stops there; the
trust-region methods (`bqnktr`, `bqnktl`) do not handle it either. The
wall-clock limit (`WALLCLOCK_STOP`) inside an evaluation ends the
optimization at TAO's last accepted point, without a rerun at it (no
output is written); inside the rerun at the solution it
stops the run as well (the output is incomplete), and a failure there fails
the run. The run's statistics (steps, iterations, cuts) are those of the
forward run; its times include the optimization. Every evaluation is a row
of
`<prefix>-optimization.log` (number, status: 0 a value, -1 stopped, -2 not
run, outside the parameters' domain, else none; f, |df/dp|, the parameters
and df/dp); the screen shows f and |df/dp|. TAO's options take the prefix
`adj_` (`-adj_tao_type nm`, `-adj_tao_test_gradient`, `-adj_tao_monitor`).
The gradient is that of the steps taken, so with adaptive steps the
objective has small jumps where the step sequence changes; fixed steps
(`-tran_ts_adapt_type none`) avoid them.

## Scope

- Costs: `TERMINAL_SPECIES_CONCENTRATION` of a primary or an immobile
  species (a sensitivity at the final time), and the misfit of
  `OBSERVATIONS` at any times (free ion or immobile species at a cell).
- Parameters: `RATE_CONSTANT` (`KINETICS MICROBIAL` with `REACTION_INDEX`,
  or `KINETICS MINERAL` with `MINERAL_NAME` and a single rate constant),
  `HALF_SATURATION_CONSTANT`, `INHIBITION_CONSTANT`, `SMOOTHSTEP_INTERVAL`,
  `YIELD`, `DECAY_RATE_CONSTANT`, and, with flow in the model,
  `PERMEABILITY` (`MATERIAL <name>`: scales the isotropic permeability of
  the material's cells; no permeability datasets). Each parameter is the
  positive value PFLOTRAN holds (internal units): `INHIBITION_CONSTANT` is
  the
  `THRESHOLD_CONCENTRATION` as given, also with `INHIBIT_ABOVE_THRESHOLD`
  (stored as its negative), `SMOOTHSTEP_INTERVAL` the log10 interval,
  `PERMEABILITY` the material's permeability [m^2].
- Everything the GIRT_DAE model supports is differentiated, including
  activity coefficients (converged in every evaluation), kinetic minerals
  with volume-fraction-dependent surface area, multirate and equilibrium
  surface complexation, transient Richards flow and transient velocity
  datasets (the backward sweep replays each waypoint interval's recorded
  forcing).
- dJ/dP(0) holds the initial molarity of the aqueous primaries fixed, as the
  initial condition does. A constraint re-equilibrated in each cell against
  its flow state (activity, mineral or gas constraints) can depend on the
  initial pressure in ways dJ/dP(0) omits.
- Integrator: the block switches the forward run from BDF to backward
  Euler, the integrator PETSc's adjoint supports. Its step sizes may adapt:
  the gradient is that of the steps taken, without their dependence on the
  parameters (the finite-difference checks use fixed steps).
- Trajectory: the default (basic, on disk) or `-ts_trajectory_type memory`
  (keeps every step, including the one that ends each waypoint interval's
  `TSSolve`).

Rejected:

- integrators other than backward Euler (BDF has no PETSc adjoint), other
  trajectory types, `-ts_trajectory_solution_only`, and with the memory
  trajectory checkpoint limits (`-ts_trajectory_max_cps_*`,
  `-ts_trajectory_max_units_*`);
- `PERMEABILITY` without flow in the model (GIRT_DAE alone has no flow);
- `PERMEABILITY` with a source/sink scaled by permeability
  (`SCALED_MASS_RATE`/`SCALED_VOLUMETRIC_RATE` but by `VOLUME`): dF/dk omits
  the scaling's derivative;
- every configuration the models themselves reject.

The adjoint is skipped, with exit code 89 (input errors exit with 87),
after a forward run with mineral regime transitions (depletion, dissolution
resuming) or one that stops before the final time
(`MAXIMUM_NUMBER_OF_TIMESTEPS`, `WALLCLOCK_STOP`).

## Flow and transport as one model

With `MODE RICHARDS_DAE` flow and `MODE GIRT_DAE` transport, one PETSc TS
integrates both (`Flow_RT_DAE_Model_class`, state `[P; W]` and `[c; N; K]`
per cell); transport reads saturation, porosity, water density, Darcy
fluxes and source rates from the flow part. On the steps of the classic
sequential scheme (RICHARDS, then GIRT), backward Euler solves the same
equations.

- Solver input: `NUMERICAL_METHODS TRANSPORT` (read as for GIRT_DAE) and the
  `-tran_` options; `-flow_ksp_*`/`-flow_pc_*` configure the flow block of
  the preconditioner; `-ts_model_flow_atol` is the flow's absolute error
  tolerance.
- Newton holds the transport update while the flow update exceeds the
  flow's error tolerance, so transport linearizes about a settled flow
  state as in the sequential scheme; `-ts_model_flow_first 0` takes full
  coupled steps.
- Rejected with flow: transverse dispersion, radioactive decay, a specified
  velocity dataset, gas transport, and what RICHARDS_DAE and GIRT_DAE reject
  alone.

## Regression tests and checks

`adjoint.cfg`, part of `make test`: each registered deck runs once, and its
regression file (the forward solution, and the Adjoint section with the
costs, the gradients and the largest initial-state sensitivities) is
compared with its gold. `adjoint-optimization` recovers a microbial rate
constant (the Optimization section too); at the optimum the costs and
gradients are round-off and only the rate constant and TAO's counts are
compared.

`adjoint_verify.cfg`, run by `make adjoint_verify` (or `make verify`, with
the GIRT_DAE checks) because it takes minutes: each case runs its deck, then
a post-processor that checks the gradients against central differences
(mostly at two step sizes) on fixed backward Euler steps.

- `check_adjoint.py` (`adjoint`, `adjoint-np2`) and `check_extensions.py`
  (`adjoint-extensions`, `adjoint-extensions-np2`): the GIRT_DAE adjoint,
  selectors, output, restoration of the final state, the memory trajectory,
  a failed stage (`-rt_dae_fail_closure`: the gradient of the steps taken),
  reruns (`-ts_model_rerun_check`), optimization of a twin experiment
  (from 3x and 0.3x the rate constant, through failed evaluations and the
  wall-clock limit; TAO's gradient test, Nelder-Mead, bounds), parallel
  runs and the rejections.
- `check_coupled.py` (`adjoint-coupled`, `adjoint-coupled-np2`): the
  coupled model's gradients, and that it reproduces the sequential scheme
  and its Jacobian entry by entry, also on an unstructured grid; reruns,
  and optimization of permeability and calcite rate constant.
- `check_applications.py` (`adjoint-applications`,
  `adjoint-applications-np2`): a microbial batch (`abcd_batch.in`), a radial
  aquifer with flow in the model (`cylindrical_coupled.in`), and a 400-cell
  vertical section with gravity-driven flow (`adjoint-applications.in`). The
  radial permeability gradient uses Richardson-extrapolated central differences:
  the injection sets the flow, so the permeability acts only through the
  water density.

Run one directly (parallel runs need `PFLOTRAN_TEST_MPIEXEC` or `mpiexec`
on the path):

```
python3 check_<name>.py --executable ../../../src/pflotran/pflotran [--np 2]
```

Debug options:

- `-adjoint_debug_jacobianp`: checks dF/dp against central differences at
  the first evaluation of the backward sweep (`-adjoint_debug_jacobianp_eps`,
  default 1e-6, relative to the parameter).
- `-rt_dae_fail_closure <n>`: the n-th preparation of the GIRT_DAE model
  (its closures at a new state) fails every local cell, a domain error (a
  test of the recovery; Newton's step guard keeps the decks' iterates in
  the domain).
- `-adj_opt_fail_eval <n>`, `-adj_opt_wallclock_eval <n>`: the rerun of
  optimization evaluation n fails at its second step, or reaches the
  wall-clock limit after its first waypoint interval (tests of the failure
  paths).
- `-ts_model_rerun_check`: after the run's sweep, reruns at the deck's
  parameters p0, at 1.1 p0 and at p0 again; both reruns at p0 must
  reproduce the run's costs, gradients, adjoint initial state and final
  state bit for bit; prints the values at 1.1 p0.
- `-ts_model_contract_test`, `-ts_model_column_check`,
  `-ts_model_adjoint_test`: check the model contract, its Jacobian entry by
  entry, and the initial-state gradient. The column check cannot resolve
  terms below 1e-4 of an entry, such as those of water compressibility.
  At hydrostatic equilibrium every face has zero flux, where the upwinded
  fluxes are not differentiable: central differences then disagree with
  the (one-sided) Jacobian, so check such states off equilibrium.
