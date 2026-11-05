Finite-resistance liquid Robin boundary (GENERAL)
================================================

Purpose and input
-----------------

``LIQUID_ROBIN_COEFFICIENT`` adds a constant hydraulic resistance between a
liquid reservoir and a domain boundary. It is a ``FLOW_CONDITION`` keyword,
outside the ``TYPE`` block. Use that condition in a ``BOUNDARY_CONDITION``::

  FLOW_CONDITION reservoir
    TYPE
      LIQUID_PRESSURE DIRICHLET
      MOLE_FRACTION DIRICHLET
      TEMPERATURE DIRICHLET
    /
    LIQUID_PRESSURE 200100 Pa
    MOLE_FRACTION 0
    TEMPERATURE 12 C
    LIQUID_ROBIN_COEFFICIENT 1.d-10
  END

  BOUNDARY_CONDITION inlet
    FLOW_CONDITION reservoir
    REGION inlet
  END

The coefficient is supplied in **m/(Pa s)**, without a unit suffix. The pressure
is the external reservoir pressure, not an imposed domain-face pressure.
``MOLE_FRACTION`` is the dissolved air mole fraction in the liquid reservoir.
The reservoir must be in the liquid state; domain cells may be liquid-only or
two-phase. The complete one-cell example is ``liquid_robin.in`` in this directory.

``LIQUID_PRESSURE HYDROSTATIC`` is also supported, with the usual ``DATUM`` and
pressure definition, provided the boundary reservoir remains in the liquid
state. Existing restrictions on GENERAL hydrostatic phase switching still apply.

Discrete boundary law
---------------------

Take positive liquid Darcy velocity into the domain. With reservoir pressure
``p_ext`` and face pressure ``p_face``, the external resistance obeys

.. math::

   q_l = \beta_l (p_{\mathrm{ext}} - p_{\mathrm{face}}).

The finite-volume cell-to-face Darcy coefficient and the effective coefficient
are

.. math::

   K_c = \frac{k_n}{d}\frac{k_{rl,\mathrm{up}}}{\mu_{l,\mathrm{up}}},
   \qquad
   K_{\mathrm{eff}} = \frac{K_c\beta_l}{K_c+\beta_l},
   \qquad
   q_l = K_{\mathrm{eff}}\Delta\Pi_l.

Here ``k_n`` is the permeability projected along the boundary connection,
``d`` is the cell-center-to-face distance, and ``Delta Pi_l`` is the existing
GENERAL pressure difference including its gravity correction. The implementation
uses the existing upwind mobility, composition, density and enthalpy. Boundary
area and component conversion are applied by the existing flux assembly.

The analytical Jacobian includes the dependence on liquid mobility ``M``:

.. math::

   \frac{\partial q_l}{\partial\Delta\Pi_l}=K_{\mathrm{eff}},
   \qquad
   \frac{\partial q_l}{\partial M}
     =\frac{k_n}{d}\left(\frac{\beta_l}{K_c+\beta_l}\right)^2\Delta\Pi_l.

The ratio is evaluated without directly multiplying ``K_c`` by ``beta_l`` or
forming their reciprocals. This avoids overflow for widely separated scales.
At zero mobility and positive beta, the one-sided mobility derivative is retained.

Limits and supported scope
--------------------------

* Omitted keyword: existing GENERAL boundary behavior is unchanged.
* Zero coefficient: liquid **Darcy** flow is zero.
* Large coefficient: the existing prescribed-pressure boundary is recovered.
* Finite positive coefficient: cell and boundary resistances act in series.
  Both inflow and outflow are allowed; this is not a one-way seepage valve.
* The coefficient already includes the external hydraulic resistance in
  pressure-based SI units. It is not interchangeable with ``CONDUCTANCE [m]``.
* Negative, NaN, infinite and duplicate coefficients are rejected.
* GENERAL Darcy flow with a liquid-state exterior is supported. Non-Darcy flow,
  a salt-bearing condition, flux/rate conditions, and pressure boundary types
  other than DIRICHLET/HYDROSTATIC are rejected.
* Use only on boundary conditions. Initial conditions and source/sink couplers
  cannot carry this keyword.

This keyword changes the liquid advective flux only. It does not impose zero
gas flow, suppress molecular diffusion, implement vapor exchange, or introduce
a bubble-entry pressure. Gas and diffusive fluxes retain their existing GENERAL
behavior. In particular, zero liquid Robin coefficient is not a fully sealed
two-phase boundary. Heat advection follows the liquid flux; the existing thermal
boundary treatment remains active. The coefficient itself is constant in time
and space for each named flow condition.

Verification
------------

The six serial ``liquid_robin*.in`` decks are registered in ``general.cfg`` under the
``standard`` and ``flow`` suites. They cover finite resistance, zero transfer,
reversed pressure, the large-coefficient limit, and two-phase cells with both
analytical and numerical Jacobians. A seventh, four-cell test
(``liquid_robin-np2``) is registered in ``standard_parallel`` and uses two MPI
ranks. Each has a native ``.regression.gold`` file.

Run them using the normal regression manager, for example from the repository
root::

  python3 regression_tests/regression_tests.py \
    -e /absolute/path/to/pflotran \
    -c regression_tests/general/general.cfg \
    -t liquid_robin liquid_robin_zero liquid_robin_reverse liquid_robin_limit \
       liquid_robin_two_phase liquid_robin_two_phase_numerical

``verify_liquid_robin.py`` provides independent physical checks in addition to
the stored regression values. It requires Python 3, NumPy and h5py. Do not run
with Python's ``-O`` option because the verifier uses assertions. Its output
directory must not already exist::

  python3 regression_tests/general/verify_liquid_robin.py \
    --executable /absolute/path/to/pflotran \
    --reference-executable /absolute/path/to/unmodified/pflotran \
    --output /tmp/liquid-robin-verification

The optional reference executable must be built from the unmodified parent
with the same compiler, dependencies and options. The verifier compares the
no-keyword field values bit for bit. It also checks:

* The steady, homogeneous, one-dimensional resistance solution
  ``q_l = Delta p / (mu L/k + 1/beta_l)`` with one and four cells.
* Zero, small and large coefficients, equilibrium and reverse flow.
* HYDROSTATIC/DIRICHLET equivalence for a horizontal connection at the
  reservoir datum elevation.
* Water and air inventory changes against integrated boundary mass fluxes.
* Finite fields and saturation closure.
* Rewetting and drainage using analytical/numerical Jacobians. PETSc's finite
  difference Jacobian is compared on the active water/air block; the fixed
  temperature row/column are constrained by the isothermal formulation.
* Nonisothermal inflow from a warmer liquid reservoir, comparing analytical
  and numerical Jacobian solutions.
* Rejection of invalid coefficients, unsupported boundary types, non-Darcy
  flow and use of the Robin condition as an initial condition.

The verifier requests 15 ASCII digits so rounding in mass-balance files does
not dominate the conservation check. For each component, the reported balance
error is the maximum absolute inventory-minus-flux discrepancy divided by
``max(initial component mass, 1 kg)``. Thus the denominator remains defined for
zero initial air mass. The JSON result records tolerances and measured errors. Build
PFLOTRAN with GNU ``-O0`` when regenerating native regression gold files, in
accordance with the project contribution guidance. Test outputs and compiled
files are not part of this contribution.
