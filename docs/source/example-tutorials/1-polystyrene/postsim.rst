Post-build simulations and analyses
-----------------------------------

Polystyrene is a linear polymer rather than a cross-linked network, so a uniaxial-deformation measurement of Young's modulus is not particularly meaningful for the short chains built here.  What we *can* do is anneal and equilibrate the system, then perform a temperature ladder to estimate the glass-transition temperature, and use ``analyze`` to measure the fractional free volume of the equilibrated melt.

The :ref:`pde_tutorial` walks through the full set of postsim and analyze techniques in detail (including deformation, which is well-suited to thermosets); this section adapts that workflow to a linear polymer.

The postsim configuration
^^^^^^^^^^^^^^^^^^^^^^^^^

Create ``postsim.yaml`` in the project base directory:

.. code-block:: yaml

    - anneal:
        input_top: systems/final-results/final.top
        input_gro: systems/final-results/final.gro
        P: 1
        T0: 300
        T1: 600
        ncycles: 2
        T0_to_T1_ps: 10
        T1_ps: 10
        T1_to_T0_ps: 10
        T0_ps: 10
    - equilibrate:
        input_top: systems/final-results/final.top
        input_gro: postsim/anneal/anneal.gro
        T: 300
        ps: 100
    - ladder:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/ladder-cool
        Tlo: 100
        Thi: 700
        deltaT: -5
        ps_per_rise: 10
        ps_per_run: 90
        warmup_ps: 10
    - ladder:
        input_top: systems/final-results/final.top
        input_gro: postsim/ladder-cool/ladder.gro
        subdir: postsim/ladder-heat
        Tlo: 100
        Thi: 700
        deltaT: 5
        ps_per_rise: 10
        ps_per_run: 90
        warmup_ps: 10

The two ladders run from 100 to 700 K, 121 rungs of 100 ps each, so
12.1 ns apiece.  That is wide enough for the whole-curve *T*:sub:`g` fit in
``plots post`` to see past both ends of the transition.

The cooling ladder runs first, starting from the equilibrated structure, and the
heating ladder starts from where it ends.  Run the other way round, the heating
ladder's 10 ps warmup would drop the 300 K structure to 100 K, a quench that distorts
the glassy branch the fit depends on.  The cooling ladder instead begins with a jump
from 300 K to 700 K in its warmup.  That was tested on example 1 (htpolynet-sweep, 2026-11-06): against cooling from
an already-hot structure, on the same builds, neither *T*:sub:`g`
moved significantly over five builds (heating -4.4 +/- 4.8 K, cooling +6.9 +/- 5.8 K),
and the top cooling rungs show no sign of an unrelaxed start.

Run it:

.. code-block:: console

    $ htpolynet postsim -cfg postsim.yaml -ocfg 1-polystyrene.yaml -proj proj-0

This produces ``proj-0/postsim/{anneal,equilibrate,ladder-heat,ladder-cool}/`` with the usual Gromacs outputs plus on-the-fly density-vs-time plots.

Density during annealing and equilibration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. admonition:: Placeholder
   :class: caution

   **TODO:** insert the ``rho_v_ns.png`` plots from ``proj-0/postsim/anneal/`` and ``proj-0/postsim/equilibrate/`` here.  Note whether the annealing has actually densified the system (it usually fails to in builds this short) and whether the post-anneal equilibration trace is flat over the second half.

Glass-transition temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After ``postsim`` has produced the heating and cooling ladders, fit T\ :sub:`g` from the density-vs-temperature traces with:

.. code-block:: console

    $ htpolynet plots post --cfg postsim.yaml --proj proj-0

This emits ``tg.png`` plus per-ladder CSVs in ``proj-0/plots/``.

.. admonition:: Placeholder
   :class: caution

   **TODO:** insert ``proj-0/plots/tg.png`` and report the heating and cooling T\ :sub:`g` values from the console output here.  In simulations this short, expect significant heating/cooling hysteresis — comment on what production-length runs would be needed to get a converged estimate.

Free volume
^^^^^^^^^^^

Use ``htpolynet analyze`` to invoke ``gmx freevolume`` on the equilibration trajectory.  Create ``fv.yaml``:

.. code-block:: yaml

    - command: freevolume

Then:

.. code-block:: console

    $ htpolynet analyze -cfg fv.yaml -proj proj-0

.. admonition:: Placeholder
   :class: caution

   **TODO:** report the fractional free volume (and standard error) from ``proj-0/analyze/freevolume/ffv.dat`` here.  Glassy polystyrene at 300 K typically has FFV in the 0.15–0.20 range.
