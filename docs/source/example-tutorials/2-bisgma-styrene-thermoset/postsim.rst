Post-build simulations and analyses
-----------------------------------

The BisGMA/styrene system is a cross-linked thermoset, so the full postsim and analyze suite applies: annealing and equilibration to densify the as-built network, a temperature ladder to estimate the glass-transition temperature, three uniaxial deformations to estimate Young's modulus, and a free-volume measurement on the equilibrated state.

The :ref:`pde_tutorial` walks through each of these stages in detail, including representative production-quality plots; this section follows the same pattern with parameters scaled to the BisGMA/styrene chemistry.

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
        Thi: 850
        deltaT: -5
        ps_per_rise: 10
        ps_per_run: 90
        warmup_ps: 10
    - ladder:
        input_top: systems/final-results/final.top
        input_gro: postsim/ladder-cool/ladder.gro
        subdir: postsim/ladder-heat
        Tlo: 100
        Thi: 850
        deltaT: 5
        ps_per_rise: 10
        ps_per_run: 90
        warmup_ps: 10
    - deform:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/deform-x
        T: 300
        P: 1
        direction: x
        edot: 0.001
        ps: 40
        replicas: 3
    - deform:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/deform-y
        T: 300
        P: 1
        direction: y
        edot: 0.001
        ps: 40
        replicas: 3
    - deform:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/deform-z
        T: 300
        P: 1
        direction: z
        edot: 0.001
        ps: 40
        replicas: 3

The two ladders run from 100 to 850 K, 151 rungs of 100 ps each, so
15.1 ns apiece.  That is wide enough for the whole-curve *T*:sub:`g` fit in
``plots post`` to see past both ends of the transition.  The heating
fit needs the top this high: on a 700 K ladder it was refused, and on 850 K ladders
five builds give 464 +/- 16 K, four of them on their own, with the transition ending
near 750 K (htpolynet-sweep, 2026-11-06).  Above about 700 K the ladder is beyond where the force field was checked for this chemistry,
so look at the top rungs' densities for anything anomalous before trusting the
rubbery line.

The cooling ladder runs first, starting from the equilibrated structure, and the
heating ladder starts from where it ends.  Run the other way round, the heating
ladder's 10 ps warmup would drop the 300 K structure to 100 K, a quench that distorts
the glassy branch the fit depends on.  The cooling ladder instead begins with a jump
from 300 K to 850 K in its warmup.  That was tested on example 1 (htpolynet-sweep, 2026-11-06): against cooling from
an already-hot structure, on the same builds, neither *T*:sub:`g`
moved significantly over five builds (heating -4.4 +/- 4.8 K, cooling +6.9 +/- 5.8 K),
and the top cooling rungs show no sign of an unrelaxed start.

Run it:

.. code-block:: console

    $ htpolynet postsim -cfg postsim.yaml -ocfg 2-bisgma-styrene-thermoset.yaml -proj proj-0

Density during annealing and equilibration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. admonition:: Placeholder
   :class: caution

   **TODO:** insert the ``rho_v_ns.png`` plots from ``proj-0/postsim/anneal/`` and ``proj-0/postsim/equilibrate/`` here, and comment on whether the density has stabilized.

Glass-transition temperature and Young's modulus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After ``postsim`` finishes, fit T\ :sub:`g` and *E* from the ladder and deformation traces with:

.. code-block:: console

    $ htpolynet plots post --cfg postsim.yaml --proj proj-0

This emits ``tg.png``, ``e.png``, ``E.csv``, and per-ladder CSVs in ``proj-0/plots/``.

.. admonition:: Placeholder
   :class: caution

   **TODO:** insert ``proj-0/plots/tg.png`` and ``proj-0/plots/e.png`` here, report the heating/cooling T\ :sub:`g` values and *E* (in GPa) from the console output, and discuss expected production-length values for a BisGMA/styrene thermoset.

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

   **TODO:** report the fractional free volume (and standard error) from ``proj-0/analyze/freevolume/ffv.dat`` here.  Cross-linked thermosets typically show FFV near 0.18–0.22 in the glassy state.
