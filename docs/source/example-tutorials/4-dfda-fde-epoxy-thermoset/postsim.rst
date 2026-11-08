Post-build simulations and analyses
-----------------------------------

DFA/FDE is a cross-linked epoxy thermoset with the same overall reaction mechanism as :ref:`example 3 <pde_tutorial>`, so the same postsim and analyze suite applies: annealing and equilibration, a temperature ladder for T\ :sub:`g`, uniaxial deformations for *E*, and a free-volume measurement.

The :ref:`pde_tutorial` walks through each stage in detail with representative production plots; this section follows the same pattern.

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
        Thi: 800
        deltaT: -5
        ps_per_rise: 10
        ps_per_run: 90
        warmup_ps: 10
    - ladder:
        input_top: systems/final-results/final.top
        input_gro: postsim/ladder-cool/ladder.gro
        subdir: postsim/ladder-heat
        Tlo: 100
        Thi: 800
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
        edot: 0.00025
        ps: 160
        replicas: 3
    - deform:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/deform-y
        T: 300
        P: 1
        direction: y
        edot: 0.00025
        ps: 160
        replicas: 3
    - deform:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/deform-z
        T: 300
        P: 1
        direction: z
        edot: 0.00025
        ps: 160
        replicas: 3

The two ladders run from 100 to 800 K, 141 rungs of 100 ps each, so
14.1 ns apiece.  That is wide enough for the whole-curve *T*:sub:`g` fit in
``plots post`` to see past both ends of the transition.  The top is 800 K because the
*T*:sub:`g` depends on it: five builds give a heating *T*:sub:`g` of 421 K on 700 K
ladders and 402 K on 800 K ones, with cooling moving the same way (htpolynet-sweep,
2026-11-06).  It keeps moving above 800 K: ladders run to 1000 K on three of the
builds give 406 K with the top cut to 800 K and 394 K at 1000 K, still falling.  So
treat this ladder's *T*:sub:`g` as an upper bound; ``plots post`` warns when it
moves with the top.  Going taller is not free: the bonded energies stay linear in
temperature and identical on heating and cooling up to 1000 K, so the force field
holds, but 4 of 10 runs on 1000 K ladders died, against 1 of 21 on 800 K ones.

The cooling ladder runs first, starting from the equilibrated structure, and the
heating ladder starts from where it ends.  Run the other way round, the heating
ladder's 10 ps warmup would drop the 300 K structure to 100 K, a quench that distorts
the glassy branch the fit depends on.  The cooling ladder instead begins with a jump
from 300 K to 800 K in its warmup.  That was tested on example 1 (htpolynet-sweep, 2026-11-06): against cooling from
an already-hot structure, on the same builds, neither *T*:sub:`g`
moved significantly over five builds (heating -4.4 +/- 4.8 K, cooling +6.9 +/- 5.8 K),
and the top cooling rungs show no sign of an unrelaxed start.

Run it:

.. code-block:: console

    $ htpolynet postsim -cfg postsim.yaml -ocfg 4-dfda-fde-epoxy-thermoset.yaml -proj proj-0

Density during annealing and equilibration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. admonition:: Placeholder
   :class: caution

   **TODO:** insert the ``rho_v_ns.png`` plots from ``proj-0/postsim/anneal/`` and ``proj-0/postsim/equilibrate/`` here, and comment on densification and equilibration adequacy.

Glass-transition temperature and Young's modulus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After ``postsim`` finishes, fit T\ :sub:`g` and *E* with:

.. code-block:: console

    $ htpolynet plots post --cfg postsim.yaml --proj proj-0

This emits ``tg.png``, ``e.png``, ``E.csv``, and per-ladder CSVs in ``proj-0/plots/``.

.. admonition:: Placeholder
   :class: caution

   **TODO:** insert ``proj-0/plots/tg.png`` and ``proj-0/plots/e.png`` here, report the heating/cooling T\ :sub:`g` values and *E* (in GPa), and contrast with the DGEBA/PACM numbers in :ref:`pde_tutorial` to highlight how chemistry choice shifts thermomechanical properties.

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

   **TODO:** report the fractional free volume (and standard error) from ``proj-0/analyze/freevolume/ffv.dat``.
