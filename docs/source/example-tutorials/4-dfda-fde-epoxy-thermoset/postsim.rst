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
        subdir: postsim/ladder-heat
        Tlo: 300
        Thi: 600
        deltaT: 5
        ps_per_rise: 10
        ps_per_run: 90
        warmup_ps: 10
    - ladder:
        input_top: systems/final-results/final.top
        input_gro: postsim/ladder-heat/ladder.gro
        subdir: postsim/ladder-cool
        Tlo: 300
        Thi: 600
        deltaT: -5
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
