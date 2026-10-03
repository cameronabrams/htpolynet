.. _badcy_postsim:

Post-build simulations and analyses
-----------------------------------

.. admonition:: Placeholder
   :class: caution

   **TODO:** everything below this box comes from a build of the pre-2.15.0 stand-in
   model (bisphenol A + pre-formed triazine + postcure repair) and will be replaced
   with results from the reference build of the cyclotrimerization configuration.

The canonical worked example for the postsim + analyze subsystems is
:ref:`tutorial 3 <tutorials_postsim_analyses>`; the workflow for
BADCy is identical save for the input filenames and a few
system-specific timings.  This page lists the postsim YAML and a few
BADCy-specific notes.

postsim.yaml
^^^^^^^^^^^^

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
        ps: 10
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
        ps: 10
    - deform:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/deform-y
        T: 300
        P: 1
        direction: y
        edot: 0.001
        ps: 10
    - deform:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/deform-z
        T: 300
        P: 1
        direction: z
        edot: 0.001
        ps: 10

Run it:

.. code-block:: console

    $ htpolynet postsim -cfg postsim.yaml -ocfg 6-cyanate-ester.yaml -proj proj-0

Density during annealing and equilibration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   The density this stage settles at is **not** an equilibrated density.
   This example's postcure anneal peaks at 500 K, and the measured
   *T*:sub:`g` for this system is about 488 K -- so the anneal spends its
   80 ps of peak time only ~12 K above the glass transition, where a
   crosslinked network barely moves.  On four independent builds the
   postcure plateau came out 2.31% below the same systems' density after a
   melt and slow re-cool, and 2.6-2.8% below experiment, where the
   re-cooled value agrees with experiment to under 1%.

   The 500 K peak is fine for this tutorial's purpose, which is to show the
   workflow.  It is not a basis for quoting a density.  Take densities from
   the cooling ladder below instead.


.. admonition:: Placeholder
   :class: caution

   **TODO:** insert the ``rho_v_ns.png`` plots from
   ``proj-0/postsim/anneal/`` and ``proj-0/postsim/equilibrate/``
   here, and note where the density of the repaired system settles.
   Expectation: ~1.0-1.1 g/cm³, somewhat below fully cyclotrimerized
   BADCy (~1.2 g/cm³) because the repaired ``-C#N`` end groups break
   the network into smaller clusters.

Glass-transition temperature and Young's modulus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After ``postsim`` finishes, fit Tg and E with:

.. code-block:: console

    $ htpolynet plots post --cfg postsim.yaml --proj proj-0

This emits ``tg.png``, ``e.png``, ``E.csv``, and per-ladder CSVs in
``proj-0/plots/``.

.. admonition:: Placeholder
   :class: caution

   **TODO:** insert ``proj-0/plots/tg.png`` and ``proj-0/plots/e.png``
   here, report the heating and cooling Tg estimates and the three
   per-direction *E* values, and discuss them in light of the
   network's actual crosslink density (which the repair stage
   reduces relative to a "complete cure" model).

Free volume
^^^^^^^^^^^

Use ``htpolynet analyze`` to invoke ``gmx freevolume`` on the
equilibration trajectory.  Create ``fv.yaml``:

.. code-block:: yaml

    - command: freevolume

Then:

.. code-block:: console

    $ htpolynet analyze -cfg fv.yaml -proj proj-0

.. admonition:: Placeholder
   :class: caution

   **TODO:** report the fractional free volume from
   ``proj-0/analyze/freevolume/ffv.dat``.  Expectation: somewhat
   higher than fully cured BADCy because the ``-C#N`` end groups are
   non-percolating chain stubs that introduce extra void.

A note on interpretation
^^^^^^^^^^^^^^^^^^^^^^^^

The thermomechanical observables on this system reflect a network at
90 % topological conversion *plus* postcure repair of the residual
artefacts.  Comparing to a real BADCy thermoset experiment, you'd
generally expect:

* a Tg shifted *down* relative to fully cyclotrimerized BADCy because
  the repaired ``-C#N`` end groups are not crosslinks;
* a Young's modulus correspondingly lower for the same reason;
* a free volume slightly higher.

These are the same kinds of shifts you'd expect from an experimental
undercured sample, which is encouraging — the topological model plus
repair stage reproduces the *direction* of incomplete-cure effects on
bulk properties, even if the absolute numbers depend on the cure
depth you targeted.  The :ref:`postcure-repair user-guide page
<postcure_repair>` discusses this further.
