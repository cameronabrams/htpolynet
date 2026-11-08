.. _htpb_postsim:

Post-build simulations and analyses
-----------------------------------

The canonical worked example for the postsim + analyze subsystems is
:ref:`tutorial 3 <tutorials_postsim_analyses>`; the workflow for
HTPB/IPDI follows the same pattern.  This page lists the postsim
YAML and a few HTPB-specific notes about what to expect from the
thermomechanical observables.

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
        subdir: postsim/ladder-cool
        Tlo: 50
        Thi: 500
        deltaT: -5
        ps_per_rise: 10
        ps_per_run: 90
        warmup_ps: 10
    - ladder:
        input_top: systems/final-results/final.top
        input_gro: postsim/ladder-cool/ladder.gro
        subdir: postsim/ladder-heat
        Tlo: 50
        Thi: 500
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

The two ladders run from 50 to 500 K, 91 rungs of 100 ps each, so
9.1 ns apiece.  That is wide enough for the whole-curve *T*:sub:`g` fit in
``plots post`` to see past both ends of the transition.  This elastomer's
*T*:sub:`g` is near 200-220 K, so it is the bottom of the ladder that matters here:
an earlier 250 K floor sat inside the transition.

The cooling ladder runs first, starting from the equilibrated structure, and the
heating ladder starts from where it ends.  Run the other way round, the heating
ladder's 10 ps warmup would drop the 300 K structure to 50 K, a quench that distorts
the glassy branch the fit depends on.  The cooling ladder instead begins with a jump
from 300 K to 500 K in its warmup.  That was tested on example 1 (htpolynet-sweep, 2026-11-06): against cooling from
an already-hot structure, on the same builds, neither *T*:sub:`g`
moved significantly over five builds (heating -4.4 +/- 4.8 K, cooling +6.9 +/- 5.8 K),
and the top cooling rungs show no sign of an unrelaxed start.

Run it:

.. code-block:: console

    $ htpolynet postsim -cfg postsim.yaml -ocfg 5-htpb-ipdi.yaml -proj proj-0

Density during annealing and equilibration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. admonition:: Placeholder
   :class: caution

   **TODO:** insert the ``rho_v_ns.png`` plots from
   ``proj-0/postsim/anneal/`` and ``proj-0/postsim/equilibrate/``
   here.  Expected: density near 0.9 g/cm³ at 300 K (HTPB-based
   polyurethanes are less dense than the bisphenol-A thermosets in
   tutorials 2-3).

Glass-transition temperature and Young's modulus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After ``postsim`` finishes, fit Tg and E with:

.. code-block:: console

    $ htpolynet plots post --cfg postsim.yaml --proj proj-0

.. admonition:: Placeholder
   :class: caution

   **TODO:** insert ``proj-0/plots/tg.png`` and ``proj-0/plots/e.png``
   here.  Expected from HTPB elastomer chemistry: a Tg significantly
   below 300 K (HTPB-based polyurethanes typically have Tg in the
   range of 200-260 K depending on crosslink density), and a Young's
   modulus in the 1-10 MPa range — three orders of magnitude lower
   than the GPa-scale moduli of the rigid epoxy / cyanate-ester
   thermosets in the other tutorials.  This is what distinguishes
   "elastomer" from "thermoset glass" in the broader thermoset
   property space.

Free volume
^^^^^^^^^^^

.. code-block:: yaml

    - command: freevolume

.. code-block:: console

    $ htpolynet analyze -cfg fv.yaml -proj proj-0

.. admonition:: Placeholder
   :class: caution

   **TODO:** report the fractional free volume from
   ``proj-0/analyze/freevolume/ffv.dat``.  Expectation: higher than
   the rigid thermoset examples — flexible HTPB backbones pack less
   efficiently than aromatic ring backbones do.

A note on system size
^^^^^^^^^^^^^^^^^^^^^

At about 111k atoms this is the largest depot example.  Postsim
times scale accordingly — expect each Tg ladder (91 temperatures ×
100 ps each = 9.1 ns) to take several hours on a single GPU,
and the deformation runs to take a few minutes each.  A useful
exercise is to compare wall-clock to the build itself: the postsim
suite is a small fraction (1-2 hours) of the ~12 hours the build
takes.
