.. _bgs_run:

Running the Build
-----------------

From inside the working directory containing
``2-bisgma-styrene-thermoset.yaml``:

.. code-block:: console

   $ htpolynet run -diag diagnostics.log 2-bisgma-styrene-thermoset.yaml &> console.log &

The stage layout under ``proj-N/systems/`` is the same as
:ref:`example 1 <ps_run>` — ``init/``, ``densification/``, ``precure/``,
``iter-K/`` (one per CURE iteration), ``capping/``, ``postcure/``,
``final-results/`` — plus ``plots/`` and ``profile.json`` at the project
root.

Setup
^^^^^

The ``setup`` stage is where the ``param``-stage reactions actually run.
After parameterizing the three primitive species, ``htpolynet`` assembles
``GM1`` and then ``GMA``, then enumerates the 16 GMA diastereomers
implied by HIE's ``stereocenters`` declaration, and finally builds the 32
chain-expanded cure templates.  An abbreviated view of the diagnostic log
during this stage:

.. code-block:: text

   INFO> Templates in proj-0/molecules/parameterized
   INFO> 32 molecules detected in 2-bisgma-styrene-thermoset.yaml
   INFO>                       explicit: 11
   INFO>     implied by stereochemistry: 21
   INFO>            implied by symmetry: 0
   INFO> AmberTools> generating GAFF parameters from BPA.mol2
   INFO> BPA: 228.28 g/mol
   INFO> AmberTools> generating GAFF parameters from HIE.mol2
   INFO> HIE: 146.18 g/mol
   INFO> AmberTools> generating GAFF parameters from STY.mol2
   INFO> STY: 106.16 g/mol
   INFO> AmberTools> generating GAFF parameters from GM1.mol2
   INFO> GM1: 372.44 g/mol
   ...
   INFO> AmberTools> generating GAFF parameters from GMA.mol2
   INFO> GMA: 516.61 g/mol
   ...
   INFO> Built GMAS-1 using topology of GMA; copying GMA.top to GMAS-1.top
   INFO> GMAS-1: 516.61 g/mol
   ...
   INFO> 32 molecules implied by chaining
   ...
   INFO> Generated 64 molecule templates
   INFO> Initial composition is STY 150, GMA 75
   INFO> 100% conversion is 300 bonds

GMA itself is fully GAFF-parameterized; the diastereomers ``GMAS-1`` …
``GMAS-15`` are *built* (atom positions reflected) but reuse GMA's
topology directly — they're geometry-only variants.

Densification
^^^^^^^^^^^^^

Eight consecutive 100 ps NPT segments compress the initial 100 kg/m³
liquid up to roughly the polymer's bulk density.  Sample trace from the
diagnostic log:

.. code-block:: text

   INFO> Initial density: 100.0 kg/m^3
   INFO> Initial box side lengths: 9.683 nm x 9.683 nm x 9.683 nm
   ...
   INFO> Repeat 6 out of 8
   INFO> Current box side lengths: 4.462 nm x 4.462 nm x 4.462 nm
   INFO> Density                      1021.60
   ...
   INFO> Repeat 8 out of 8
   INFO> Current box side lengths: 4.465 nm x 4.465 nm x 4.465 nm
   INFO> Density                      1020.05

By the end of densification the liquid is at roughly 1.02 g/cm³.

CURE iterations
^^^^^^^^^^^^^^^

With ``desired_conversion: 0.95`` and a max of 300 bonds, the target is
285 bonds.  Early iterations grab dozens of bonds at a time at the
default 0.5 nm search radius; as conversion approaches the target, the
remaining reactive pairs get sparser and ``htpolynet`` has to grow the
radius (and sometimes drag pairs in) to find anything.  A typical run
needs ~10–15 iterations at the default ``min_bonds_per_iteration: 10``.
Excerpted from a representative log:

.. code-block:: text

   INFO> ~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~
   INFO> Bond search using radius 0.5 nm initiated
   INFO> Iteration 1 will generate 82 new bonds
   ...
   INFO> Iteration 1 current conversion 0.273 or 82 bonds
   INFO> ~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~
   INFO> Bond search using radius 0.5 nm initiated
   INFO> Iteration 2 will generate 67 new bonds
   ...
   INFO> Iteration 2 current conversion 0.497 or 149 bonds
   ...
   INFO> ~~~~~~~~~~~~~~ Iteration 12 begins ~~~~~~~~~~~~~~
   INFO> Bond search using radius 0.5 nm initiated
   INFO> Radius increased to 0.75 nm
   INFO> Iteration 12 will generate 3 new bonds
   ...
   INFO> Iteration 12 current conversion 0.947 or 284 bonds
   INFO> ~~~~~~~~~~~~~~ Iteration 13 begins ~~~~~~~~~~~~~~
   INFO> Bond search using radius 0.5 nm initiated
   INFO> Radius increased to 0.75 nm
   INFO> Iteration 13 will generate 1 new bond
   INFO> Step "cure_drag" initiated on 1 distance (max 0.691 nm)
   ...
   INFO> Iteration 13 current conversion 0.950 or 285 bonds

.. tip::

   This is exactly the workflow ``CURE.controls.min_bonds_per_iteration``
   was designed to compress.  Setting it to e.g. ``10`` would make
   ``htpolynet`` grow the radius further on the late iterations so each
   one picks up at least 10 bonds before relaxation.  Whether that's a
   win depends on whether your relax cascade can absorb the larger
   strained-bond batches — try it on this example and watch how many
   iterations are needed.

Capping and postcure
^^^^^^^^^^^^^^^^^^^^

After the last cure iteration, ``htpolynet`` looks for any unreacted
``C1–C2`` sites:

.. code-block:: text

   INFO> Capping begins
   INFO> Capping will generate 0 new bonds
   INFO> Connect-Update-Relax-Equilibrate (CURE) ends

At 95% conversion in this run, all radicals had paired up — but on a
build that hits ``max_iterations`` before reaching ``desired_conversion``,
expect the cap reactions (``STY-restore-C1-C2`` and ``HIE-restore-C2-C1``;
the recorded run above predates them and called them ``styCC`` and
``hieCC``) to do non-trivial work.

Postcure runs two more 300/600 K annealing cycles plus a 100 ps NPT at
300 K, then ``htpolynet`` writes ``final.gro``, ``final.top``,
``final.tpx``, ``final.grx`` plus the VMD viz pair into
``proj-N/systems/final-results/`` and emits the run profile.

Now to look at some :ref:`results <bgs_results>`.
