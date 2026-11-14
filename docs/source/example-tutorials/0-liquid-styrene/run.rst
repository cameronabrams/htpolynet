.. _liquid_styrene_run:

Running the Build
-----------------

From inside the working directory containing ``0-liquid-styrene.yaml``:

.. code-block:: console

   $ htpolynet run -diag diagnostics.log 0-liquid-styrene.yaml &> console.log &

This kicks off the full workflow:

1. **Monomer setup.**  ``htpolynet`` reads the SMILES from
   ``constituents.STY.smiles``, builds ``lib/molecules/inputs/STY.mol2``,
   then invokes ``antechamber``/``parmchk2``/``tleap`` to produce the
   parameterized monomer in ``lib/molecules/parameterized/``.
2. **Initial pack.**  1000 copies of styrene are placed into a box sized
   for ``initial_density: 300 kg/m³``.
3. **Densification.**  The cascade under ``densification`` runs:
   minimization, 10 ps NVT, then 200 ps NPT at 10 bar *repeated until the
   density settles*.  On a 24-core run of the earlier version of this
   example (1000 ethylbenzenes, 18,000 atoms, rather than today's 1000
   styrenes, 16,000 atoms) that took two extra segments, so 600 ps of NPT
   rather than 200, and densification as a whole took about 5.5 minutes.
4. **Anneal cascade.**  The cascade under ``precure`` runs:
   pre-equilibration, two 300/600 K cycles, post-equilibration.

Each stage writes into its own subdirectory under ``proj-N/systems/``,
where ``N`` is the next available index (the line ``Working in new
project proj-N`` near the top of ``diagnostics.log`` tells you which
one was chosen).  For this no-cure build you'll see:

.. code-block:: text

   proj-0/
   └── systems/
       ├── init/             # initial packed box
       ├── densification/    # min + nvt + npt cascade output
       ├── precure/          # preequilibration, anneal, postequilibration
       ├── postcure/         # empty -- no postcure block in this example
       └── final-results/    # the final structure + topology + viz files

There are no ``capping/`` or ``iter-*/`` subdirectories — those are only
produced when a ``CURE`` block is configured.

If you tail ``diagnostics.log`` while the run is going, the major stage
transitions are clearly logged.  A clean finish ends with an
``htpolynet runtime ends`` line and the ``final-results`` directory
populated.
