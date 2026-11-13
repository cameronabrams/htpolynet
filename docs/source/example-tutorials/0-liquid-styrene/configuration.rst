.. _liquid_styrene_configuration:

Configuration
-------------

The full ``0-liquid-styrene.yaml`` is short — no ``reactions``, no
``CURE``, no ``postcure``.  It is just enough to parameterize a single
monomer, build a low-density box, and run the densification + anneal
cascade.

.. literalinclude:: ../../../../src/htpolynet/resources/example_depot/0-liquid-styrene.yaml
   :language: yaml

A block-by-block walk-through:

``Title``
   A free-text label for the run; appears in logs and in the title
   field of generated topology/structure files.

``gromacs``
   Tells ``htpolynet`` how to invoke GROMACS.  ``gmx_options`` are
   passed to every ``gmx`` invocation; ``mdrun_options`` are passed to
   ``mdrun`` only.  The ``gpu_id: 0`` line is left in by default and
   the runtime will automatically strip it (with a warning) if the
   ``gmx`` build is not GPU-capable or no GPU is visible.

``ambertools``
   Selects the partial-charge method for ``antechamber``; ``gas`` is
   Gasteiger, which is fast.  It is not the charge model GAFF was
   parameterized with, so every example's numbers are tutorial numbers; see
   the note on :ref:`the tutorials page <example_tutorials>`.

``constituents``
   Declares the molecular species in the box.  Here, one species
   (``STY``, 1000 copies) generated from SMILES.  See
   :ref:`monomer <liquid_styrene_monomer>` for the atom-mapping
   details.

``densification``
   The post-pack equilibration cascade.  Starts at
   ``initial_density: 300 kg/m³`` (deliberately loose so the pack
   succeeds without overlaps) and runs minimization → 10 ps NVT at
   300 K → 200 ps NPT at 300 K, 10 bar.  The deliberately high pressure
   accelerates compaction.

``precure``
   Despite the name, ``precure`` is not chemically a curing step — it
   is the post-densification equilibration that ``htpolynet`` runs
   whether or not a cure is configured.  Here it is:

   * ``preequilibration`` — 200 ps NPT at 300 K, 1 bar (relax the
     pressure jolt from densification);
   * ``anneal`` — two 300 K → 600 K → 300 K cycles, 80 ps per cycle,
     useful for shaking the liquid out of any pack-induced local
     traps;
   * ``postequilibration`` — 100 ps NPT at 300 K, 1 bar.

The comment block at the top of the YAML also sketches an
``htpolynet analyze`` configuration (``ck.yaml``) you can run after the
build to compute free volume and similar properties; see
:ref:`results <liquid_styrene_results>` and the :doc:`postsim` section.
