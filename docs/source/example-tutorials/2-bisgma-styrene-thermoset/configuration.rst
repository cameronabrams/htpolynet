.. _bgs_configuration:

The Configuration File
----------------------

The complete ``2-bisgma-styrene-thermoset.yaml`` from the depot:

.. literalinclude:: ../../../../src/htpolynet/resources/example_depot/2-bisgma-styrene-thermoset.yaml
   :language: yaml

The ``Title``, ``gromacs``, ``ambertools``, ``precure``, ``CURE``, and
``postcure`` blocks follow the same pattern as
:ref:`example 1 <ps_configuration>`.  The two interesting differences are
``constituents`` and ``densification``.

``constituents``
^^^^^^^^^^^^^^^^

.. code-block:: yaml

   constituents:
     STY:
       smiles: "c1ccccc1[CH2:1][CH3:2]"
       reactive_atoms: {1: C1, 2: C2}
       count: 1200
     BPA:
       smiles: "CC(c1ccc([OH:1])cc1)(c1ccc([OH:2])cc1)C"
       reactive_atoms: {1: O1, 2: O2}
     HIE:
       smiles: "[CH3:2][CH:1]([CH3])C(=O)OC[CH:3](O)[CH3:4]"
       reactive_atoms: {1: C1, 2: C2, 3: C3, 4: C4}
       stereocenters:
         - C1
         - C3
     GMA:
       count: 600

Four species, three of them generated from SMILES.  Three things worth
noting:

* **BPA and HIE have no ``count``.**  That's because they're consumed by
  the ``param``-stage reactions to assemble GMA; they don't appear in the
  initial liquid as their own species.  Their ``mol2`` files are still
  generated and parameterized so the param reactions have something to
  consume.
* **GMA has no ``smiles``.**  ``htpolynet`` knows GMA from the
  ``param``-stage reactions in the ``reactions`` block; it builds the
  template by running ``B1`` and ``B2`` at startup.  GMA's ``count: 600``
  is what actually populates the initial liquid (paired with 1200 STY,
  giving a 1:2 BisGMA:styrene ratio).
* **``stereocenters`` is declared on HIE, not GMA.**  HIE has two chiral
  carbons (``C1`` and ``C3``).  Each GMA contains two HIE moieties, so
  the racemic pool ``htpolynet`` builds for filling the initial box is
  16 GMA diastereomers (2 chiral centres × 2 HIE per GMA = 4 bits → 16
  configurations).  The diagnostics log at the start of a run lists the
  full set (e.g. ``GMAS-1`` through ``GMAS-15`` plus the canonical
  ``GMA``).

A note on the HIE SMILES
^^^^^^^^^^^^^^^^^^^^^^^^

The HIE SMILES uses ``[CH:1]`` (with an explicit ``H``) rather than
``[C:1]``.  This matters: SMILES bracket atoms without an explicit
hydrogen count default to **zero** implicit H, which would leave atom
``C1`` at valence 3.  ``antechamber`` then mistypes it as ``c2`` (sp2),
mistypes the adjacent carbonyl carbon as ``c2`` instead of ``c``, and
``tleap`` later fails with "no angle parameter for o-c2-os".  Writing
``[CH:1]`` keeps the valence at 4 and the typing correct.  The same
applies to ``[CH:3]``.

``densification``
^^^^^^^^^^^^^^^^^

.. code-block:: yaml

   densification:
     initial_density: 100.0  # kg/m3
     equilibration:
       - ensemble: min
       - ensemble: nvt
         temperature: 300
         ps: 10
       - ensemble: npt
         temperature: 300
         pressure: 10
         ps: 100
         repeat: 8

GMA is a long, awkward-shaped molecule; ``gmx insert-molecules`` has
trouble packing it cleanly at higher densities, so the initial density
is set to a deliberately loose ``100 kg/m³``.  The NPT compression step
then uses ``repeat: 8`` to run eight 100 ps NPT segments back-to-back —
serial restarts let ``gmx`` rebuild its neighbour lists between
segments, which avoids the "box size changed too much" warning that
would otherwise trip a single long run.
