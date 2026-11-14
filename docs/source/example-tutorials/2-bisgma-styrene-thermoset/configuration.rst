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
       smiles: "c1ccccc1[CH:1]=[CH2:2]"     # styrene
       reactive_atoms: {1: C1, 2: C2}
       activate: [[1, 2]]                 # the C=C the cure consumes
       count: 1200
     BPA:
       smiles: "CC(c1ccc([OH:1])cc1)(c1ccc([OH:2])cc1)C"
       reactive_atoms: {1: O1, 2: O2}
     HIE:
       # Glycidyl methacrylate.  `activate` opens the methacrylate C=C (2-1) and
       # the epoxide's O-C4 bond (5-4), where a BPA hydroxyl attacks.
       smiles: "[CH2:2]=[C:1]([CH3])C(=O)OC[CH:3]1[O:5][CH2:4]1"
       reactive_atoms: {1: C1, 2: C2, 3: C3, 4: C4}
       activate: [[2, 1], [5, 4]]
       stereocenters:
         - C1
         - C3
     GMA:
       count: 600

Four species, three of them generated from SMILES.  STY and HIE are
written as they really are, styrene and glycidyl methacrylate, and each
has an ``activate`` key naming the bonds the cure consumes; ``htpolynet``
derives their active forms from that (see :ref:`example 1 <ps_monomer>`).
HIE has two: its methacrylate C=C, and the epoxide bond a BPA hydroxyl
opens to build GMA.  Its epoxide O carries a label (``:5``) only so
``activate`` can name it; it is not in ``reactive_atoms``.  Four things
worth noting:

* **HIE's active form is the old HIE.**  With both bonds opened, HIE is a
  2-hydroxypropyl ester of isobutyric acid, the molecule earlier versions
  of this example wrote out by hand, atom for atom.

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

SMILES bracket atoms take an **explicit** hydrogen count: an atom written
without one, like ``[C:1]``, has **zero** implicit H.  In HIE's real form
that is right for ``C1``, which is the sp2 carbon of the C=C and carries
no H; ``activate`` then adds the one it has in the active form.  It is
wrong for ``C3``, the epoxide CH, which must be written ``[CH:3]``.  When
an active form is written by hand, the same rule bites the other way:
earlier versions of this example wrote ``[CH:1]``, because ``[C:1]`` would
have left the saturated ``C1`` at valence 3, ``antechamber`` would have
typed it ``c2`` (sp2) and the adjacent carbonyl carbon ``c2`` instead of
``c``, and ``tleap`` would have failed with "no angle parameter for
o-c2-os".

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
