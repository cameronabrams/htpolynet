.. _dfe_configuration:

The Configuration File
----------------------

The complete ``4-dfda-fde-epoxy-thermoset.yaml`` from the depot:

.. literalinclude:: ../../../../src/htpolynet/resources/example_depot/4-dfda-fde-epoxy-thermoset.yaml
   :language: yaml

Most blocks follow the same pattern as
:ref:`example 3 <pde_configuration>` — same ``densification``,
``precure``, ``CURE``, and ``postcure`` cascades — so the discussion
here focuses on what is specific to this build.

``constituents``
^^^^^^^^^^^^^^^^

A 2:1 stoichiometric ratio (200 FDE : 100 DFA) gives at most 400 C–N
bonds at full conversion — the same maximum-bond count as example 3.

.. code-block:: yaml

   constituents:
     FDE:
       smiles: "N(Cc1occc1)(C[CH:3]([OH:5])[CH3:1])C[CH:4]([OH:6])[CH3:2]"
       reactive_atoms: {1: C1, 2: C2, 3: C3, 4: C4, 5: O1, 6: O2}
       count: 1600
       symmetry_equivalent_atoms: [[C1, C2], [C3, C4], [O1, O2]]
       stereocenters: [C3]
     DFA:
       smiles: "C(c1oc([CH2][NH2:1])cc1)c1oc([CH2][NH2:2])cc1"
       reactive_atoms: {1: N1, 2: N2}
       count: 800
       symmetry_equivalent_atoms: [[N1, N2]]

The mapping onto example 3 is direct: FDE.{C1,C2,C3,C4,O1,O2}
corresponds to DGE.{C1,C2,C3,C4,O1,O2}; DFA.{N1,N2} corresponds to
PAC.{N1,N2}.  The downstream behaviour — reaction expansion, the
16-diastereomer pool for FDE, the cap reaction set — also carries over.

Two small differences from example 3 worth noting:

* DFA carries no ``stereocenters`` declaration.  The aromatic furan
  rings are planar; the central methylene is not chiral.  So unlike
  PACM (which contributes 4 diastereomers via its two cyclohexyl
  stereocenters), DFA contributes only a single canonical form.
* The furan rings are aromatic five-membered heterocycles.  GAFF/AMBER
  handle them, but during ``setup`` the AmberTools chain spends a bit
  more time on type-perception than for the all-carbocycle systems —
  watch for that in the ``profile.json``.

What the expansion produces
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the same reasons spelled out in
:ref:`example 3 <pde_configuration>`, the 2 cure + 1 cap reactions in
the YAML expand into 4 + 8 = 12 cure templates plus 2 cap templates.
Together with the 16 FDE diastereomers and the canonical DFA, the
diagnostic log at the start of a run will report something close to:

.. code-block:: text

   INFO> Templates in proj-0/molecules/parameterized
   INFO> N molecules detected in 4-dfda-fde-epoxy-thermoset.yaml
   INFO>                       explicit: ...
   INFO>     implied by stereochemistry: ...
   INFO>            implied by symmetry: ...

(Exact counts will appear when you run the build.)

Now we can turn to actually :ref:`running the build <dfe_run>`.
