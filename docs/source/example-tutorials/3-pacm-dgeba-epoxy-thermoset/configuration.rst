.. _pde_configuration:

The Configuration File
----------------------

The complete ``3-pacm-dgeba-epoxy-thermoset.yaml`` from the depot:

.. literalinclude:: ../../../../src/htpolynet/resources/example_depot/3-pacm-dgeba-epoxy-thermoset.yaml
   :language: yaml

The ``Title``, ``gromacs``, ``ambertools``, ``densification``,
``precure``, ``CURE``, and ``postcure`` blocks follow the same structure
as :ref:`example 1 <ps_configuration>`.  The two interesting items here
are ``constituents`` (specifically the ``symmetry_equivalent_atoms`` and
``stereocenters`` keys) and what those mean for the template set
``htpolynet`` ends up generating.

``constituents``
^^^^^^^^^^^^^^^^

.. code-block:: yaml

   constituents:
     DGE:
       smiles: "CC(C)(c1ccc(OC[CH:3]([OH:5])[CH3:1])cc1)c1ccc(OC[CH:4]([OH:6])[CH3:2])cc1"
       reactive_atoms: {1: C1, 2: C2, 3: C3, 4: C4, 5: O1, 6: O2}
       count: 1600
       symmetry_equivalent_atoms: [[C1, C2], [C3, C4], [O1, O2]]
       stereocenters: [C3]
     PAC:
       smiles: "C1C[CH:1](CCC1CC2CC[CH:2](CC2)[NH2:3])[NH2:4]"
       reactive_atoms: {1: C1, 2: C1, 3: N1, 4: N1}
       count: 800
       symmetry_equivalent_atoms: [[N1, N2], [C1, C2]]
       stereocenters: [C1]

A 2:1 stoichiometric ratio (200 DGE : 100 PAC) gives at most 400 C–N
bonds at full conversion.

Reaction expansion via ``symmetry_equivalent_atoms``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DGEBA has two chemically equivalent reactive carbons (C1, C2), two
equivalent stereo carbons (C3, C4), and two equivalent oxirane oxygens
(O1, O2); PACM has two equivalent amine nitrogens (N1, N2) and the two
carbons they're attached to (C1, C2).  Declaring these pairings via
``symmetry_equivalent_atoms`` tells ``htpolynet`` to **expand** each
user-supplied reaction into the full set of symmetry-equivalent
variants.  The two reactions in this YAML expand as follows:

* The single primary-to-secondary reaction ``PAC.N1 + DGE.C1`` expands
  into **four**: ``PAC.{N1,N2} × DGE.{C1,C2}``, with products
  ``PAC~N1-C1~DGE``, ``PAC~N1-C2~DGE``, ``PAC~N2-C1~DGE``,
  ``PAC~N2-C2~DGE``.
* Because each of those products is itself a reactant in the
  secondary-to-tertiary reaction, *that* reaction expands into **eight**
  variants — one per (primary-stage product) × (which DGEBA carbon
  reacts next).  See the diagnostic log at the start of a run for the
  full list.

So the user wrote 2 cure reactions and 1 cap reaction in the YAML, and
``htpolynet`` produces 12 cure templates + 2 cap templates (one per
DGEBA carbon).

Diastereomer expansion via ``stereocenters``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``stereocenters: [C3]`` (for DGE) and ``stereocenters: [C1]`` (for PAC)
look like a single chirality declaration each — but combined with the
symmetry-equivalence above, they expand:

* On DGEBA, declaring ``C3`` as chiral implicitly makes ``C4`` chiral
  too (they're symmetry-equivalent).  ``C1``/``C2`` are not chiral in
  the *active* form (each has two methyl ligands), but they become
  chiral after one of those methyls forms a C–N bond — so we treat them
  as *potential* stereocenters as well, giving 4 stereocenters per
  DGEBA → **16 diastereomers**.
* On PACM, ``C1`` chirality plus the C1↔C2 symmetry pairing gives
  2 stereocenters → **4 diastereomers**.

``htpolynet`` builds the entire diastereomer pool at startup by flipping
stereocenters, then samples uniformly from the pool when placing
molecules into the initial liquid.  This guarantees the starting state
is racemic; polymerization cannot then introduce a net handedness as an
artefact.

This dance — 2 user reactions → 12 cure templates + 2 cap templates,
and 1 + 1 stereocenter declarations → 16 + 4 = 20 diastereomers —
explains why the diagnostic log for this example reports 22 molecules
detected (5 explicit + 6 implied by stereochemistry + 11 implied by
symmetry).

Now we can turn to actually :ref:`running the build <pde_run>`.
