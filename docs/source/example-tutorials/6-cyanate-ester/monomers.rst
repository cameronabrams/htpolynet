.. _badcy_monomers:

Monomers
--------

One constituent:

* **BDC** -- bisphenol-A dicyanate, N#C-O-C\ :sub:`6`\ H\ :sub:`4`-C(CH\ :sub:`3`)\ :sub:`2`-C\ :sub:`6`\ H\ :sub:`4`-O-C#N.
  The bisphenol-A core of :ref:`example 2 <bgs_tutorial>` with each phenolic
  hydrogen replaced by a cyanate group.

BDC
^^^

.. admonition:: Placeholder
   :class: caution

   **TODO:** insert a 2D structure render at ``pics/BDC.png``.

The atom-mapped SMILES:

.. code-block:: yaml

   BDC:
     smiles: "[N:1]#[C:2]Oc1ccc(cc1)C(C)(C)c1ccc(O[C:3]#[N:4])cc1"
     reactive_atoms: {1: N1, 2: C1, 3: C2, 4: N2}
     symmetry_equivalent_atoms: [[C1, C2], [N1, N2]]
     count: 360

* ``N1#C1`` and ``C2#N2`` are the two cyanate groups.  **Both atoms of each group are
  reactive**: around a triazine, the carbon of one group bonds to the nitrogen of the
  next, so a group contributes its C to one ring bond and its N to another.
* ``symmetry_equivalent_atoms`` declares the two groups chemically equivalent, so the
  one reaction on the :ref:`next page <badcy_reactions>`, which is written against
  ``N1`` and ``C1``, applies to either group of any monomer.
* There are **no sacrificial hydrogens**.  Cyclotrimerization is an addition: the
  triazine's atoms are the cyanate groups' own C and N atoms, and nothing leaves.

Counting groups
^^^^^^^^^^^^^^^

360 BDC × 2 cyanate groups = 720 groups, so full conversion would be 240 triazine
rings and 720 new bonds.  The build does not reach full conversion -- no real cure
does -- and the groups left over stay as intact -O-C#N end groups, which need no
further treatment.  The diagnostic log reports the group count when the ring cure
starts::

   Ring cure begins: 720 reactive BDC group(s), (('N1', 'C1'), ('N2', 'C2')), template BDC3

The next page walks through the :ref:`reaction <badcy_reactions>`.
