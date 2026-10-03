.. _badcy_reactions:

Reactions
---------

One reaction appears in the YAML:

.. code-block:: yaml

   - name: cyclotrimerize
     stage: cure
     reactants: {1: BDC, 2: BDC, 3: BDC}
     product: BDC3
     probability: 1.0
     atoms:
       N1: {reactant: 1, resid: 1, atom: N1, z: 1}
       C1: {reactant: 1, resid: 1, atom: C1, z: 1}
       N2: {reactant: 2, resid: 1, atom: N1, z: 1}
       C2: {reactant: 2, resid: 1, atom: C1, z: 1}
       N3: {reactant: 3, resid: 1, atom: N1, z: 1}
       C3: {reactant: 3, resid: 1, atom: C1, z: 1}
     bonds:
       - atoms: [N1, C2]
         order: 1
         sacrificial_h: false
       - atoms: [N2, C3]
         order: 1
         sacrificial_h: false
       - atoms: [N3, C1]
         order: 1
         sacrificial_h: false

Read it as: three BDC molecules each contribute one cyanate group, and the nitrogen of
each group bonds to the carbon of the next -- N1-C2, N2-C3, N3-C1 -- which closes the
six-membered C\ :sub:`3`\ N\ :sub:`3` ring.  The keys ``N1`` ... ``C3`` on the left are
labels local to the reaction; the ``atom:`` values are the monomer's atom names.

A few things to notice:

* **Three reactants and three bonds that form a cycle among them.**  htpolynet
  recognizes a cure reaction whose inter-reactant bonds close a loop -- here 1-2, 2-3,
  3-1 -- as *ring-closing*, and hands it to the ring cure.  ``CURE``'s pairwise search
  skips it.
* **``sacrificial_h: false`` on every bond.**  Cyclotrimerization is an addition; the
  ring is made of the cyanate groups' own atoms and no hydrogen is removed.
* **One product template, ``BDC3``.**  htpolynet builds the trimer once at setup:
  three BDC molecules with their rings closed.  Before it computes the trimer's charges
  it closes the ring under restraints, so antechamber sees a triazine with ordinary
  bond lengths rather than one stretched across the gap the reactants started at.  The
  setup log accordingly reports two templates, BDC and BDC3.
* **Either cyanate group can react.**  The reaction names only ``N1`` and ``C1``, but
  the ring cure finds every unreacted group -- ``(N1, C1)`` or ``(N2, C2)`` -- and maps
  whichever group of each monomer took part onto the trimer template by name, using the
  ``symmetry_equivalent_atoms`` declared on the :ref:`monomers page <badcy_monomers>`.

How the ring cure uses it
^^^^^^^^^^^^^^^^^^^^^^^^^

Each iteration of the ring cure

1. enumerates candidate triples: groups whose C-to-N gaps around a would-be ring are
   all within ``search_radius``, scored by the total length of the three new bonds;
2. chooses greedily, shortest first, a set of triples that share no group -- by default
   no two groups from one molecule either (``same_molecule: false``);
3. pulls the chosen triples together through ``closure.nstages`` stages of harmonic
   restraints, each followed by a short minimization and MD;
4. declines a triple whose gaps are still longer than ``closure.max_accept`` after the
   ladder, or whose new ring would thread an existing bond or leave a dragged
   monomer's bond threading an existing ring -- the groups stay unreacted and can be
   offered again;
5. closes the accepted rings, applies the ``BDC3`` template's types and charges to the
   atoms of the three residues, registers each new triazine as a ring so later
   threading checks see it, and runs the relax and equilibration stages.

When an iteration finds too few rings it widens the search by ``radial_increment``.
The cure ends at ``desired_conversion``, after ``max_iterations``, or when an iteration
finds nothing even at the largest search radius.  The next page shows the :ref:`full YAML <badcy_configuration>`.
