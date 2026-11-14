.. _bgs_reactions:

Reactions
---------

The YAML has three groups of reactions: two ``param``-stage reactions that
assemble bisGMA from BPA + 2×HIE at runtime, four ``cure``-stage reactions
that polymerize STY and HIE radicals, and two ``cap``-stage reactions,
written by ``htpolynet`` itself, that restore vinyl double bonds on any
monomer that did not react.

Building bisGMA at runtime
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::

    * - .. figure:: pics/BPA.png

           Bisphenol-A (BPA)

      - .. figure:: pics/HIE.png

           HIE (glycidyl methacrylate), active form

      - .. figure:: pics/GMA.png

           bisGMA (active form)

The first reaction esterifies one HIE onto one of BPA's two phenolic
hydroxyls:

.. code-block:: yaml

  - name: B1
    stage: param
    reactants:
      1: BPA
      2: HIE
    product: GM1
    atoms:
      A:
        reactant: 1
        resid: 1
        atom: O1
        z: 1
      B:
        reactant: 2
        resid: 1
        atom: C4
        z: 2
    bonds:
      - atoms: [A, B]
        order: 1

The intermediate ``GM1`` then reacts with a second HIE on the other phenolic
hydroxyl:

.. code-block:: yaml

  - name: B2
    stage: param
    reactants:
      1: GM1
      2: HIE
    product: GMA
    atoms:
      A:
        reactant: 1
        resid: 1
        atom: O2
        z: 1
      B:
        reactant: 2
        resid: 1
        atom: C4
        z: 1
    bonds:
      - atoms: [A, B]
        order: 1

Both reactions are ``stage: param`` because the product needs to be
GAFF-parameterized.  Once both run, the system has a fully parameterized
``GMA`` template ready to be inserted into the initial liquid.

Why build GMA at runtime rather than supply a pre-made ``GMA.mol2``?  The
reason is that the two HIE moieties in a GMA each house a reactive C=C
double bond, and we want ``htpolynet`` to treat them as *equivalent* HIE
sites during cure.  By keeping HIE as the canonical reactive species in
the ``cure`` reactions, the chain-expanded template set is much smaller
than it would be if GMA itself appeared on either side of a cure bond.

Cure reactions
^^^^^^^^^^^^^^

Cure proceeds by C=C double-bond opening: ``C1`` of one radical bonds to
``C2`` of another, with each side donating a sacrificial H.  On both STY
and HIE, ``C1`` is the radical carbon and ``C2`` is the terminal methyl.
We need to cover all four reactant combinations:

.. table:: Monomer–monomer reactions in the bisGMA/styrene system
    :widths: auto

    ===================  =====================
    Attacker (C1 owner)  Attackee (C2 owner)
    ===================  =====================
    HIE                  HIE
    STY                  STY
    STY                  HIE
    HIE                  STY
    ===================  =====================

Each row gets a single ``C1``-attacks-``C2`` cure reaction in the YAML
(four reactions total; see ``2-bisgma-styrene-thermoset.yaml`` for the
full block — they all have the same structure as the polystyrene example,
just with the appropriate reactant pair).

Cap reactions
^^^^^^^^^^^^^

Any STY or HIE monomer that didn't react during CURE still carries its
saturated active form.  Cap reactions convert the leftover ``C1–C2``
single bonds back to the natural double bond.  The configuration does
not contain them: ``htpolynet`` writes one for each bond an ``activate``
key opens whose two atoms are both named, so STY's ``[[1, 2]]`` and HIE's
``[2, 1]`` each give one, equivalent to these.  HIE's other pair,
``[5, 4]``, gives none, because its O is unnamed, and needs none: every
HIE epoxide is opened building GMA.

.. code-block:: yaml

  - name:         STY-restore-C1-C2
    stage:        cap
    reactants:    {1: STY}
    product:      STY-capped-C1-C2
    probability:  1.0
    atoms:
      A: {reactant: 1, resid: 1, atom: C1, z: 1}
      B: {reactant: 1, resid: 1, atom: C2, z: 1}
    bonds:
      - atoms: [A, B]
        order: 2
  - name:         HIE-restore-C2-C1
    stage:        cap
    reactants:    {1: HIE}
    product:      HIE-capped-C2-C1
    probability:  1.0
    atoms:
      A: {reactant: 1, resid: 1, atom: C2, z: 1}
      B: {reactant: 1, resid: 1, atom: C1, z: 1}
    bonds:
      - atoms: [A, B]
        order: 2

Chain expansion
^^^^^^^^^^^^^^^

``htpolynet`` automatically generates the trimer and tetramer templates
needed to cover the dihedral environments that arise as polymer chains
grow.  Because the four cure reactions span 2 monomer types × 2 attacker
roles, there are 8 distinct trimer sequences, with the templated bond
in one of two positions in each sequence → **16 trimer templates**.  For
tetramers, only the middle (second–third) bond position is needed (the
other positions are already covered by the trimer templates), giving
**16 tetramer templates**.  Total chain-expanded cure templates: **32**.

You do not list these explicitly; ``htpolynet`` reports the expanded set
in the diagnostic log at the start of a run.

The next page walks through the full :ref:`configuration file <bgs_configuration>`.
