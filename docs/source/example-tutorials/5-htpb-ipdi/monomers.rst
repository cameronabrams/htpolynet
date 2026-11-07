.. _htpb_monomers:

Monomers
--------

Six constituents:

* Four small **building blocks** (SMILES-defined):

  * ``OB`` — 1-butene
  * ``TB`` — trans-2-butene
  * ``TBO`` — 1-hydroxy-trans-2-butene (the chain terminator)
  * ``IPD`` — isophorone diisocyanate, in its urethane-active form

* Two long-chain **assembled monomers** (built from the small blocks
  by param- and build-stage reactions):

  * ``DHT`` — linear HTPB chain
  * ``THT`` — three-arm branched HTPB chain

OB: 1-butene
^^^^^^^^^^^^

A four-carbon vinyl monomer; the vinyl C=C provides the bond-forming
site for chain extension.

.. code-block:: yaml

   OB:
     smiles: "[CH2:1]=[CH:2][CH2:3][CH3:4]"
     reactive_atoms: {1: C1, 2: C2, 3: C3, 4: C4}

All four carbons are named explicitly because the chain-assembly
reactions reference both ``C3`` (in the ``OBT`` dimer) and ``C4``
(propagated through ``OB3``).  ``OB`` acts as the chain initiator on
each HTPB chain — the chain assembly grows out of one ``OB`` per
chain plus a long tail of ``TB`` segments and a ``TBO`` end-cap.

TB: trans-2-butene
^^^^^^^^^^^^^^^^^^

The repeat unit.  HTPB chains are mostly polymerized trans-2-butene.

.. code-block:: yaml

   TB:
     smiles: "[CH3:1]/[CH:2]=[CH:3]/[CH3:4]"
     reactive_atoms: {1: C1, 2: C2, 3: C3, 4: C4}

``C1`` and ``C4`` are the two methyl carbons that lose a sacrificial H
to extend the chain; ``C2`` and ``C3`` are the internal vinyl carbons
that retain their double bond after extension.  In the assembled HTPB
chain, the conjugated C=C bonds remain — HTPB is famously elastomeric
because of this internal unsaturation.

TBO: 1-hydroxy-trans-2-butene
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The chain terminator: a 2-butene unit with one methyl replaced by a
hydroxymethyl.  Provides the reactive ``-OH`` group at the chain end
that bonds to IPDI during cure.

.. code-block:: yaml

   TBO:
     smiles: "[OH:5][CH2:1]/[CH:2]=[CH:3]/[CH3:4]"
     reactive_atoms: {1: C1, 2: C2, 3: C3, 4: C4, 5: O1}

``C4`` is the chain-extension site (bonds to a ``TB.C4`` during
build); ``O1`` is the cure site (bonds to an ``IPD.C1`` or
``IPD.C2`` during cure).  Each HTPB chain ends in a ``TBO`` residue
on each terminus, so linear ``DHT`` has 2 hydroxyls and branched
``THT`` has 3.

IPD: isophorone diisocyanate (active form)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Isophorone diisocyanate
<https://en.wikipedia.org/wiki/Isophorone_diisocyanate>`_ is a
cyclohexyl ring carrying two non-equivalent ``-N=C=O`` isocyanate
groups.  htpolynet needs an active form where each reactive carbon
carries a sacrificial H, so the YAML draws each ``-N=C=O`` as the
hydrated/formamide ``-N(H)-C(H)=O``:

.. code-block:: yaml

   IPD:
     smiles: "CC1(CN[CH:2]=O)C[CH:3](N[CH:1]=O)CC(C)(C)C1"
     reactive_atoms: {1: C1, 2: C2, 3: C3}
     stereocenters:
       - C3
     count: 250

* ``C1, C2`` — the two formyl carbons.  Each loses its sacrificial H
  when forming a urethane bond with a ``TBO.O1``.  The two are
  distinguished because IPDI's chemistry makes the primary (``C1``,
  on the secondary carbon ``C3``) and secondary (``C2``, on the
  methylene of the ring substituent) isocyanates non-equivalent.  The
  YAML therefore declares ``urethane-1`` (TBO.O1 + IPD.C1) and
  ``urethane-2`` (TBO.O1 + IPD.C2) as two independent cure reactions
  — see :ref:`htpb_reactions`.
* ``C3`` — the ring stereocenter.  Declared in ``stereocenters:``
  so the initial liquid samples both diastereomers uniformly; the
  one user-declared center expands to a 2-isomer pool that htpolynet
  draws from when placing the 250 IPD molecules in the initial box.

DHT and THT: assembled HTPB chains
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``DHT`` (linear) and ``THT`` (three-arm branched) are not declared
with a SMILES.  htpolynet assembles them at runtime from the four
small constituents using a sequence of param- and build-stage
reactions, detailed on the :ref:`reactions page <htpb_reactions>`.
The constituent declarations only specify counts and conformer-
generation parameters:

.. code-block:: yaml

   DHT:
     count: 100
     conformers:
       count: 6
       generator:
         name: gromacs
         params:
           ensemble: nvt
           temperature: 900
           ps: 100
           pad: 10.0
   THT:
     count: 100
     conformers:
       count: 6
       ...

The ``conformers`` block tells htpolynet to generate **6 thermalized
conformers per chain** at 900 K via GROMACS MD before placement, so
the initial liquid has chain-shape diversity rather than 100 copies of
a single rigid conformer.  This is a meaningful cost (~10 min per
chain × 100 chains × 6 conformers) but is essential for getting a
plausible starting state of a polymer melt.

Once assembled, each ``DHT`` chain has ~960 atoms and each ``THT``
chain has ~1100 atoms; combined with 250 IPDs, the initial system
totals roughly 111,000 atoms.

The next page walks through the :ref:`param/build/cure reactions
<htpb_reactions>` that assemble the chains and cure the urethane
network.
