.. _dfe_monomers:

Monomers
--------

FDE
^^^

FDE — *furfurylmethylenediepoxide*, name ``FDE`` — is the diepoxide
hardener.  Its skeleton is a tertiary nitrogen with three arms:

* one furfuryl arm (a furan ring with a methylene bridge to the N);
* two opened-oxirane arms (each a ``-CH2-CH(OH)-CH3`` chain).

The native, *closed* oxirane form would have each arm end in a
three-membered ``C-O-C`` ring.  As in :ref:`example 3 <pde_monomers>`,
``htpolynet`` expects the **active** form with the oxirane already
opened into ``-CH(OH)-CH3``, so each reactive carbon carries a
sacrificial H that gets removed when a new C–N bond forms.  The YAML's
``constituents.FDE`` block carries that active-form SMILES with
atom-mapping tokens:

.. code-block:: yaml

   FDE:
     smiles: "N(Cc1occc1)(C[CH:3]([OH:5])[CH3:1])C[CH:4]([OH:6])[CH3:2]"
     reactive_atoms: {1: C1, 2: C2, 3: C3, 4: C4, 5: O1, 6: O2}
     count: 1600
     symmetry_equivalent_atoms: [[C1, C2], [C3, C4], [O1, O2]]
     stereocenters: [C3]

The six named atoms map one-for-one onto the DGE atoms from
:ref:`tutorial 3 <pde_monomers>`:

* **C1, C2** — the two reactive (sacrificial-methyl) carbons.
* **C3, C4** — the two stereocenters (the carbons bearing the
  hydroxyl).
* **O1, O2** — the two hydroxyl oxygens, used by the cap reaction to
  re-form the oxirane ring on any unreacted carbon.

The ``symmetry_equivalent_atoms`` and ``stereocenters`` declarations
have the same effect as in DGEBA: the two opened-oxirane arms are
treated as chemically equivalent, and the diastereomer pool is built
from flipping the chiral carbons.

DFA
^^^

DFA — *difurfuryl diamine* — is the diamine that plays the role of
PACM.  Two furfurylmethylamine units are joined by a central methylene
bridge.  Each amine is primary, so each can react twice — a single DFA
can crosslink up to four FDE molecules through its two nitrogens.

.. code-block:: yaml

   DFA:
     smiles: "C(c1oc([CH2][NH2:1])cc1)c1oc([CH2][NH2:2])cc1"
     reactive_atoms: {1: N1, 2: N2}
     count: 800
     symmetry_equivalent_atoms: [[N1, N2]]

Only the two nitrogens need explicit names; the cure reactions
reference them as ``N1`` (and via the symmetry pairing, ``N2``).  No
stereocenters: the furan rings are planar aromatic and the central
methylene is not chiral.

A note on the SMILES
^^^^^^^^^^^^^^^^^^^^

For the FDE arms, ``[CH:3]`` and ``[CH:4]`` are written with an
explicit ``H`` rather than the bare ``[C:3]`` / ``[C:4]``.  This is the
same valence-conservation gotcha noted in :ref:`example 2
<bgs_configuration>`: a SMILES bracket atom without an explicit H
count defaults to **zero** implicit Hs, which would leave the carbon
at valence 3, and ``antechamber`` would then mistype it (and the
adjacent atoms) so badly that ``tleap`` fails.  Writing ``[CH:3]``
keeps the carbon at valence 4 and AmberTools handles it cleanly.

Historically, this example was awkward because openbabel's mol2 writer
produced inconsistent atom numbering for fused heterocycles, so the
original workflow used PDB inputs with manual atom renaming via
``sed``.  Atom-mapped SMILES sidesteps that entirely: the mapping
tokens encode which atom gets which name, so the mol2 atom-ordering
never matters.

The next page describes the :ref:`reaction dictionaries <dfe_reactions>`.
