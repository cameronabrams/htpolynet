.. _badcy_reactions:

Reactions
---------

Two reactions appear in the YAML:

1. ``etherify`` — a **cure-stage** reaction.  Forms a single
   ``BPA-O—C(triazine)`` aryl ether bond.  Each side loses one
   hydrogen.
2. ``cap_with_cyanate`` — a **repair-stage** reaction.  Never executed
   as a normal reaction; ``htpolynet`` uses it at setup time to
   parameterize the ``BPA~O1-C1~CYN`` linked-product template that the
   postcure repair driver then splices into the system for every cap
   it forms.

The cure reaction
^^^^^^^^^^^^^^^^^

.. code-block:: yaml

   - name: etherify
     stage: cure
     reactants: {1: BPA, 2: TAZ}
     product: BPA~O1-C1~TAZ
     probability: 1.0
     atoms:
       A: {reactant: 1, resid: 1, atom: O1, z: 1}
       B: {reactant: 2, resid: 1, atom: C1, z: 1}
     bonds:
       - atoms: [A, B]
         order: 1

Read this as "an ``O1`` on BPA with one remaining crosslink site
(``z: 1``) bonds to a ``C1`` on TAZ with one remaining crosslink site;
one sacrificial H disappears on each side; the resulting linked
product is ``BPA~O1-C1~TAZ``."

The symmetry-expander applies both monomers'
``symmetry_equivalent_atoms`` groups to this single reaction and
produces six (O, C) variants — ``BPA.{O1,O2} × TAZ.{C1,C2,C3}`` — at
setup time, each as a separately parameterized linked-product
template.  The cure machinery treats them as six distinct reactions in
the iterative bond search, but the user only had to write one.

Six templates are not enough to give the ring the right charges,
though.  Each describes a triazine with one BPA attached, so it shows
the other two ring carbons unreacted, and a template splice rewrites
the charges of every atom within two bonds of the new bond, which
includes those carbons.  Without more context, the second and third
bonds on a ring would reset the carbons that reacted before them to
their unreacted charge.  So htpolynet also builds, at setup time, one
template for each bonding carbon, BPA oxygen, and set of ring carbons
already bonded: 18 more, with names like
``BPA~O1-C3~BPA~O1-C2~BPA~O1-C1~TAZ``.  At each bond, the cure picks
the one that matches the ring as it stands.  Nothing in the
configuration asks for these.  htpolynet generates them for any
reactive atom whose ``symmetry_equivalent_atoms`` partners sit two
bonds away, as a triazine's ring carbons do.

What the cure stage produces, then, is an A2+B3 step-growth ether
network: every formed bond is the same chemistry (an aryl O to a
triazine ring C, replacing one H on each side with a new C-O bond),
and triazine rings act as trifunctional crosslink nodes that exist in
the monomer before cure starts.

The repair-stage reaction
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

   - name: cap_with_cyanate
     stage: repair
     reactants: {1: BPA, 2: CYN}
     product: BPA~O1-C1~CYN
     probability: 1.0
     atoms:
       A: {reactant: 1, resid: 1, atom: O1, z: 1}
       B: {reactant: 2, resid: 1, atom: C1, z: 1}
     bonds:
       - atoms: [A, B]
         order: 1

This looks identical to ``etherify`` save for two things:

* ``stage: repair``.  The new ``repair`` stage value (introduced in
  htpolynet 2.1) tells the runtime that this reaction is **never
  executed** during the cure-loop or capping pass; it exists only so
  the same setup-time machinery that builds the cure-product
  templates also builds the ``BPA~O1-C1~CYN`` linked-product
  template.  See the :ref:`postcure-repair user-guide page
  <postcure_repair>` for the architecture details.

* ``reactant 2`` is ``CYN`` (hydrogen cyanide), not ``TAZ``.  The
  ``BPA~O1-C1~CYN`` product is a single BPA residue bonded through one
  of its phenolic oxygens to a ``-C#N`` group — the BADCy residual
  end-group we want every unreacted BPA-OH to end up with after
  repair.

At setup time ``htpolynet`` parameterizes 11 templates in total: BPA,
TAZ, CYN (the three constituents), the six ``BPA~O1-Cn~TAZ`` variants
from symmetry-expanding ``etherify``, and the two ``BPA~On-C1~CYN``
variants from symmetry-expanding ``cap_with_cyanate`` (only O1 and O2
of BPA expand because CYN has no symmetry-equivalent atoms).

The diagnostic-log report at startup makes this concrete:

.. code-block:: text

   INFO> 11 molecules detected in 6-cyanate-ester.yaml
   INFO>                       explicit: 5    # BPA, TAZ, CYN, BPA~O1-C1~TAZ, BPA~O1-C1~CYN
   INFO>     implied by stereochemistry: 0
   INFO>            implied by symmetry: 6    # 5 more from etherify's symmetry; 1 more from cap_with_cyanate

The repair driver
^^^^^^^^^^^^^^^^^

The actual repair operation — find incomplete triazines, dismantle
them, attach the resulting fragments — is driven by a
``postcure_repair`` block in the YAML, **not** by the
``cap_with_cyanate`` reaction directly.  The repair block configures a
``triazine_to_cyanate_cap`` driver that runs after cure completes:

.. code-block:: yaml

   postcure_repair:
     - type: triazine_to_cyanate_cap
       crosslinker:
         residue: TAZ
         ring_carbon_atoms: [C1, C2, C3]
         ring_nitrogen_atoms: [N1, N2, N3]
         full_bond_count: 3
       bridge:
         residue: BPA
         reactive_oxygen_atoms: [O1, O2]
       cap_residue: CYN
       cap_template: BPA~O1-C1~CYN
       cap_search_radius: 0.6   # nm

The fields are exactly what they look like:

* ``crosslinker`` — which residue and which atoms make up the ring to
  be dismantled.  ``full_bond_count: 3`` is the threshold below which
  a ring is considered incomplete (and thus dismantled).
* ``bridge`` — which residue and atoms carry the reactive O that
  bonds to a ring C during cure.  Free bridge-Os (still carrying their
  H after cure) are the recipients of donated free-cap fragments.
* ``cap_residue`` / ``cap_template`` — the residue name to assign to
  each newly-formed cap, and the parameterization template that
  supplies its atom types, charges, and bonded interactions.  The
  ``cap_template`` is the product name of the ``cap_with_cyanate``
  repair-stage reaction above.
* ``cap_search_radius`` — the radius (nm) within which a free-cap
  fragment looks for an unreacted bridge-O to attach to.  Expanded up
  to 10× on miss; falls back to globally nearest as a last resort.

The :ref:`postcure-repair user-guide page <postcure_repair>` walks
through the dismantle-and-donate algorithm in detail, including the
atom-conservation argument and how the within-ring C-N matching is
chosen.  The next page covers the :ref:`full YAML
<badcy_configuration>`.
