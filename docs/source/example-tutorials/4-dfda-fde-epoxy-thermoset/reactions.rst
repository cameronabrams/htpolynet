.. _dfe_reactions:

Reactions
---------

Because the cure chemistry here is identical to
:ref:`example 3 <pde_reactions>` — epoxy + amine step-growth — the
reactions block is a near-rename of tutorial 3's: ``DFA`` plays the role
of ``PAC``, ``FDE`` plays the role of ``DGE``, and the bond label
``N1-C1`` is the same.

The three directives in this YAML are:

Primary-to-secondary amine
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

   - name:        'Primary-to-secondary-amine'
     stage:       cure
     reactants:
       1: DFA
       2: FDE
     product:     DFA~N1-C1~FDE
     probability: 1.0
     atoms:
       A: {reactant: 1, resid: 1, atom: N1, z: 2}
       B: {reactant: 2, resid: 1, atom: C1, z: 1}
     bonds:
       - atoms: [A, B]
         order: 1

``z: 2`` on the nitrogen says it must still be primary (two H's
available); ``z: 1`` on the carbon says it has one sacrificial methyl
H to donate.

Secondary-to-tertiary amine
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

   - name:        'Secondary-to-tertiary-amine'
     reactants:
       1: DFA~N1-C1~FDE
       2: FDE
     product:     DFA~N1-C1~FDE-C1~FDE
     stage:       cure
     probability: 0.5
     atoms:
       A: {reactant: 1, resid: 1, atom: N1, z: 1}
       B: {reactant: 2, resid: 1, atom: C1, z: 1}
     bonds:
       - atoms: [A, B]
         order: 1

The ``probability: 0.5`` is the same intrinsic-rate weighting used in
example 3 — secondary-to-tertiary amine formation is empirically
slower than primary-to-secondary.

Oxirane-formation cap
^^^^^^^^^^^^^^^^^^^^^

The configuration contains no cap reaction; ``htpolynet`` writes this
one from FDE's ``activate`` key (one for the pair of symmetric
epoxides, whose copy it generates as it does any symmetric reaction):

.. code-block:: yaml

   - name:        FDE-restore-O1-C1
     reactants:
       1: FDE
     product:     FDE-capped-O1-C1
     stage:       cap
     probability: 1.0
     atoms:
       A: {reactant: 1, resid: 1, atom: O1, z: 1}
       B: {reactant: 1, resid: 1, atom: C1, z: 1}
     bonds:
       - atoms: [A, B]
         order: 1

Same intra-monomer ring-closure as in example 3: any FDE that still
has a sacrificial-methyl H at the end of CURE gets its three-membered
oxirane ring re-formed by bonding the hydroxyl oxygen ``O1`` to the
reactive carbon ``C1``.

How these expand
^^^^^^^^^^^^^^^^

Through the ``symmetry_equivalent_atoms`` declarations on
:ref:`the constituents <dfe_configuration>`, these three reactions
expand into:

* 4 primary-to-secondary cure templates (``FDE.{C1,C2} × DFA.{N1,N2}``);
* 8 secondary-to-tertiary cure templates (one per primary-stage product
  × which FDE carbon reacts next);
* 2 cap templates (one per FDE carbon).

The diagnostic log at the start of a run lists the full expanded set.

The next page covers the :ref:`configuration file <dfe_configuration>`.
