.. _liquid_styrene_monomer:

The Styrene Monomer
-------------------

Styrene is a vinyl monomer, and since this example does no curing it is
simulated as it really is, vinyl group and all (``c1ccccc1C=C``).

.. figure:: ../1-polystyrene/pics/STYCC.png
   :align: center
   :width: 240px

   Styrene.  ``C1`` and ``C2`` are the vinyl carbons.

Examples that cure styrene (:ref:`example 1 <ps_tutorial>` and beyond)
add ``activate: [[1, 2]]`` to the same entry.  That has ``htpolynet``
work with styrene's *active* form, ethylbenzene, whose vinyl carbons
each carry the sacrificial hydrogen a cure reaction strips when it forms
a new C–C bond.  Here there is no cure, so there is nothing to activate.

.. note::

   Through htpolynet 2.17, this example simulated ethylbenzene: it used
   the active form, which with no cure stayed saturated.  Its density is
   not the density of liquid styrene.

In-config SMILES generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``constituents`` block of ``0-liquid-styrene.yaml`` carries the
SMILES string and, via the SMILES atom-mapping syntax (``[CH:1]``,
``[CH2:2]``), names the two vinyl carbons, as example 1 does:

.. code-block:: yaml

   constituents:
     STY:
       smiles: "c1ccccc1[CH:1]=[CH2:2]"     # styrene
       reactive_atoms: {1: C1, 2: C2}
       count: 1000

.. warning::

   Run each example in its own directory.  Examples 0 and 1 both call
   their molecule ``STY``, but example 0's is styrene and example 1's is
   activated (ethylbenzene).  A generated ``lib/molecules/inputs/STY.mol2``
   is reused rather than regenerated, so in a shared directory the second
   example would build from the first one's molecule; ``htpolynet`` warns
   when that happens.

When ``htpolynet run`` starts, it sees the ``smiles`` key, invokes
RDKit to generate ``lib/molecules/inputs/STY.mol2`` with the atom-name
mapping (``:1`` → ``C1``, ``:2`` → ``C2``), and then proceeds with the
normal AmberTools parameterization workflow (``antechamber`` →
``parmchk2`` → ``tleap``).  Users without RDKit can instead supply a
pre-built mol2 in ``lib/molecules/inputs/`` and use the
``rename_atoms`` key to assign atom names — see
:ref:`molecular_structure_inputs`.
