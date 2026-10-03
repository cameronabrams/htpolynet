.. _badcy_introduction:

Introduction
------------

Set up a clean working directory and pull the example YAML:

.. code-block:: console

   $ mkdir my_badcy
   $ cd my_badcy
   $ htpolynet fetch-example 6
   Fetched 6-cyanate-ester.yaml  (run with: htpolynet run 6-cyanate-ester.yaml)
   $ ls
   6-cyanate-ester.yaml

Self-contained YAML as in the earlier examples: the one monomer is generated from an
atom-mapped SMILES string at startup.

What's new in this example
^^^^^^^^^^^^^^^^^^^^^^^^^^

**The cure closes rings.**  BADCy thermosets cure by `cyclotrimerization
<https://en.wikipedia.org/wiki/Cyanate_ester>`_: three cyanate groups come together
and close into a 1,3,5-triazine ring with three aryl-ether arms.  No atom is lost.
Every earlier example cures by forming one bond at a time between two reactive atoms,
which ``CURE`` does by searching for close pairs.  A triazine cannot be built that way,
because its three bonds only make sense together.  The ``ring_cure`` block replaces
``CURE`` for a reaction like this one, and each iteration of it

1. searches for **triples** of unreacted cyanate groups close enough to close a ring,
   and chooses a set that shares no group;
2. pulls each chosen triple together through a ladder of harmonic restraints, the
   three-body analogue of ``CURE``'s drag;
3. **declines** any ring the ladder could not pull shut, and any ring that would close
   around an existing bond or leave a monomer's bond threaded through an existing ring
   -- a threaded ring is a permanent topological defect, so it is refused rather than
   repaired;
4. closes the accepted rings, re-types and re-charges their atoms from the trimer
   template, and relaxes and re-equilibrates the box, as ``CURE`` does after each batch
   of bonds.

When the target conversion is reached, the network is settled once -- a minimization
and a short 300 K NVT run -- before postcure.

**Conversion means cyanate conversion.**  Three groups per ring, counted directly.
Fang and Shimp's end-group arithmetic caps *acyclic* growth at 75 % cyanate
conversion, so a network that reaches the 97 % this example asks for necessarily
contains macrocycles.  That is the real chemistry, and it is why no cycle-length
filter is applied.

What you'll see in the build
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. admonition:: Placeholder
   :class: caution

   **TODO:** iteration count, ring count, final conversion and wall time from the
   reference build of the 2.15.0 configuration.

The remaining pages walk through the monomer SMILES, the cure reaction, the YAML in
full, and what to look for in the diagnostic log and plots.
