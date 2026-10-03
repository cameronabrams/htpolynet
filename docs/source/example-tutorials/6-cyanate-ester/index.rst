.. _badcy_tutorial:

Bisphenol-A Dicyanate (BADCy) Thermoset
=======================================

This tutorial builds a **cyanate-ester thermoset** by the chemistry that actually
cures one: **cyclotrimerization**.  Three cyanate groups (-O-C#N) on three different
bisphenol-A dicyanate monomers close into one 1,3,5-triazine ring, and those rings
are the network's trifunctional crosslink junctions.  It is the only depot example
whose cure reaction has three reactants and closes a ring, so it is driven by a
``ring_cure`` block rather than ``CURE``.

Two things set it apart from the earlier examples:

* **A three-body cure reaction.**  ``CURE`` forms one bond at a time between two
  reactive atoms.  A triazine needs three bonds among three molecules at once, and a
  ring that is two-thirds closed is not a chemical species.  The ring cure therefore
  finds *triples* of cyanate groups, pulls each triple together under a ladder of
  restraints, and closes all three bonds of a ring in one step.

* **Conversion is the cyanate conversion.**  Each ring consumes three cyanate groups,
  and ``desired_conversion`` is the fraction of groups consumed -- the quantity FTIR
  measures and the literature reports.  A group that never finds two partners stays an
  intact -O-C#N end group, which is exactly what real under-cured BADCy contains.

.. note::

   Releases before 2.15.0 built this example with a stand-in model: bisphenol A plus
   a pre-formed 1,3,5-triazine, joined by aryl-ether bonds, with a
   :ref:`postcure repair stage <postcure_repair>` that dismantled the incomplete rings
   into -O-C#N caps afterward.  Its triazine rings existed before cure began, so the
   way the network formed -- and its gel point -- were not those of the chemistry.
   That model's conversion was also a bond conversion, which maps to the cyanate
   conversion only roughly as its cube.  The repair stage remains available for
   reproducing builds made that way.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction
   monomers
   reactions
   configuration
   run
   results
   postsim
