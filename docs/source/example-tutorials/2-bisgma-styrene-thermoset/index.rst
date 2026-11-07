.. _bgs_tutorial:

BisGMA-Styrene Thermoset
========================

This tutorial covers a free-radical cure of bisGMA with styrene comonomers.
What makes it interesting beyond :ref:`example 1 <ps_tutorial>` is that
``htpolynet`` assembles the bisGMA molecule itself at runtime by reacting
2-hydroxypropyl isopropyl ester (HIE) onto each hydroxyl of bisphenol A
(BPA) via two **param**-stage reactions, then uses the assembled GMA
(and free STY) as the building blocks for the cure.  It also exercises
stereochemistry handling: HIE has two chiral carbons, so each GMA in the
initial liquid is drawn from a racemic pool of 16 diastereomers.

.. note::

   The configuration was enlarged in November 2026 for showcasing the analysis: it
   doubles the box side, 8 times the monomers of the earlier one.  Recorded output, timings and
   measured values on these pages are from the earlier, smaller builds until they are
   re-recorded from builds at the new size.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction
   reactions
   configuration
   run
   results
   postsim
