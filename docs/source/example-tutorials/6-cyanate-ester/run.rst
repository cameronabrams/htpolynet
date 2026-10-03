.. _badcy_run:

Running the Build
-----------------

From inside the working directory containing ``6-cyanate-ester.yaml``:

.. code-block:: console

   $ htpolynet run -diag diagnostics.log 6-cyanate-ester.yaml &> console.log &

The stage layout under ``proj-N/systems/`` matches the earlier examples: ``init/``,
``densification/``, ``precure/``, one ``iter-K/`` per ring-cure iteration,
``postcure/`` and ``final-results/``, plus ``plots/`` and ``profile.json`` at the
project root.  The ring cure also keeps its progress in ``systems/ring_state.yaml``,
which is what lets ``-restart`` resume it.

Setup
^^^^^

``htpolynet`` parameterizes the two templates, BDC and BDC3, discussed on the
:ref:`reactions page <badcy_reactions>`.  The trimer is the expensive one -- three
monomers, 105 atoms -- but it is built once and cached in the library, so later runs
skip it.  Because the only cure reaction closes a ring, the log says that conversion
is counted by the ring cure rather than quoting a bond count::

   INFO> Generated 2 molecule templates
   INFO> Initial composition is BDC 360
   INFO> Conversion is counted in reactive groups by the ring cure; see "Ring cure begins" for the total

Densification + precure
^^^^^^^^^^^^^^^^^^^^^^^

200 kg/m³ initial density and NPT segments at 300 K and 10 bar, extended until the
mean density is known to within 1 kg/m³, bring the box to liquid density.  Precure
runs a 200 ps NPT preequilibration at 300 K and 1 bar and two anneal cycles between
300 and 500 K to relax high-energy contacts from the random initial placement.

Ring cure
^^^^^^^^^

The log brackets the ring cure, reports each iteration's rings and the running group
count, and reports the settle::

   INFO>  Ring cure begins: 720 reactive BDC group(s), (('N1', 'C1'), ('N2', 'C2')), template BDC3
   INFO> Iteration 1: ... ring(s), ... of 720 groups consumed (conversion ...)
   ...
   INFO> Ring cure reached conversion ...
   INFO> Ring cure settled after ... ring(s); the network is relaxed before postcure rather than at 500 K
   INFO> ********** Ring cure ends: ... ring(s), conversion ... **********

.. admonition:: Placeholder
   :class: caution

   **TODO:** the per-iteration table (rings, cumulative conversion, wall time) and the
   declined-ring counts from the reference build of the 2.15.0 configuration.

Expect the familiar long tail.  Early iterations close dozens of rings at once while
groups are plentiful and close together; late ones close a handful, because each
remaining group must find two partners within reach that are not already used up,
and the ring cure must pull three separate molecules into place inside a network that
is already crosslinked.

Postcure
^^^^^^^^

Postcure runs the standard anneal -- two cycles between 300 and 500 K -- followed by a
100 ps NPT postequilibration at 300 K and 1 bar.

.. admonition:: Placeholder
   :class: caution

   **TODO:** final density and the end-of-run stage profile from the reference build.

Next is the :ref:`results page <badcy_results>`.
