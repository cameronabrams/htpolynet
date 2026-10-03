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

``htpolynet`` parameterizes the two templates, BDC and BDC3 (three monomers, 105
atoms), discussed on the :ref:`reactions page <badcy_reactions>`.  With ``gas``
charges the whole setup stage takes seconds.  Because the only cure reaction closes a
ring, the log says that conversion is counted by the ring cure rather than quoting a
bond count::

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
   INFO> Iteration 1: 97 ring(s), 291 of 720 groups consumed (conversion 0.404)
   ...
   INFO> Iteration 21: 3 ring(s), 702 of 720 groups consumed (conversion 0.975)
   INFO> Ring cure reached conversion 0.975
   INFO> Ring cure settled after 234 ring(s); the network is relaxed before postcure rather than at 500 K
   INFO> ********** Ring cure ends: 234 ring(s), conversion 0.975 **********

On the reference build (one V100 GPU, 16 cores) the ring cure took 21 iterations and
55 minutes:

.. list-table::
   :header-rows: 1
   :widths: 20 20 30 30

   * - Iteration
     - Rings closed
     - Groups consumed (of 720)
     - Cyanate conversion
   * - 1
     - 97
     - 291
     - 0.404
   * - 2
     - 30
     - 381
     - 0.529
   * - 3
     - 17
     - 432
     - 0.600
   * - 4
     - 8
     - 456
     - 0.633
   * - 5
     - 7
     - 477
     - 0.662
   * - 6
     - 4
     - 489
     - 0.679
   * - 7
     - 4
     - 501
     - 0.696
   * - 8
     - 6
     - 519
     - 0.721
   * - 9
     - 7
     - 540
     - 0.750
   * - 10
     - 6
     - 558
     - 0.775
   * - 11
     - 5
     - 573
     - 0.796
   * - 12
     - 5
     - 588
     - 0.817
   * - 13
     - 6
     - 606
     - 0.842
   * - 14
     - 7
     - 627
     - 0.871
   * - 15
     - 4
     - 639
     - 0.887
   * - 16
     - 3
     - 648
     - 0.900
   * - 17
     - 6
     - 666
     - 0.925
   * - 18
     - 4
     - 678
     - 0.942
   * - 19
     - 2
     - 684
     - 0.950
   * - 20
     - 3
     - 693
     - 0.963
   * - 21
     - 3
     - 702
     - 0.975

It stopped at iteration 21 because 702 groups is the first count at or above
``desired_conversion: 0.97``.  Six candidate rings were declined along the way, all in
iterations 16-21, each left 0.305-0.319 nm from closing by the closure ladder against
the 0.300 nm limit; their groups were offered again later.  None was declined for
threading.

Expect the familiar long tail.  Early iterations close dozens of rings at once while
groups are plentiful and close together; late ones close a handful, because each
remaining group must find two partners within reach that are not already used up,
and the ring cure must pull three separate molecules into place inside a network that
is already crosslinked.

Postcure
^^^^^^^^

Postcure runs the standard anneal -- two cycles between 300 and 500 K -- followed by a
100 ps NPT postequilibration at 300 K and 1 bar.

The reference build ended at 1149.5 kg/m³.  The :ref:`results page <badcy_results>`
has the stage profile, and the :ref:`postsim page <badcy_postsim>` explains why that
density is not the one to quote.

Next is the :ref:`results page <badcy_results>`.
