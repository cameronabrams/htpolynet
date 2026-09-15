.. _configuration_run:

Configuration Files for ``htpolynet run``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``htpolynet run`` configuration file tells ``htpolynet`` what it needs in order to generate a polymerized system, beginning with structures of the individual monomers and a description of the polymerization chemistry.   ``htpolynet run`` expects at most ten distinct sections in a configuration file:

=================   =====================  
Section name        Role 
=================   =====================  
``Title``           Just a descriptive title
``gromacs``         Directives for interacting with Gromacs
``ambertools``      Directives for interacting with AmberTools
``GAFF``            Directives for handling inconsistencies in the General Amber Force Field
``constituents``    Directives for molecular constituents in the initial system
``densification``   Directives for how the densification phase is run
``precure``         Directives for how precure equilibration and annealing is run
``CURE``            Directives for how the CURE is run
``postcure``        Directives for how postcure equilibration and annealing is run
``reactions``       List of all reactions needed to build and cure the system
=================   =====================  

Sections can appear in any order (since the whole YAML file is like a nested python dictionary).  The ``gromacs`` and ``ambertools`` sections are mainly for specifying system-specific commands for running Gromacs and AmberTools executables; they have defaults that work for simple Linux workstations. ``densification``, ``precure``, and ``postcure`` all specify series of MD simulations on the system.  ``constituents`` specifies the initial make-up of the system, and ``reactions`` describe both how to build those constituents from input monomers (if necessary) as well as the types of bonds that you want to occur during polymerization.  Both ``constituents`` and ``reactions`` sections are where ``htpolynet`` extracts the names of molecular species in the system and how they are chemically interrelated.  ``CURE`` is the most complicated section, and it desribes how the CURE algorithm is to be run.  We consider each of these eight sections (minus ``Title``) below.

.. _configuration_details:

All the details
!!!!!!!!!!!!!!!

In this section we show all subdirectives for each of the five main directives in the configuration file.  Please note that, although default values are quoted for some of these parameters, these default values are *not* guaranteed to be optimal or even functional in all possible scenarios.

.. _gromacs_yaml_directive:

* ``gromacs``:  This directive specifies parameters ``htpolynet`` uses when invoking the Gromacs executable.
  
    =====================================    ==============  =====================
    ``gromacs`` subdirective                 Type            Description (default)
    =====================================    ==============  =====================
    ``gmx``                                  str             ``gmx`` command (default ``gmx``)
    ``gmx_options``                          quoted string   options to pass to all ``gmx`` calls (default ``-quiet -nobackup``)
    ``mdrun``                                quoted string   ``mdrun`` command (default ``gmx (options) mdrun``)
    ``mdrun_single_molecule``                quoted string   version of ``mdrun`` to use for any single-molecule Gromacs runs
    ``mdrun_options``                        dict            command-line arguments to pass to ``mdrun`` (none)
    =====================================    ==============  =====================

    If you are running on a supercomputer with a native installation of Gromacs, it is likely you should point the parameter ``gmx`` to the fully resolved pathname of ``gmx_mpi`` (or load the appropriate module), and use the ``mdrun`` parameters to specify the ``mpirun`` or ``mpiexec`` syntax needed to launch ``gmx_mpi mdrun``.  The ``gromacs_single_molecule`` subdirective allows you to specify a particular form of ``mdrun`` appropriate for single-molecule simulations.  These are most often used as part of parameterization or conformer generation.  Typically, it's best to run these on a single processor without domain decomposition.

    The ``gromacs`` directive is optional; if none is specified the default values are used.

* ``ambertools``:  This directive specifies parameters ``htpolynet`` uses when working with the AmberTools suite.

    =====================================    ==============  =====================
    ``ambertools`` subdirective              Type            Description (default)
    =====================================    ==============  =====================
    ``charge_method``                        string          charge model used by ``antechamber`` (default ``gas``)
    ``net_charge``                           integer         formal charge of each molecule, passed to ``antechamber -nc`` (default ``0``)
    ``atom_type``                            string          atom-type set used by ``antechamber``, ``parmchk2`` and ``tleap`` (default ``gaff``)
    =====================================    ==============  =====================

    For now, you can choose any charging method compatible with ``antechamber``.  The ``antechamber`` directive is optional.

    These three directives determine the parameters ``antechamber`` produces, so each parameterization records them in a ``.parm`` file checked into the library alongside its ``gro``, ``top``, ``itp``, ``tpx`` and ``grx`` files.  Because the library is keyed on molecule name, that record is the only thing distinguishing a molecule parameterized one way from the same molecule parameterized another.  A run whose directives disagree with the record re-parameterizes rather than reusing the cached files, and says which directives differed.

    A parameterization checked in before htpolynet 2.3 has no such record.  It is still used **as-is**, which means the build takes whatever charges that entry holds even if they were produced with a different ``charge_method``; the run warns, per molecule and again in a block at the end of the parameterization stage, that it could not verify what produced them.  This is deliberate -- it keeps an existing library working instead of re-parameterizing it wholesale -- but it does mean that asking for ``bcc`` against a pre-2.3 library can still yield ``gas`` numbers.  Rebuild those entries with ``--force-parameterization`` if a build's charges need to be certain.  Note that re-parameterizing after a mismatch does **not** replace the library entry unless you also pass ``--force-checkin``, so a library holding a ``gas`` parameterization keeps it even after a ``bcc`` run.

* ``GAFF``

    In very rare instances, AmberTools will generate GAFF atom types and parameters that are internally inconsistent, or at least are not understandable by the ``parmed`` package that translates them into Gromacs topology files.  Directives in this section instruct ``htpolynet`` how to resolve these inconsistencies.  

    =====================================    ==================  =====================
    ``GAFF`` subdirective                    Type                Description (default)
    =====================================    ==================  =====================
    ``resolve_type_discrepancies``           **resolve record**  Directives for resolving atom type discrepancies
    =====================================    ==================  =====================


    ==============================================    =========================================
    ``GAFF.resolve_type_discrepancies`` directives    Description (default)
    ==============================================    =========================================
    ``typename``                                      Name of interaction type in ``top`` file (``dihedraltypes``)
    ``funcidx``                                       Function index assigned to this type in the ``top`` file
    ``rule``                                          What rule to apply to resolve discrepancies (``stiffest``)
    ==============================================    =========================================

    Currently, the only type of discrepancy that can be handled is one in which four atom types define one dihedral type in one ``top`` file but a different one in another.  If these to topologies are merged to make a composite system, Gromacs will flag this as an error in the topology file.  This directive allows you to check for all such discrepancies before Gromacs does, and keep only the one that follows the stipulated ``rule``.  For dihedrals, the ``stiffest`` rule means that of the conflicting types, the one with the largest energetic parameters is kept.

* ``constituents``
  
    This directive is a set of one or more "key":"record" pairs in which each key is the name of a molecule (here, "STY") and the record is a dictionary of keyword:value pairs.  The allowable keywords in a ``constituent`` record are as follows.

    =====================================    =================  =====================
    ``constituents`` record keyword          Type               Description (default)
    =====================================    =================  =====================
    ``count``                                int                (optional) number of these molecules in the system
    ``stereocenters``                        list               (optional) list of names of chiral carbon atoms if any
    ``symmetry_equivalent_atoms``            list               (optional) list of sets of symmetry equivalent atom names, if any
    ``conformers``                           conformers record  (optional) parameters specifying if and how initial conformers are generated
    ``smiles``                               string             (optional) SMILES for in-config monomer generation; htpolynet runs ``obabel`` (or RDKit) to materialize ``lib/molecules/inputs/<NAME>.mol2`` before parameterization
    ``rename_atoms``                         dict               (optional) ``{1-based-mol2-index: atom-name}`` map applied after obabel writes the mol2 (obabel path; pairs with ``smiles``)
    ``reactive_atoms``                       dict               (optional) ``{smiles-atom-map-label: atom-name}`` map (RDKit path; pairs with a ``smiles`` that uses ``[*:N]`` atom-mapping tokens).  ``rdkit`` is a required runtime dependency and is pulled in by any normal ``htpolynet`` install
    ``frcmod``                               string             (optional) path, relative to the configuration file, to an Amber ``frcmod`` of parameter overrides for this molecule; see below
    =====================================    =================  =====================

    ``frcmod`` corrects GAFF for one monomer without touching the others, for
    example a torsion refit against quantum-chemical scans.  ``tleap`` loads
    the file after GAFF and ``parmchk2``, so its entries win.  It is loaded
    for the molecule itself and for every template that contains it: cure
    products, capping and repair templates, and the multiply-substituted
    templates htpolynet derives.  The file's contents are part of each of
    those parameterizations' cache records, so editing it re-parameterizes
    them on the next run rather than reusing the old numbers.

    .. code-block:: yaml

       constituents:
         BAF:
           smiles: "FC(F)(F)C(c1ccc([OH:1])cc1)(c1ccc([OH:2])cc1)C(F)(F)F"
           reactive_atoms: {1: O1, 2: O2}
           frcmod: baf-torsion.frcmod

    One restriction follows from how GROMACS topologies work.  A system holds
    one parameter table per bonded *type* (four atom types, for a dihedral),
    shared by every molecule in it.  So an override cannot apply to a type in
    one molecule and leave the same type alone in another: if a molecule
    built without the frcmod also contains the overridden type, with
    different parameters, setup stops and names it.  Either load the same
    frcmod for every constituent that contains the type, or give the atoms
    you want to treat differently their own atom type.

    In the example below, we are requesting a system of 100 styrene molecules.  The key ``STY`` signals to ``htpolynet`` that it should look for either ``STY.mol2`` or ``STY.pdb`` in ``./lib/molecules/inputs`` **or** it should look for ``STY.gro``, ``STY.itp``, ``STY.top``, and ``STY.grx`` in ``./lib/molecules/parameterized``.  The latter is the case if either ``htpolynet run`` or ``htpolynet parameterized`` has already been run with ``STY.mol2`` or ``STY.pdb``.  If neither set of files is found but the constituent carries a ``smiles`` key, ``htpolynet`` generates the input ``mol2`` itself.  Multiple records in ``constituents`` should all have the "key":"record" syntax and be separated by commas.

    Any entry in ``constituents`` for which the ``count`` is 0 or missing is treated as an intermediate for which stereocenters or symmetry-equivalent atoms are specified.

    ``htpolynet`` allows you the option to use multiple conformers of flexible molecules to build the initial liquid system.  It can use either ``obabel``'s ``confomers`` capability or an MD simulation via ``gromacs`` to generate these.  The ``conformers`` record has two subdirectives:

    =====================================    ==========================  =====================
    ``conformers`` record keyword            Type                        Description (default)
    =====================================    ==========================  =====================
    ``count``                                int                         (required) number of unique conformers to generate (per stereoisomer)
    ``generator``                            conformer generator record  (optional) parameters specifying how conformers are generated
    =====================================    ==========================  =====================

    The ``conformers.generator`` record has several subdirectives:

    =====================================    ===========================  =====================
    ``generator`` record keyword             Type                         Description (default)
    =====================================    ===========================  =====================
    ``name``                                 str                          (required) ``obabel`` or ``gromacs``
    ``params``                               generator parameters record  (optional) parameters specifying the generator's operation (only relevant for ``gromacs``)
    =====================================    ===========================  =====================

    The ``conformers.generator.params`` record has several subdirectives:

    =====================================    ===========================  =====================
    ``params`` record keyword                Type                         Description (default)
    =====================================    ===========================  =====================
    ``ensemble``                             str                          ``nvt`` is the only option that makes sense
    ``temperature``                          float                        (optional) Temperature of the conformer-generating MD simulation
    ``ps``                                   float                        (optional) Duration of the conformer-generating MD simulation
    ``pad``                                  float                        (optional) Box-size padding for the vacuum MD simulation
    =====================================    ===========================  =====================

    If an entry in ``constituents`` has no ``confomers`` member directive, then the confomer used for it in building the system is whatever one is in either the input ``mol2`` file (after AmberTools and ``parmed`` convert it to a ``gro`` file) **or** the ``gro`` file of the constructed molecule.

* ``densification``

    This directive instructs ``htpolynet`` how to run the initial densification of the fresh simulation system.  It has two subdirectives:

    =====================================    ==============  =====================
    ``densification`` subdirective           Type            Description (default)
    =====================================    ==============  =====================
    ``initial_density``                      float           density in kg/m^3 at which molecules are placed randomly into a box to make the initial coordinates (default 200.0)
    ``equilibration``                        list            list of **MD records** 
    =====================================    ==============  =====================

    The ``equilibration`` subdirective should contain one or more *MD records*. An MD record is a dictionary of keyword:value pairs:

    =====================================    ==============  =====================
    MD record keyword                        Type            Description
    =====================================    ==============  =====================
    ``ensemble``                             string          (required) min (minimization), npt, or nvt
    ``temperature``                          float           (required if ``ensemble`` is nvt or npt) Temperature in K assigned to ``ref_t`` in Gromacs ``mdp`` file
    ``pressure``                             float           (required if ``ensemble`` is npt) Pressure in bar assigned to ``ref_p`` in Gromacs ``mdp`` file
    ``nsteps``                               int             (optional; required if ``ps`` not provided) Duration of MD simulation in number of time steps
    ``ps``                                   float           (optional; required if ``nsteps`` not set) Duration of MD simulation in picoseconds
    ``repeat``                               int             (optional) number of times to repeat this simulation in series; default is 0 (i.e., run once)
    =====================================    ==============  =====================

    The ``repeat`` subdirective is especially useful for densifications that start at very low initial densities.  It is better to run several short NPT simulations than a single long one so that the box size shrinkage doesn't overwhelm Gromacs' domain decomposition algorithm.

    An NPT record may also carry a ``converge`` subdirective, which repeats the stage until the density settles instead of running for however many steps you guessed.  It is off unless configured, so existing configurations reproduce exactly.

    =====================================    ==============  =====================
    ``converge`` subdirective                Type            Description (default)
    =====================================    ==============  =====================
    ``tolerance``                            float           largest acceptable standard error of the mean density, in kg/m^3 (default 1.0)
    ``max_repeats``                          int             ceiling on extension segments; reaching it is reported as a *failure* to converge, not as success (default 10)
    ``min_samples``                          int             refuse to judge a window with fewer energy frames than this (default 50)
    ``drift_sems``                           float           how far the window's two halves may differ, in units of that difference's own standard error (default 2.0)
    =====================================    ==============  =====================

    Two things are worth understanding before setting a ``tolerance``.  The criterion is an **autocorrelation-corrected** standard error, :math:`\sigma/\sqrt{N/\tau_\mathrm{int}}`: an NPT cell density is correlated over hundreds of steps, so successive energy frames are not independent samples and a naive standard error is optimistic by roughly a factor of 1.5 -- a gate built on it stops early.  The correction also makes the tolerance mean the same thing for a small box as for a large one, since :math:`\sigma/\mu` falls as :math:`1/\sqrt{N_\mathrm{atoms}}` while :math:`\tau_\mathrm{int}` is roughly size-independent.  Second, a small standard error alone is not taken as convergence, because a trace that is still climbing can look tight in every short window; the two halves of the window are compared as well, and the stage is extended if they differ by more than ``drift_sems`` times that difference's own standard error.

    Each extension is a fresh ``mdrun`` and regenerates velocities, exactly as ``repeat`` does, so the criterion is applied to the last segment rather than to a concatenation of segments that do not share a velocity history.

    .. note::

       This is deliberately offered for **densification** and not for the CURE relaxation stages.  Densification is one-shot and involves a large volume change, which is what the machinery is for.  Gating the cure loop was measured and rejected: past the gel point the unreacted species are bonded into the network and are topologically constrained rather than merely slow, so the effective diffusion exponent falls to about 0.14 and restoring the mobility a gate would wait for takes of order :math:`10^4` times the relaxation time.  A gate there would correctly report "not settled" on every late iteration with nothing to be done about it.

    .. warning::

       The gate reports convergence of the **Berendsen** barostat that ``npt.mdp`` currently uses, which does not sample a correct NPT ensemble.  A converged number here means the box has stopped changing under that barostat, not that it is the correct equilibrium density for the force field.  Treat the tolerance as a reproducibility criterion rather than a physical one until the barostat is changed.

* ``precure``
    
    The ``precure`` directive instructs ``htpolynet`` on running a series of MD simulations after densification but before the cure.  There are three allowable subdirectives for ``precure``: 

    =====================================    =================    =====================
    ``precure`` subdirective                 Type                 Description (default)
    =====================================    =================    =====================
    ``preequilibration``                     MD record            optional MD simulation
    ``anneal``                               **Anneal record**    Description of an annealing simulation after the optional ``preequilibration``
    ``postequilibration``                    MD record            optional MD simulation         
    =====================================    =================    =====================

    Both the ``preequilibration`` and ``postequilibration`` directives contain MD records described above.  The *Anneal record* has the following subdirectives:

    =====================================    =================    =====================
    Anneal record subdirective               Type                 Description (default)
    =====================================    =================    =====================
    ``ncycles``                              int                  number of annealing cycles
    ``initial_temperature``                  float                (optional) Initial temperaure in K, really only sets the ``gen-temp`` ``mdp`` parameter 
    ``cycle_segments``                       list                 list of **cycle records**
    =====================================    =================    =====================

    A **cycle record** corresponds to an "annealing-point" in the Gromacs ``mdp`` file.  

    =====================================    =================    =====================
    Cycle record subdirective                Type                 Description (default)
    =====================================    =================    =====================
    ``T``                                    float                Targe temperature in K 
    ``ps``                                   float                cycle duration; if prior ``T`` is different, simulation is *brought to* this ``T`` in this amount of time; if prior ``T`` is the same, simulation is *held at* this ``T`` for this amount of time.
    =====================================    =================    =====================

    Each cycle consists of one pass through the cycle segments.  In the example here, one cycle consists of Gromacs taking the system from 300 to 600 K in the first 20 ps, then holding at 600 for 20 pm, then reducing to 300 K over 20 ps and holding it there for 20 ps.

* ``CURE``
   
    This directive contains all instructions governing the :ref:`CURE algorithm <cure_section>`.  There are five possible subdirectives:

    =====================================    =================    =====================
    ``CURE`` subdirective                    Type                 Description (default)
    =====================================    =================    =====================
    ``controls``                             list                 Control parameter values
    ``drag``                                 list                 Dragging parameter values
    ``relax``                                list                 Bond relaxation parameter values
    ``equilibrate``                          MD record            CURE iteration equilibration parameters
    ``gromacs``                              list                 any ``mdp`` keyword:value pairs to include in all ``mdp`` files in the ``CURE`` sequence
    =====================================    =================    =====================

.. _cure.controls:

    * ``CURE.controls`` parameters

        ==================================    =================   ======================
        ``CURE.controls`` parameter           Type                Description (default)
        ==================================    =================   ======================
        ``search_radius``                     float               initial search radius in nm (default 0.5)
        ``radial_increment``                  float               increment in nm by which search radius is increased if not enough bonds are found at current radius (default 0.05)
        ``min_bonds_per_iteration``           int                 grow the radius until at least this many bonds have been identified, then proceed with whatever the count is at the current radius (default 10). Effective floor is clamped against the remaining bonds needed to reach ``desired_conversion`` and against ``max_conversion_per_iteration``, so demanding e.g. 50 near end-of-cure never stalls the build. Setting to 1 reproduces the un-batched legacy behavior.
        ``max_iterations``                    int                 absolute maximum number of allowed iterations (default 100)
        ``desired_conversion``                float [0-1]         target **bond** conversion -- bonds formed over bonds possible -- between 0 and 1.0 (default 0.5). This is not the crosslinker conversion an experiment measures, which is lower; see :ref:`what repair reports <postcure_repair_reporting>`.  Note the default is a half-cure, not the 0.95 the examples use; every shipped example that cures sets the value explicitly rather than relying on it.
        ``late_threshold``                    float [0-1]         conversion above which bond probabilities are ignored (default 0.85)
        ``max_conversion_per_iteration``      float [0-1]         upper limit, as a fraction of total reactable bonds, on the new bonds formed in any single iteration (default 1.0)
        ``min_allowable_bondcycle_length``    int                 minimum number of C atoms allowed in a cycle of C-C bonds that form via polymerization (default -1; any value <= 0 disallows all such cycles)
        ``completion_bias``                   bool                rank bond candidates by how many bonds their ``B``-side residue already carries, and only then by distance, so that partly-reacted crosslinkers are completed before untouched ones are started (default ``False``)
        ==================================    =================   ======================

      The ``min_bonds_per_iteration`` knob batches small late-stage finds together: instead of accepting whichever 1-2 bonds happen to lie within the initial 0.5 nm radius and immediately moving on to relax + equilibrate, the search grows the radius until at least 10 (the default) eligible bonds have been gathered.  On the DGEBA/PACM example, raising it from 1 to 10 cuts the cure iteration count from 41 to 15; raising further to 20 only buys one more iteration.  The default of 10 is roughly the diminishing-returns sweet spot.

      ``completion_bias`` changes which bonds a CURE iteration chooses, and with it what ``desired_conversion`` means physically, so it is off by default and every shipped example runs without it.  Candidates are ranked by separation alone by default.  Separation is uncorrelated with how many bonds a crosslinker already carries, so in an A2+B3 system the bonds spread evenly across every ``B`` in the box rather than finishing any of them: at a bond conversion of 0.90, a trifunctional crosslinker placed at random is complete with probability 0.90\ :sup:`3`, about 0.73.  Turning ``completion_bias`` on makes the number of bonds already on the candidate's ``B``-side residue the primary sort key, with separation as the tie-break inside each group, so a crosslinker with two of three sites filled is finished before an untouched one is started.  The ``B`` side is where htpolynet's A2+B3 idiom puts the crosslinker; declare it as ``A`` instead and the bias completes your *bridges*, so the first cure iteration compares the two sides' functionality and warns if the more functional one is on ``A``.  Nothing else about the search changes -- same radius growth, same dragging and relaxation, same probabilities, same cycle handling.

      Which conversion this matters for depends on the chemistry.  Where a partly-reacted crosslinker is a real, stable species, the unbiased ranking is the more faithful one.  Where it is a reactive intermediate -- a cyclotrimerizing cyanate ester, whose triazine ring either closes or does not -- the biased ranking is, and the difference shows up directly in what :ref:`postcure repair <postcure_repair>` has to dismantle.  Note that the two settings reach different structures at the same ``desired_conversion``: a config that yields a crosslinker conversion of 0.76 at ``desired_conversion: 0.90`` unbiased will yield about 0.90 biased.  Comparisons across the setting have to be made at matched crosslinker conversion, not matched ``desired_conversion``.

      The ``min_allowable_bondcycle_length`` refers to the fact that in systems that polymerize via activation of carbon-carbon double bonds, it is possible in the htpolynet implementation that the "head" of a chain of C-C bonds can attack the "tail" and form a cycle, because those represent atom types that can react.  It is unclear whether such cycles actually form; if a monomer remains bound to a radical initiator it is hard to see how the head of the growing chain could attack it, but maybe it could.  Setting ``min_allowable_bondcycle_length`` to zero or less (the default is -1) disallows any bonds that would form cycles involving only atoms that were once part of C=C double bonds.  (Think about the backbone of polystyrene, for example.)  In a given CURE iteration, htpolynet tests the full set of suggested bonds to see if together they result in any cycles, and for each nascent cycle longer than ``min_allowable_bondcycle_length``, htpolynet will disallow the nascent bond that has the longest initial length.

.. _cure.drag:

    * ``CURE.drag`` parameters:  Dragging refers to a series of MD simulations (called "stages") in which harmonic restraints are applied to each pair of atoms assigned to form a bond, but **before** the bonds actually form.  Dragging is useful to reduce 1-4 distances that ultimately arise when bonds form.  Each stage in the series uses a specially modified topology file in which "new" bonds of type 6 are added, one for each pair of to-be-bonded atoms. Each of these bonds has a parameter ``kb``, the spring constant, and ``b0``, the equilibrium length.  The ``drag`` directive governs how those ``b0`` parameters are linearly decreased through the set of stages to slowly bring the atoms closer together.   The ``limit`` parameter is the target distance of dragging, and ``increment`` determines the number of stages it will take to get there.

        =================================    =================   ======================
        ``CURE.drag`` parameter              Type                Description (default)
        =================================    =================   ======================
        ``increment``                        float               minimum amount by which target ``drag`` distance is decreased in steps (default 0.08)
        ``limit``                            float               distance in nm to which all bonds are dragged (default 0.3)
        ``equilibration``                    MD record           describes the MD simulations used to equilibrate at each stage 
        =================================    =================   ======================

.. _cure.relax:

    * ``CURE.relax`` parameters:  Relaxation refers to a series of MD simulations (also called "stages") in which the ``kb`` and ``b0`` parameters of each new bond are "attenuated" from a weak (low ``kb``), long (large ``b0``) state to the state dictated by the force field.  The ``increment`` determines the number of stages are performed.

        =================================    =================   ======================
        ``CURE.relax`` parameter             Type                Description (default)
        =================================    =================   ======================
        ``increment``                        float               minimum amount by which ``b0`` parameters are decreased in steps (default 0.08)
        ``equilibration``                    MD record           describes the MD simulations used to equilibrate at each stage 
        =================================    =================   ======================

    * ``gromacs`` parameters:  These parameters govern modification to ``mdp`` files used in the dragging and relaxation MD simulations.  ``htpolynet`` adjusts the cutoff distances to conform to the longest unrelaxed bond in the system, and the ``rdefault`` parameter provides the floor below which it will not go any lower.

        =================================    =================   ======================
        ``CURE.gromacs`` parameter           Type                Description (default)
        =================================    =================   ======================
        ``rdefault``                         float               minimum cutoff radius (default 0.9)
        =================================    =================   ======================


* ``postcure`` 

    The ``postcure`` directive instructs ``htpolynet`` on running a series of MD simulations after cure.  Its form is identical to that of ``precure``, namely with optional ``preequilibration``, ``anneal``, and ``postequilibration`` subdirectives.

    .. warning::

       **The density at the end of ``postequilibration`` is not an equilibrated density, and should not be quoted as one.**  A crosslinked network only relaxes appreciably above its glass transition, and a postcure anneal whose peak sits near *T*:sub:`g` spends very little time where the network can move -- so the postcure plateau reports an under-relaxed structure no matter how long the plateau itself runs.

       Measured on four independent BPA cyanate-ester builds: the postcure NPT plateau gave 1.1712 ± 0.0029 g/cm³ at 300 K, while the same systems melted and slowly re-cooled gave 1.1983 ± 0.0037 g/cm³ at 293 K -- a 2.31% difference between two extractions of the same quantity.  Against experiment (1.205 by mercury dilatometry, 1.202 by flotation, same monomer), the plateau value is 2.6-2.8% low and the re-cooled value 0.3-0.6% low.  A build reporting the plateau density would have read as a force-field failure when it was a protocol artifact.

       If you need a density, extract it from a melt-and-slow-cool trajectory (an :doc:`htpolynet postsim <configs-for-postsim>` ladder does this) rather than from the postcure plateau.  The same caution applies to any property that depends on the structure having relaxed.

.. _reactions:

* ``reactions``

    The ``reactions`` directive contains a list of **reaction records** that specify the chemisty of any bonds that form to either build molecular constituents or polymers/crosslinks. A reaction is defined by the precise pairs of atom types that become new covalent bonds.  To precisely define each such pair, the reaction must also identify one or more reactant molecules.  Each reaction also names a single product molecule.  htpolynet will build oligomer templates using these reactions and then GAFF-parameterize them.  The parameterizations are used during CURE to re-type atoms and reset charges after each new bond is formed.

    ==============================  ==========  =================
    ``reaction`` record directives  Type        Description
    ==============================  ==========  =================
    ``name``                        str         descriptive name
    ``stage``                       str         one of ``cure``, ``cap``, ``build``, or ``param``
    ``probability``                 float       probability that bond will form in one iteration if identified (1.0)
    ``reactants``                   dict        keyword: reactant key, value: reactant molecule name
    ``product``                     str         name of product molecule
    ``atoms``                       dict        keyword: atom key, value: **atom record**
    ``bonds``                       list        list of **bond records**, one item per bond formed in reaction
    ==============================  ==========  =================

    The ``stage`` value signifies how ``htpolynet`` uses the reaction.  It will generate GAFF parameters and topologies for any product of a reaction with stage ``cure``, ``cap``, or ``param``.  ``cure`` reactions are those assigned to take place during CURE.  ``cap`` reactions are optional and take place once the CURE has finished; these can be used to revert the active form of any unreacted monomers back to their proper forms.  ``param`` reactions are only performed in the beginning when molecular constituents are being built.  If you want to build the molecular constituents out of simpler monomers, you will likely want to use ``param`` reactions.  
    
    .. 
        Finally, if you have constituents that are themselves made of repeating monomeric components, you need parameterize on one such reaction, and others can be specified to be state ``build``, for which no parameters are generated, only a bond is formed.

    The ``atoms`` directive is a dictionary of atom records where the key is an atom "key", which is referenced in bond record.

    * Atom records uniquely identify atoms in reactants, assigning them a shorthand key that is used in subsequent bond records.
        
        ======================== ============== =================
        Atom record subdirective type           Description
        ======================== ============== =================
        ``reactant``             arb.           Reactant key that references the ``reactants`` directive of the reaction
        ``resid``                int            Residue number in the reactant containing this atom
        ``atom``                 str            Atom name (originates in monomer ``mol2`` or ``pdb`` file)
        ``z``                    int            Number of possible bonds atom can participate in
        ======================== ============== =================

    * Bond records specify the bond(s) that form during this reaction.

        ======================== ============== =================
        Bond record subdirective type           Description
        ======================== ============== =================
        ``atoms``                list           The two atom keys that define the atoms that form the bond
        ``order``                int            Order (1=single, 2=double) of resulting bond
        ======================== ============== =================

    In the example here, we define two unique reactions.  One is the C1-C2 bond that links two styrene monomers, and the other is the *intramolecular* C1-C2 double bond that "reverts" the active form of a monomer back to its "proper" form.  Since that reaction's ``stage`` is ``cap``, this signifies that it is formed only **after** CURE has finished.

.. _configuration_run_example:

A Simple Configuration Example:  Polymerizing styrene
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Consider the following methodological description:

A system of 100 molecules of styrene was prepared by initially placing 100 styrene molecules in a box at a density of 0.3 g/cc.  This system was minimized and subjected to 10 ps of NVT MD at 300 K and then 200 ps of NPT MD at 300 K and 10 bar.  The pressure was reduced to 1 bar and NPT MD was run for an additional 200 ps.  The system was then put through two cycles of constant-pressure annealing at 1 bar in which a cycle is 20 ps at 300 K, 20 ps of linear temperature increase to 600 K, 20 ps at 600 K, and 20 ps of linear temperature decrease back to 300 K.  This was followed by 100 ps of equilibration at 1 bar and 300 K.  Then the CURE algorithm was used with an initial search radius of 0.5 nm and increment of 0.25 nm to create new intermonomer bonds until a conversion of 95% was reached.  At each CURE iteration, the system was driven using a series of steered MD simulations to first force all pairs of atoms designated to be bonded but are more than 0.6 nm apart to approach one another to within 0.3 nm. When all such designated pairs were within 0.3 nm of each other, any sacrificial H's were deleted and new bonds and accompanying angles and dihedrals were introduced into the topology.  The bond lengths were then linearly reduced to their equilibrium values through a series of MD simulations.  Once all bonds were at their equilibrium lengths, the iteration was completed by a 100 ps NPT MD simulation at 1 bar and 300 K, before beginning the next CURE iteration.  Following CURE, the system was again annealed following the previous protocol and then NPT equilibrated at 1 bar and 300 K for 200 ps.

A simple YAML ``htpolynet run`` configuration file that encodes this description might look like this.

.. code-block:: YAML

  Title: polystyrene
  constituents:
    STY:
      count: 100
  densification:
    initial_density: 300.0  # kg/m3
    equilibration:
      - ensemble: min
      - ensemble: nvt
        temperature: 300
        ps: 10
      - ensemble: npt
        temperature: 300 
        pressure: 10
        ps: 200
  precure:
    preequilibration:
      ensemble: npt
      temperature: 300        # K
      pressure: 1             # bar
      ps: 200
    anneal:
      ncycles: 2
      initial_temperature: 300
      cycle_segments:
      - T:  300
        ps: 0
      - T:  600
        ps: 20
      - T:  600
        ps: 20
      - T:  300
        ps: 20
      - T:  300
        ps: 20
    postequilibration:
      ensemble: npt
      temperature: 300        # K
      pressure: 1             # bar
      ps: 100
  CURE:
    controls:
      search_radius: 0.5         # nm
      radial_increment: 0.25     # nm
      min_bonds_per_iteration: 10
      max_iterations: 150
      desired_conversion: 0.95
      late_threshold: 0.85
    drag:
      trigger_distance: 0.6   # nm
      increment: 0.08         # nm
      limit: 0.3              # nm
      equilibration:
        - ensemble: min
        - ensemble: nvt
          temperature: 600
          nsteps: 1000
        - ensemble: npt
          temperature: 600
          pressure: 1
          nsteps: 2000
    relax:
      increment: 0.08         # nm
      equilibration:
        - ensemble: min
        - ensemble: nvt
          temperature: 600
          nsteps: 1000
        - ensemble: npt
          temperature: 600
          pressure: 1
          nsteps: 2000
    equilibrate:
      ensemble: npt
      temperature: 300       # K
      pressure: 1            # bar
      ps: 100
      gromacs:
      rdefault: 0.9 # nm
  postcure:
    anneal:
      ncycles: 2
      initial_temperature: 300
      cycle_segments:
        - T:  300
          ps: 0
        - T:  600
          ps: 20
        - T:  600
          ps: 20
        - T:  300
          ps: 20
        - T:  300
          ps: 20
    postequilibration:
      ensemble: npt
      temperature: 300       # K
      pressure: 1            # bar
      ps:  100
  reactions:
    - name:        sty1_1
      stage:       cure
      reactants:
        1: STY
        2: STY
      product:     STY~C1-C2~STY,
      probability: 1.0
      atoms:
        A:
          reactant: 1
          resid: 1
          atom: C1
          z: 1
        B:
          reactant: 2
          resid: 1
          atom: C2
          z: 1
      bonds:
        - atoms: [A, B]
          order: 1
    - name:         styCC,
      stage:        cap,
      reactants:
        1: STY
      product:      STYCC,
      probability:  1.0,
      atoms:
        A:
          reactant: 1
          resid: 1
          atom: C1
          z: 1
        B:
          reactant: 1
          resid: 1
          atom: C2
          z: 1
      bonds:
        - atoms: [A, B]
          order: 2

A note on YAML format: notice here that the left-most indents are all words; each of these is a key that refers to a specification for the run.  Here is what this configuration specifies.  First, we are starting with 100 styrene molecules.  So, ``htpolynet`` expects to find an **input** file name ``STY.mol2`` **or** ``STY.pdb`` in the system or user ``molecules/inputs`` library.  A liquid system of these 100 styrenes is densified starting from an initial density of 300 kg/m\ :sup:`3`\ 
using first an energy minimzation, then an NVT MD simulation at 300 K for 10 picoseconds, then an NPT MD simulation at 300 K and 10 bar for 200 ps.  (The 10 bar helps to ensure rapid densification, but don't feel pressured to use it.)  To prepare for cure, the system is brought to 1 bar in the precure stage, where it is also annealed by raising the temprature from 300 to 600 K and bringing it back down, for two cycles.  The CURE is run such that 

1. The search radius begins at 0.5 nm and goes up in increments of 0.25 nm;
2. The algorithm bails out after 150 iterations; and
3. The desired cure conversion is 95%;

Prebond dragging is permitted if any newly identified bond is more than 0.6 nm in length, and the dragging happens in increments of 0.08 nm and each increment involves an energy minimization, an NVT MD simulation, and an NPT MD simulation, all at 600 K. (I find curing at elevated temperature keeps the system from jamming up, but don't feel forced to use this temperature.)  Bond relaxation takes place using a similar series of MD stages.  Remember that dragging is performed on the system **before** bonds are formed and atoms deleted, while bond relaxation occurs **after** the bonds are formed and the sacrificial, valence-conserving H atoms are deleted.  Finally, when all new bonds are relaxed, a single NPT MD simulation is performed to end an iteration.  Postcure involves an annealing simulation much like the precure stage, followed by an NPT MD simulation.

Finally, we stipulate the reactions.  In this system, there is really only one reaction: the one in which the C1 of one styrene bonds to the C2 of another.  The reaction named ``sty1_1`` specifies this reaction, and causes ``htpolynet`` to parameterize the dimeric product named ``STY~C1-C2~STY``.  This molecule provides a template for atom types, charges, and new bonded interactions that must be merged into a system if such a bond forms.  The other reaction, ``styCC``, specifies a ``cap`` reaction that reverts any unreacted styrene back to its proper form (with the C-C double bond). Capping reactions are 100\% optional; don't feel forced to use them.
