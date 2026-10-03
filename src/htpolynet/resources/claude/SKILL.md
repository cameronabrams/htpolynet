---
name: htpolynet
description: Build crosslinked polymer systems with htpolynet — authoring a run configuration, parameterizing monomers, running a cure, and reading the results. Use when the task involves an htpolynet YAML config, the `htpolynet` command (run, input-check, gen-slurm-script, plots, postsim, analyze, make-viz), a monomer's active form or sacrificial hydrogens, a cure that stalls below its target conversion, or sizing an htpolynet build for a cluster.
---

# Working with htpolynet

htpolynet builds crosslinked polymer systems: it parameterizes monomers with
GAFF, packs them, and then *cures* the system by repeatedly forming bonds
between reactive atoms until it reaches a target conversion. The output is a
Gromacs topology and coordinate pair plus the analysis to characterize it.

## The working loop

1. Fetch the bundled example closest in **reaction topology** and run it
   unmodified once.
2. Edit it into your system. Describe monomers in their **active form**.
3. `htpolynet input-check <cfg>` — fast, touches nothing, reports the atom
   count you size a core request against.
4. `htpolynet run <cfg>`, in the background, or via a batch script.
5. `htpolynet plots build --proj proj-0` to see what you got.

## Start from an example, never from an empty file

A working configuration is a few hundred lines and most of it is not
chemistry-specific. Writing one from scratch is the slowest route and the one
that produces silent errors.

```bash
htpolynet fetch-example 1     # linear chain growth (polystyrene)
htpolynet fetch-example 2     # A2 + B4 thermoset (bisGMA/styrene)
htpolynet fetch-example 3     # amine + epoxy (PACM/DGEBA)
htpolynet fetch-example 6     # cyclotrimerization (cyanate ester)
```

Pick by **reaction topology**, not by chemistry. Building a step-growth
thermoset from a difunctional and a tetrafunctional monomer? Example 2 beats
an example that shares a functional group but polymerizes by chain growth.

Run it unmodified first. That confirms the toolchain and gives you a
known-good log to diff against when your own config misbehaves.

## Monomers are described in their active form

The most common conceptual error, and it produces a build that completes and
is wrong rather than one that fails.

htpolynet conserves valence: when two atoms bond, each gives up a
**sacrificial hydrogen**. So you do not describe the textbook monomer — you
describe the monomer with its reactive sites already saturated. Styrene is
described as ethylbenzene.

Prefer SMILES written directly into `constituents` with RDKit atom-map
labels, which name reactive atoms by chemical identity rather than by file
index:

```yaml
constituents:
  STY:
    smiles: "c1ccccc1[CH2:1][CH3:2]"
    reactive_atoms: {1: C1, 2: C2}
```

In SMILES bracket atoms the hydrogen count is explicit: `[C:1]` means *zero*
implicit hydrogens. For an sp³ carbon you almost always want `[CH2:1]` or
`[CH3:1]`. A wrong count usually surfaces much later as a missing GAFF angle
parameter in `tleap`, a long way from its cause.

## Check before spending compute

`htpolynet input-check <cfg>` reports the initial atom count in about a
second. Run it every time you change monomer counts.

Know what it does **not** check. It does not verify that
`symmetry_equivalent_atoms` groups really are topologically equivalent, that
A2 + B3 site counts balance, or that `desired_conversion` is reachable given
the monomer counts. Those are silent failure modes: a wrong symmetry group
generates reaction templates the cure stage never matches, and the build
spends its full wall-clock reaching a conversion it could never have reached.

## Running

Builds are slow — production systems are hours to a day. Run them in the
background, or generate a batch script rather than writing one:

```bash
htpolynet gen-slurm-script <cfg>
```

Size the core request from the `input-check` atom count.

**If you are using the container image, target CPU partitions.** Its Gromacs
comes from conda-forge, built against OpenCL, which Gromacs no longer uses to
drive NVIDIA devices. The image cannot use a GPU: `--gres=gpu` and `--nv` buy
nothing and may cost queue time.

Parameterization runs first and is cached across projects and directories.
The first run of a new chemistry pays several minutes per species; later runs
reuse it. Editing a monomer's structure **without renaming it** silently
reuses the old parameterization — `--force-parameterization` rebuilds.

## Builds are not reproducible run to run

There is no seed control. Conformer selection, the per-bond probability test,
and Gromacs velocity generation are all unseeded, so two runs of one config
diverge. Convenient in one direction — independent replicas of a
quenched-disorder ensemble come free, and each `htpolynet run` in the same
base directory makes its own project directory — but an exact rebuild is not
available. Record the commit (`htpolynet info`) when the result matters.

## Reading the results

```bash
htpolynet plots build --proj proj-0     # traces, reaction graph, cluster sizes
htpolynet make-viz -top final.top -gro final.gro
htpolynet postsim postsim.yaml --proj proj-0
htpolynet plots post --proj proj-0      # E and Tg
htpolynet analyze analyze.yaml --proj proj-0
```

**The postcure NPT plateau is not an equilibrated density.** It is an
under-relaxed structure, and on the cyanate-ester example as it shipped before
2.15.0 (the pre-formed-triazine stand-in) it sat about 2.3% below what the same system gives after a melt and slow re-cool.
Do not report a plateau density as a force-field result.

## Subcommand routing

| Intent | Command |
|---|---|
| Build a system | `htpolynet run <cfg>` |
| Size it before queueing | `htpolynet input-check <cfg>` |
| Stop after parameterization | `htpolynet run --param-only <cfg>` |
| Submit to a cluster | `htpolynet gen-slurm-script <cfg>` |
| Plot a finished build | `htpolynet plots build --proj <dir>` |
| Post-build MD / analysis | `htpolynet postsim`, `htpolynet analyze` |
| Rebuild VMD viz files | `htpolynet make-viz` |
| Report the environment | `htpolynet info` |

## When something goes wrong

- **The cure stalls below the target conversion.** Usually a reaction
  template that never matches: suspect `symmetry_equivalent_atoms` and
  reactive-atom names before suspecting the chemistry.
- **A `tleap` error about a missing parameter.** Usually an atom typed
  wrongly from a hydrogen-count error in a SMILES bracket atom.
- **Results that do not match the config you are reading.** Suspect a stale
  cache entry: a structure edit without a rename.
- **A build that dies during cure.** `-restart` is experimental and currently
  broken at the cure stage, which is the worst place to lose work. Treat a
  dead cure as a restart from the beginning.
- **Anything else.** `htpolynet info` reports the versions of every external
  tool the build depends on; include its output in any bug report.

## Do not

Do not write workarounds into this file. If the tool surprises a user in a
way that requires a ritual to avoid, that is a bug: fix it in code, or file
it with enough context to act on later. A skill that teaches people to route
around a defect keeps the defect.

## Full reference

<https://htpolynet.readthedocs.io/> — in particular the user guide's
*Building a System, Start to Finish*, which this file condenses, and
*Molecular Structure Inputs*, which explains the active form with pictures.
