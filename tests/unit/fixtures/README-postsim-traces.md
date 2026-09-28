# Real postsim traces, kept as test fixtures

Two outputs of actual htpolynet runs, not synthetic data.

- `ladder-FDE-DFDA.csv` — a temperature ladder, 61 points from 300 to 599 K, with
  density and the glassy/rubbery lines that were fitted to it at the time.  From the
  furan series, system 1 (FDE-DFDA).
- `deform-x-FDE-FCPDA.csv` — a uniaxial deformation, 1001 frames, written by
  `PostSimDeform` itself: `Box-X`, `Pres-XX` and the `-strain` / `-stress` columns the
  stage derives.  From the furan series, system 3 (FDE-FCPDA), replicate 3.

They exist because the synthetic tests cannot show what these fits do to real data.
A clean straight line recovers its slope to eight digits; an MD trace does not.  On
these two traces the answer moves further with the fit window than most people would
guess, which is the thing worth regression-testing.

There is deliberately **no shear fixture**: no shear simulation has ever been run.
When one is, a `shear-xy.csv` belongs here beside these.
