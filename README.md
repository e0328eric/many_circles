# many_circle

Glosso port of `../many_circles_c`, solving the two-dimensional Helmholtz
transmission problem on circular inclusions. The port uses Glosso's installed
standard library for complex Bessel functions and complex arithmetic; it has
no Fortran, GSL, or external numerical-library dependency.

## Build and run

From the project root, with the installed `glosso` compiler on `PATH`:

```sh
glosso first.glo
glosso first.glo -- run
glosso first.glo -- test
glosso first.glo -- clean
```

`first.glo` runs at compile time using Glosso's native build API; Make is
not needed. The default command (also `-- build` or `-- all`) builds
`build/many_circle` at `-O2`, with the child compiler cache in `build/cache`.
`run` builds and runs the example, `test` builds and runs all four test
executables, and `clean` removes `build/`. Use `-- help` for usage.
Builds invoke the compiler each time, without Make-style timestamp checks.

To compile the application directly:

```sh
mkdir -p build
glosso -O2 --cache-dir build/cache src/main.glo -o build/many_circle
./build/many_circle
```

The compiler discovers its standard library under `~/.local/glosso/std`.
To use another compiler, invoke `/path/to/glosso first.glo`; child builds use
that compiler too. Child build options are set in `first.glo` independently
of the outer compiler's flags. For a custom standard-library path, pass
`--std /path/to/std` to the outer compiler and add the same option with
`Compiler.add_compiler_arg` in `build_executable`.

The default example retains the C program's frequency 5, 11 Fourier modes,
unit circle, unit density and bulk modulus, and incident coefficient
`c_5 = 5`. It writes `datas_real.dat`, `datas_imag.dat`, and `datas_abs.dat`
in the current directory, replacing existing files with those names. Each
contains 31,815 `(x, y, value)` samples, with five decimal places and the same
sampling order as the C demo. The original gnuplot scripts are included:

```sh
gnuplot plotting_abs.gp
```

Run the two-circle material example with:

```sh
./build/many_circle --multi
```

This uses frequency 2.3, 21 modes, radius-0.6 circles centered at `(-1.3, 0)`
and `(1.3, 0)`, and a global incident `J_0` wave. Edit `src/main.glo` to change
the demonstration geometry, material parameters, incident wave, or grid.

## Solver

For each circle, `k = omega * sqrt(rho / kappa)` uses the principal complex
square root. The background density and bulk modulus are both 1. The field
and its density-weighted normal derivative are continuous at each boundary;
scattered waves satisfy the outgoing radiation condition.

The C code represents fields through single-layer densities. This port
stores the corresponding physical expansion coefficients directly:

```text
outside: u_incident + sum_n a_n H_n^(1)(omega*r) exp(i*n*theta)
inside:              sum_n b_n J_n(k*r)         exp(i*n*theta)
```

This cancels the interior Hankel normalization from the C formulation.
Glosso provides complex `J_n`, but its installed standard library does not
provide complex `Y_n` or Hankel functions. Only real positive exterior Hankel
arguments remain, evaluated as `J_n(x) + i Y_n(x)`. Derivatives use the
adjacent-order recurrence. No replacement complex-Bessel implementation or
Fortran binding is needed.

The single-circle solver solves one complex 2-by-2 system per Fourier mode.
The original C multi-circle source is unfinished and excluded from its
build; this port completes that path using Graf translation between disjoint
circles and a dense complex solve with column scaling and scaled row
pivoting. Interior expansions belong only to their own circles.

Frequency and radii must be finite and positive; density and bulk modulus
must be finite and nonzero. Mode counts must be positive and odd. Circles
must not overlap or touch. Invalid inputs, reported special-function errors,
and numerically singular systems stop with a diagnostic. Accuracy depends
on Fourier truncation; increase the mode count and check convergence for
larger frequencies or closely spaced circles.

## Library use

Load `src/module.glo`, or the individual components. `Circle` contains
`center: c128`, `radius: f64`, `rho: c128`, and `kappa: c128`.

- `solve_single_circle(circle, omega, data)` returns an owned
  `Single_Solution`. Its incident coefficients use polar coordinates about
  that circle's center. Evaluate with `get_solution_value(*solution, r, theta)`
  in those local coordinates, or `get_solution_at(*solution, point)` in
  global coordinates. Release with `single_solution_free(*solution)`.
- `solve_multi_circle(circles, omega, data)` returns an owned
  `Multi_Solution`. Here the incident coefficients are about the global
  origin. Evaluate with `get_solution_value_multi(*solution, point)` and
  release with `multi_solution_free(*solution)`.

Both APIs take slices and copy their inputs. For `N` entries, index `j`
represents order `j - N/2`. Each solution's `outside` and `inside` arrays
contain the physical `a_n` and `b_n` above, not the C layer densities.
Multi-circle arrays are grouped by circle, then mode. Owning solution
structs must not be shallow-copied and freed twice. Evaluation supports circle
centers and takes the interior trace at an exact boundary.

## Verification and C build

`glosso first.glo -- test` checks Bessel references and identities, complex-material
transmission conditions, transparent media, translated centers, multi-circle
boundary and flux continuity, circle-order invariance, and agreement with
the single-circle solver. It also runs fixed numerical fixtures generated
from the original C solver using AMOS and GSL; these fixtures require only
Glosso to run.

The C project now has its own flake, retaining its Fortran dependency:

```sh
cd ../many_circles_c
nix build path:.
./result/bin/many-circles
```

From this Glosso project, compare the full default output grids:

```sh
python3 tests/compare_grids.py ../many_circles_c/result/bin/many-circles build/many_circle
```

The comparison uses temporary output directories. To rebuild and rerun the
complex-material reference fixtures against the sibling C source, see
`tests/check_c_reference.sh`. These optional regression checks require Nix
and/or the C executable; the Glosso application itself does not.

The original project's MIT license is retained in `LICENSE`.

## Measured comparison

See [COMPARISON.md](COMPARISON.md) for the C/Glosso build-time, runtime,
binary-size, and numerical-correctness results, including fixes in the sibling
Glosso compiler. [benchmarks/README.md](benchmarks/README.md) explains how to
repeat the measurements and the resonance accuracy check.

See [GENERAL_NUMERICS.md](GENERAL_NUMERICS.md) for the follow-up fixes to shared
complex division and numerical error reporting, with tests beyond Bessel functions.
