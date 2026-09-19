# many_circles
Numeric solution for Helmholtz equation on circles (master homework)

## Build with Nix

The flake builds the C single-circle solver, GSL, and the bundled Fortran
complex Bessel library on Linux. Initialize the library submodule first:

```sh
git submodule update --init --recursive
nix build path:.
./result/bin/many-circles
```

The executable writes `datas_real.dat`, `datas_imag.dat`, and `datas_abs.dat`
in the current directory. It locates its installed Bessel library independently
of that directory. `nix run path:.` builds and runs the same executable, and
`nix flake check path:.` also runs a smoke test of all three output files.
Use `nix develop path:.` for C, Fortran, CMake, GSL, and gnuplot development tools.
After the new flake files are tracked in Git, the `path:` prefix is optional;
the flake automatically includes the `complex_bessel` submodule.

For a local development build inside `nix develop path:.`:

```sh
cmake -S complex_bessel -B complex_bessel/build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_LIBDIR=lib
cmake --build complex_bessel/build
cc -std=gnu11 -O2 -Wall -Wextra src/main.c src/bessel_func.c src/circle.c \
  src/single_layer.c src/single_circle_solv.c -o many_circles \
  -lgsl -lgslcblas -lm -ldl
./many_circles
```

The multi-circle code in `src/many_circle_solv.c` is unfinished and is not part
of the original demo or this package.

## Reference
- Calculation of Complex Valued Bessel Functions: https://dl.acm.org/doi/10.1145/7921.214331
- Several Properties for Bessel Functions: https://dlmf.nist.gov/10.6
- Graf's Addition Theorem: https://wikiwaves.org/Graf%27s_Addition_Theorem
