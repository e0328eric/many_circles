{
  description = "Circular Helmholtz solver in C with GSL and complex Bessel functions";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    self.submodules = true;
  };

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" ];
      forAllSystems = nixpkgs.lib.genAttrs systems;
    in
    {
      packages = forAllSystems (system:
        let
          pkgs = import nixpkgs { inherit system; };
          package = pkgs.stdenv.mkDerivation {
            pname = "many-circles";
            version = "0.1.0";

            src = pkgs.lib.cleanSourceWith {
              src = self;
              filter = path: type:
                let name = baseNameOf path;
                in pkgs.lib.cleanSourceFilter path type
                  && !(builtins.elem name [ "result" "build" ".cache" ])
                  && !(builtins.elem name [ "main" "nob" "nob.old" ])
                  && !(pkgs.lib.hasSuffix ".dat" name);
            };

            nativeBuildInputs = [ pkgs.cmake pkgs.gfortran ];
            buildInputs = [ pkgs.gsl ];
            dontConfigure = true;

            # The original demo resolves the Fortran library relative to cwd.
            # An installed executable must also work in an output directory.
            postPatch = ''
              substituteInPlace src/main.c \
                --replace-fail './complex_bessel/build/lib/libcomplex_bessel.so' \
                  "$out/lib/libcomplex_bessel.so"
            '';

            buildPhase = ''
              runHook preBuild
              cmake -S complex_bessel -B bessel-build \
                -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_LIBDIR=lib \
                -DBUILD_TESTING=OFF
              cmake --build bessel-build --parallel "$NIX_BUILD_CORES"
              $CC -std=gnu11 -O2 -Wall -Wextra \
                src/main.c src/bessel_func.c src/circle.c \
                src/single_layer.c src/single_circle_solv.c \
                -o many-circles -lgsl -lgslcblas -lm -ldl
              runHook postBuild
            '';

            installPhase = ''
              runHook preInstall
              install -Dm755 many-circles "$out/bin/many-circles"
              install -Dm755 bessel-build/lib/libcomplex_bessel.so \
                "$out/lib/libcomplex_bessel.so"
              install -Dm644 README.md "$out/share/doc/many-circles/README.md"
              install -Dm644 plotting_abs.gp plotting_imag.gp plotting_real.gp \
                -t "$out/share/many-circles"
              runHook postInstall
            '';

            doInstallCheck = true;
            installCheckPhase = ''
              runHook preInstallCheck
              mkdir smoke-test
              cd smoke-test
              "$out/bin/many-circles"
              for data in datas_real.dat datas_imag.dat datas_abs.dat; do
                test -s "$data"
                awk 'NR > 1 { if (NF != 3 || tolower($0) ~ /nan|inf/) exit 1; count++ }
                     END { if (count < 30000) exit 1 }' "$data"
              done
              cd ..
              runHook postInstallCheck
            '';

            meta = {
              description = "Single-circle Helmholtz transmission solver";
              mainProgram = "many-circles";
              platforms = systems;
            };
          };
        in {
          default = package;
          many-circles = package;
        });

      apps = forAllSystems (system: {
        default = {
          type = "app";
          program = "${self.packages.${system}.default}/bin/many-circles";
          meta.description = "Run the single-circle Helmholtz demo";
        };
      });

      checks = forAllSystems (system: {
        solver = self.packages.${system}.default;
      });

      devShells = forAllSystems (system:
        let pkgs = import nixpkgs { inherit system; };
        in {
          default = pkgs.mkShell {
            inputsFrom = [ self.packages.${system}.default ];
            packages = [ pkgs.gnuplot ];
          };
        });
    };
}
