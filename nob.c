#define NOB_IMPLEMENTATION
#include "nob.h"

int main(int argc, char** argv) {
    NOB_GO_REBUILD_URSELF(argc, argv);

    Nob_Cmd cmd = {0};
    nob_cmd_append(&cmd, "cc", "-std=gnu11", "-Wall", "-Wextra",
        "-IC:/msys64/ucrt64/include", "-LC:/msys64/ucrt64/lib", "-lgsl", "-lgslcblas", "-lm",
        "-o", "main",
        "src/main.c",
        "src/bessel_func.c",
        "src/circle.c",
        "src/single_layer.c",
        "src/sigle_circle_solv.c");

    if (!nob_cmd_run_sync(cmd)) return 1;
    return 0;
}
