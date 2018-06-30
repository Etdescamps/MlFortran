#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static int _print_usage(int out, char *name) {
	fprintf(stderr, "Usage: %s -m model.so (-p m_parameter.dat) (-c checkpoint.h5)\n", name);
  fprintf(stderr, "  -m model.so           Model used as optimisation function\n");
  fprintf(stderr, "  -p m_parameter.dat    Parameter file used for model initialisation\n");
  fprintf(stderr, "  -c checkpoint.h5      Checkpoint file that save the optimisation\n");
  fprintf(stderr, "                        Restart at the last saved state if provided\n");
  return out;
}

int main(int argc, char **argv) {
  return 0;
}

