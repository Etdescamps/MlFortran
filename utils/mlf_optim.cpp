#include <iostream>
#include <cstdlib>
#include <string>
#include <unistd.h>

using namespace std;



static int _print_usage(int out, char *name) {
  cerr << "Usage: "<< name <<" -m model.so (-p m_parameter.dat) (-c checkpoint.h5)" << endl;
  cerr << "  -m model.so           Model used as optimisation function" << endl;
  cerr << "  -p m_parameter.dat    Parameter file used for model initialisation" << endl;
  cerr << "  -c checkpoint.h5      Checkpoint file that save the optimisation" << endl;
  cerr << "                        Restart at the last saved state if provided" << endl;
  cerr << "  -i number_of_input    Specify (if required) the number of input (dimension of the search dataspace)" << endl;
  cerr << "  -o number_of_output   Specify (if required) the number of objective" << endl;
  cerr << "  -h                    Display this help message" << endl;
  return out;
}

int main(int argc, char **argv) {
  string nmodel, nparameter, ncheckpoint;
  int opt, ninput = -1, noutput = -1;
  while((opt = getopt(argc, argv, "m:p:c:hi:o:")) != -1) {
    switch(opt) {
      case 'm':
        nmodel = optarg;
        break;
      case 'p':
        nparameter = optarg;
        break;
      case 'c':
        ncheckpoint = optarg;
        break;
      case 'i':
        ninput = atoi(optarg);
        break;
      case 'o':
        noutput = atoi(optarg);
        break;
      case 'h':
        return _print_usage(0, argv[0]);
      default:
        return _print_usage(-1, argv[0]);
    }
  }
  if(nmodel.size() == 0)
    return _print_usage(-1, argv[0]);
  return 0;
}

