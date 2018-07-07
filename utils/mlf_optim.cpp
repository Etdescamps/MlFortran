#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include "tools_dlopen.hpp"

using namespace std;
using namespace ToolsDlopen;
using namespace MlFortran;

int _proceed_optim(const char *nalg, MlfObject &fobj, int lambda, int mu, double sigma) {
  MlfObject obj_optim(mlf_getoptimobj(nalg, fobj.get(), nullptr, lambda, mu, sigma));
  return 0;
}

static int _print_usage(int out, char *name) {
  cerr << "Usage: "<< name <<" -m model.so (-X prefix) (-p m_parameter.dat) (-c checkpoint.h5)" << endl;
  cerr << "  -m model.so           Model used as optimisation function" << endl;
  cerr << "  -X prefix_name        Prefix name used by the function" << endl;
  cerr << "  -p m_parameter.dat    Parameter file used for model initialisation" << endl;
  cerr << "  -c checkpoint.h5      Checkpoint file that save the optimisation" << endl;
  cerr << "  -A algorithm          Select an algorithm for optimization" << endl;
  cerr << "                        Possible algorithms: cmaes (default), maes" << endl;
  cerr << "                        Restart at the last saved state if provided" << endl;
  cerr << "  -i number_of_input    Specify (if required) the number of input (dimension of the search dataspace)" << endl;
  cerr << "  -o number_of_output   Specify (if required) the number of objective (default 1)" << endl;
  cerr << "  -L lambda             Specify the offspring population size" << endl;
  cerr << "  -M mu                 Specify the selected number of individuals" << endl;
  cerr << "  -S sigma              Initial step size (default 1.0)" << endl;
  cerr << "  -h                    Display this help message" << endl;
  return out;
}

int main(int argc, char **argv) {
  string nmodel, nparameter, ncheckpoint;
  string nalg = "cmaes", nprefix = "fun";
  int opt, ninput = -1, noutput = -1;
  int mu = -1, lambda = -1;
  double sigma = 1.0;
  try {
    while((opt = getopt(argc, argv, "m:X:p:c:A:hi:o:L:M:S")) != -1) {
      switch(opt) {
        case 'm':
          nmodel = optarg;
          break;
        case 'X':
          nprefix = optarg;
          break;
        case 'p':
          nparameter = optarg;
          break;
        case 'c':
          ncheckpoint = optarg;
          break;
        case 'A':
          nalg = optarg;
          break;
        case 'i':
          ninput = stoi(optarg);
          break;
        case 'o':
          noutput = stoi(optarg);
          break;
        case 'L':
          lambda = stoi(optarg);
          break;
        case 'M':
          mu = stoi(optarg);
          break;
        case 'S':
          sigma = stod(optarg);
        case 'h':
          return _print_usage(0, argv[0]);
        default:
          cerr << "Unknown argument -"<< opt << endl;
          return _print_usage(-1, argv[0]);
      }
    }
  }
  catch(invalid_argument& e) {
    cerr << "Error option -" << opt << " invalid argument " << optarg << endl;
    return _print_usage(-1, argv[0]);
  }
  catch(out_of_range& e) {
    cerr << "Error option -" << opt << " value too huge " << optarg << endl;
    return -1;
  }
  if(nmodel.size() == 0)
    return _print_usage(-1, argv[0]);
  LibraryFun lib;
  MlfObject obj(lib.init(nmodel, nprefix, LibraryFunType::OptimFun, nparameter, ninput, noutput));
  _proceed_optim(nalg.c_str(), obj, lambda, mu, sigma);
  return 0;
}

