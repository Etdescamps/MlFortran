/*-
 * Copyright (c) 2017-2019 Etienne Descamps
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation and/or
 *    other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be
 *    used to endorse or promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <limits>
#include "mlf_dlopen.hpp"

using namespace std;
using namespace MlFortran;

int _proceed_optim(string &nalg, MlfObject &fobj, int64_t niter, double target, int lambda, int mu, double sigma, MlfHdf5 &handler, int cEvery) {
  MlfOptimObject obj_optim;
  obj_optim.setAlgo(nalg, fobj, handler, target, lambda, mu, sigma);
  obj_optim.printFields(cout);
  for(int64_t i = 0; i < niter; ++i) {
    obj_optim.printLine(cout);
    if(handler.isInit() && ((i %cEvery) == 0))
      handler.pushState(obj_optim);
    if(obj_optim.step() != 0)
      break;
  }
  if(handler.isInit())
    handler.pushState(obj_optim);
  obj_optim.printLine(cout);
  obj_optim.printMinX(cout);
  return 0;
}

static int _print_usage(int out, char *name) {
  cerr << "Usage: "<< name <<" -m model.so (-X prefix) (-p m_parameter.dat) (-c checkpoint.h5)" << endl;
  cerr << "  -m model.so           Model used as optimisation function" << endl;
  cerr << "  -X prefix_name        Prefix name used by the function" << endl;
  cerr << "  -p m_parameter.dat    Parameter file used for model initialisation" << endl;
  cerr << "  -c checkpoint.h5      Checkpoint file that save the optimisation" << endl;
  cerr << "  -C checkpointEvery    Set the number of steps between each checkpoint" << endl;
  cerr << "  -A algorithm          Select an algorithm for optimization" << endl;
  cerr << "                        Possible algorithms: cmaes (default), maes" << endl;
  cerr << "                        Restart at the last saved state if provided" << endl;
  cerr << "  -i number_of_input    Specify (if required) the number of input (dimension of the search dataspace)" << endl;
  cerr << "  -o number_of_output   Specify (if required) the number of objective (default 1)" << endl;
  cerr << "  -L lambda             Specify the offspring population size" << endl;
  cerr << "  -M mu                 Specify the selected number of individuals" << endl;
  cerr << "  -N nstep              Maximum number of iterations" << endl;
  cerr << "  -S sigma              Initial step size (default 1.0)" << endl;
  cerr << "  -T target             Target value" << endl;
  cerr << "  -h                    Display this help message" << endl;
  return out;
}

int main(int argc, char **argv) {
  string nmodel, nparameter, ncheckpoint;
  string nalg = "cmaes", nprefix = "fun";
  int opt, ninput = -1, noutput = -1, cEvery = 16;
  int64_t niter = numeric_limits<int64_t>::max();
  int mu = -1, lambda = -1;
  double sigma = 1.0, target = -HUGE_VAL;
  try {
    while((opt = getopt(argc, argv, "m:X:p:c:C:A:i:o:L:M:N:T:S:h")) != -1) {
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
        case 'C':
          cEvery = stoi(optarg);
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
        case 'N':
          niter = stoll(optarg);
          break;
        case 'S':
          sigma = stod(optarg);
          break;
        case 'T':
          target = stod(optarg);
          break;
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
  mlf_init();
  {
    MlfDlLoader lib(nmodel);
    MlfHdf5 handler;
    if(ncheckpoint.size() > 0)
      handler.readWOrCreate(ncheckpoint);
    MlfFunObject obj(lib, nprefix, MlfLibraryFunType::OptimFun, nparameter, ninput, noutput);
    _proceed_optim(nalg, obj, niter, target, lambda, mu, sigma, handler, cEvery);
    handler.finalize();
  }
  mlf_quit();
  return 0;
}

