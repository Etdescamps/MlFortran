/*-
 * Copyright (c) 2017-2018 Etienne Descamps
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

#include <dlfcn.h>
#include <cstdlib>
#include "tools_dlopen.hpp"

namespace ToolsDlopen {
  const char *DlException::what() const noexcept {
    switch(errorDl) {
      case DlErrorType::FilePathError:
        return "DlException: Incorrect file path";
      case DlErrorType::InvalidDll:
        return "DlException: Invalid library file";
      case DlErrorType::MissingFunctions:
        return "DlException: Library file is missing functions";
      case DlErrorType::DataError:
        return "DlException: Error when handling data";
      case DlErrorType::InvalidFunctionType:
        return "DlException: Invalid function type";
      default:
        return MlfException::what();
    }
  }

  void DlLoader::init(const string &path) {
    if(handle) {
      dlclose(handle);
      handle = nullptr;
    }
    handle = dlopen(path.c_str(), RTLD_LAZY);
    if(!handle)
      throw DlException(DlErrorType::FilePathError);
  }

  DlLoader::~DlLoader(){
    if(handle) {
      dlclose(handle);
    }
  }

  MLF_OBJ *LibraryFun::getFunObj(const string &path, const string &funPrefix, LibraryFunType typeFun, const string &fileName, int nIn, int nOut) {
    MLF_OBJ *object = nullptr;
    DlLoader::init(path);
    mlf_init_fun finit = DlLoader::getSymOrNull<mlf_init_fun>(funPrefix+"_init");
    ffree = DlLoader::getSymOrNull<mlf_free_fun>(funPrefix+"_free");
    mlf_getinfo_fun finfo = DlLoader::getSymOrNull<mlf_getinfo_fun>(funPrefix+"_getinfo");
    if(finit)
      data = finit(fileName.c_str());
    void *info = nullptr;
    if(finfo) {
      for(int i=mlf_NAME; i <= mlf_FIELDS; i++)
        description[i] = (char *) finfo(data, i);
      finfo(data, mlf_FUNINFO);
    }
    switch(typeFun) {
      case LibraryFunType::OptimFun:
        {
          mlf_objective_fun fobj = DlLoader::getSym<mlf_objective_fun>(funPrefix+"_objfun");
          mlf_objective_fun fcstr = DlLoader::getSymOrNull<mlf_objective_fun>(funPrefix+"_cstrfun");
          MLF_OBJFUNINFO nfo = { nIn, -1, nOut};
          if(info) {
            // Get default values of the model
            nfo = *(MLF_OBJFUNINFO*) info;
            // Use input dimension if provided
            if(nIn>0)
              nfo.nDimIn = nIn;
            if(nOut>0)
              nfo.nDimOut = nOut;
          }
          if(nfo.nDimOut<0) // Not initialised number of output dimensions
            nfo.nDimOut = 1;
          object = mlf_objfunction(fobj, data, fcstr, &nfo);
        }
        break;
      case LibraryFunType::BasisFun:
        {
          mlf_basis_fun fbasis = DlLoader::getSym<mlf_basis_fun>(funPrefix+"_basisfun");
          object = mlf_basisfunction(fbasis, data);
        }
        break;
      default:
        throw DlException(DlErrorType::InvalidFunctionType);
    }
    return object;
  }

  LibraryFun::~LibraryFun() {
    if(data) {
      if(ffree)
        ffree(data);
      else
        free(data);
    }
  }
}

