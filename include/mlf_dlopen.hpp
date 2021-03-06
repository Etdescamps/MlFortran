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

#pragma once
#include <cstdint>
#include <string>
#include <dlfcn.h>
#include "mlf.hpp"


namespace MlFortran {
  enum class MlfDlErrorType {FilePathError, InvalidDll, MissingFunctions, DataError, InvalidFunctionType};
  class MlfDlException : MlfException {
    public:
      MlfDlErrorType errorDl;
      MlfDlException(MlfDlErrorType e) : MlfException(mlf_OTHERERROR), errorDl(e) {}
      const char* what() const noexcept override;

  };
  
  enum class MlfLibraryFunType {OptimFun, BasisFun};

  class MlfDlLoader {
    protected:
      void *handle = nullptr;
    public:
      MlfDlLoader() {}
      MlfDlLoader(const string &path) {
        init(path);
      }
      void init(const string &path);
      template<typename FType>
      FType getSymOrNull(const string &name) {
        return (FType) dlsym(handle, name.c_str());
      }
      template<typename FType>
      FType getSym(const string &name) {
        void *address = dlsym(handle, name.c_str());
        if(!address)
          throw MlfDlException(MlfDlErrorType::MissingFunctions);
        return (FType) address;
      }
      ~MlfDlLoader();
  };

  class MlfFunObject : public MlfObject {
    protected:
      void *data = nullptr;
      mlf_free_fun ffree = nullptr;
      const char *description[mlf_FIELDS+1] = {nullptr};
    public:
      MlfFunObject(MlfDlLoader &dl, const string &funPrefix, MlfLibraryFunType typeFun, const string &fileName = "", int nIn = -1, int nOut = -1);
      ~MlfFunObject();
  };
}

