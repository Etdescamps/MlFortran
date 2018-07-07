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
// MlFortran C++ interface
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <string>
#include <exception>
#include <memory>
#include <string_view>
#include <dlfcn.h>
#include "mlf_cintf.h"
#include "mlf_funintf.h"

namespace MlFortran {
  using std::string;
  using std::string_view;
  using std::exception;
  using std::shared_ptr;

  template<typename T>
  struct MlfDataType {
    static const MLF_DATATYPE t = mlf_RAW;
  };
  template<> struct MlfDataType<int> {
    static const MLF_DATATYPE t = mlf_INT;
  };
  template<> struct MlfDataType<int64_t> {
    static const MLF_DATATYPE t = mlf_INT64;
  };
  template<> struct MlfDataType<bool> {
    static const MLF_DATATYPE t = mlf_BOOL;
  };
  template<> struct MlfDataType<size_t> {
    static const MLF_DATATYPE t = mlf_SIZE;
  };
  template<> struct MlfDataType<float> {
    static const MLF_DATATYPE t = mlf_FLOAT;
  };
  template<> struct MlfDataType<double> {
    static const MLF_DATATYPE t = mlf_DOUBLE;
  };

  template<typename T, bool readOnly>
  struct MlfRefType {
    typedef T Ref;
  };

  template<typename T>
  struct MlfRefType<T, false> {
    typedef T& Ref;
  };


  class MlfException : exception {
    public:
      MLF_ERRORTYPE errorType;
      MlfException(MLF_ERRORTYPE e) : errorType(e) {}
      MlfException() : errorType(mlf_OTHERERROR) {}
      const char* what() const noexcept override;
  };

  class MlfOutOfBounds : MlfException {
    public:
      MlfOutOfBounds() : MlfException() {}
      const char* what() const noexcept override;
  };

  enum class MlfRscErrorType : short {NotWritable, NotReadable, InvalidAccessType};
  class MlfRessourceError :  MlfException {
    public:
      MlfRscErrorType rscError;
      MlfRessourceError(MlfRscErrorType t) : MlfException(), rscError(t) {}
      const char* what() const noexcept override;
  };


  class MlfLogger {
    public:
      virtual void reportError(const char *msg, int error) = 0;
  };

  class MlfDeleter {
    protected:
      MlfLogger *logger = nullptr;
    public:
      void operator()(MLF_OBJ *obj) const {
        int k = 0;
        if(obj)
          k = mlf_dealloc(obj);
        if(k<0 && logger)
          logger->reportError("Error deallocating", k);
      }
  };

  typedef shared_ptr<MLF_OBJ> MlfShared;

  enum class MlfDataAccessType : short {Direct, InDirect, NotLoaded, WriteOnly};

  template<typename T, bool readOnly = true>
  class MlfData {
    public:
      typedef T Type;
      typedef typename MlfRefType<T, readOnly>::Ref Ref;
    protected:
      Type *data = nullptr;
      bool isAllocated = false;
      MLF_ACCESSTYPE accessType;
      size_t size;
      
      void _cleanAlloc() {
        if(isAllocated && data)
          free(data);
      }
      MlfDataAccessType _putData(Type *d, size_t s, MLF_ACCESSTYPE at) {
        _cleanAlloc();
        accessType = at;
        size = s;
        if(!d) {
          return _allocData(s);
        }
        switch(at) {
          case mlf_READONLY:
            if(!readOnly)
              throw MlfRessourceError(MlfRscErrorType::NotWritable);
          case mlf_INDIRECT:
            if(!readOnly)
              return _allocData(s, d);
          case mlf_DIRECT:
            data = d;
            return MlfDataAccessType::Direct;
          case mlf_WRITEONLY:
            if(readOnly)
              throw MlfRessourceError(MlfRscErrorType::NotReadable);
            data = d;
            return MlfDataAccessType::WriteOnly;
          case mlf_COPYONLY:
            return _allocData(s, d);
          default:
            throw MlfRessourceError(MlfRscErrorType::InvalidAccessType);
        }
      }
      MlfDataAccessType _allocData(size_t s, Type *d = nullptr) {
        data = calloc(s, sizeof(Type));
        isAllocated = true;
        if(d) {
          memcpy(data, d, s*sizeof(Type));
          return MlfDataAccessType::InDirect;
        }
        else
          return MlfDataAccessType::NotLoaded;
      }
    public:
      static const MLF_DATATYPE dt = MlfDataType<Type>::t;
      Ref operator[](const size_t i) {
        return data[i];
      }
      Type *getData() {
        return data;
      }
      ~MlfData() {
        _cleanAlloc();
      }
  };

  template<typename Type, bool readOnly = true>
  class MlfDataVector : public MlfData<Type,readOnly> {
    protected:
      size_t dim;
    public:
      const static int rank = 1;
      size_t getSize() {return dim;}
      MlfDataAccessType putData(Type *d, int dim0[], MLF_ACCESSTYPE at) {
        dim = dim0[0];
        return _putData(d, getSize(), at);
      }
      // WARNING:
      // Fortran geometry (from 1 to N)
      // Not the same as this of the operator[]
      auto operator()(const size_t i) {
        if(i <= 0 || i > dim)
          throw MlfOutOfBounds();
        return MlfData<Type>::operator[](i-1);
      }
  };

  // WARNING:
  // THIS CLASS USE FORTRAN ARRAY CONVENTION!!!!!
  template<typename Type, bool readOnly = true>
  class MlfDataMatrix : public MlfData<Type,readOnly> {
    protected:
      size_t dims[2];
    public:
      const static int rank = 2;
      size_t getSize() {return dims[0]*dims[1];}
      MlfDataAccessType putData(Type *d, int dim0[], MLF_ACCESSTYPE at) {
        dims[0] = dim0[0]; dims[1] = dim0[1];
        return _putData(d, getSize(), at);
      }
      auto operator()(const size_t i, const size_t j) {
        if(i <= 0 || i >= dims[0])
          throw MlfOutOfBounds();
        if(j <= 0 || j >= dims[1])
          throw MlfOutOfBounds();
        return MlfData<Type>::operator[](i-1+(j-1)*dims[0]);
      }
  };

  // WARNING:
  // THIS CLASS USE FORTRAN ARRAY CONVENTION!!!!!
  template<typename Type, bool readOnly = true>
  class MlfData3DMatrix : public MlfData<Type,readOnly> {
    protected:
      size_t dims[3];
    public:
      const static int rank = 3;
      size_t getSize() {return dims[0]*dims[1]*dims[2];}
      MlfDataAccessType putData(Type *d, int dim0[], MLF_ACCESSTYPE at) {
        dims[0] = dim0[0]; dims[1] = dim0[1]; dims[2] = dim0[2];
        return _putData(d, getSize(), at);
      }
      auto operator()(const size_t i, const size_t j, const size_t k) {
        if(i <= 0 || i >= dims[0])
          throw MlfOutOfBounds();
        if(j <= 0 || j >= dims[1])
          throw MlfOutOfBounds();
        if(k <= 0 || k >= dims[2])
          throw MlfOutOfBounds();
        return MlfData<Type>::operator[](i-1+(j-1)*dims[0]);
      }
  };

  class MlfObject {
    protected:
      MlfShared obj;
    public:
      MlfObject(MlfShared &obj) : obj(obj) {}
      MlfObject(MLF_OBJ *obj) : obj(obj,
          [](MLF_OBJ *obj) {if(obj) mlf_dealloc(obj);}) {}
      MlfObject(MLF_OBJ *obj, MlfDeleter &deleter) : obj(obj, deleter) {}
      MLF_OBJ *get() {
        return obj.get();
      }
      const char *getName(int id) {
        return mlf_getinfo(obj.get(), id, mlf_NAME);
      }
      const char *getDesc(int id) {
        return mlf_getinfo(obj.get(), id, mlf_DESC);
      }
      const char *getFields(int id) {
        return mlf_getinfo(obj.get(), id, mlf_FIELDS);
      }
      void *getRsc(int id, MLF_DT &dt, int &rank, int dims[], void *data=nullptr) {
        return mlf_getrsc(obj.get(), id, &dt, &rank, dims, data);
      }
      int updateRsc(int id, void *data) {
        return mlf_updatersc(obj.get(), id, data);
      }
      template<typename MData>
      MlfDataAccessType getRsc(int id, MData &md) {
        typedef typename MData::Type T;
        int rank, dims[MLF_MAXRANK];
        MLF_DT dt;
        void *data = getRsc(id, dt, rank, dims);
        if(dt.dt != MData::t)
          throw MlfException(mlf_WRONGTYPE);
        if(rank != MData::rank)
          throw MlfException(mlf_WRONGRANK);
        MlfDataAccessType at = md.putData((T*) data, dims, dt.access);
        
        if(at == MlfDataAccessType::NotLoaded) {
          data = getRsc(id, dt, rank, dims, (void*) md.getData());
        }
        return at;
      }
  };
}

