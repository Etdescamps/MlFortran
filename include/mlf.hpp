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

#pragma once
// MlFortran C++ interface
// Implementation is on src/cpp_intf/mlfcpp.cpp

// Standard C Library headers
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <dlfcn.h>

// Standard C++ Library headers
#include <iostream>
#include <string>
#include <exception>
#include <memory>
#include <string_view>

// MLF C interface headers
#include "mlf_cintf.h"
#include "mlf_funintf.h"
#include "mlf_hdf5.h"

namespace MlFortran {
  using std::string;
  using std::string_view;
  using std::exception;
  using std::shared_ptr;
  using std::unique_ptr;
  using std::ostream;


  // Associate C++ datatypes to an MLF_DATATYPE enum
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

  // Return type (return by reference if readOnly)
  template<typename T, bool readOnly>
  struct MlfRefType {
    typedef T Ref;
  };

  template<typename T>
  struct MlfRefType<T, false> {
    typedef T& Ref;
  };

  // Main Exception class of MlFortran
  // It uses the error codes from the Fortran interface
  class MlfException : exception {
    public:
      MLF_ERRORTYPE errorType;
      MlfException(MLF_ERRORTYPE e) : errorType(e) {}
      MlfException() : errorType(mlf_OTHERERROR) {}
      const char* what() const noexcept override;
  };

  class MlfFinalizeCopy :  MlfException {
    public:
      MlfFinalizeCopy() : MlfException() {}
      const char* what() const noexcept override;
  };


  class MlfOutOfBounds : MlfException {
    public:
      MlfOutOfBounds() : MlfException() {}
      const char* what() const noexcept override;
  };

  enum class MlfRscErrorType : short {NotWritable, NotReadable, InvalidAccessType, NotAllocated, NotFound};
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
        // std::cerr << "unallocate object" << std::endl;
        if(obj)
          k = mlf_dealloc(obj);
        if(k<0 && logger)
          logger->reportError("Error deallocating", k);
      }
  };

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
        data = (Type*) calloc(s, sizeof(Type));
        isAllocated = true;
        if(d) {
          memcpy((void*) data, (void*) d, s*sizeof(Type));
          return MlfDataAccessType::InDirect;
        }
        else
          return MlfDataAccessType::NotLoaded;
      }
    public:
      static const MLF_DATATYPE dataType = MlfDataType<Type>::t;
      Ref operator[](const size_t i) const {
        return data[i];
      }
      Type *getData() const {
        return data;
      }
      ~MlfData() {
        _cleanAlloc();
      }
      bool associated() {
        return data != nullptr;
      }
  };

  template<typename Type, bool readOnly = true>
  class MlfDataVector : public MlfData<Type,readOnly> {
    protected:
      size_t dim;
    public:
      typedef MlfData<Type,readOnly> Super;
      const static int rank = 1;
      size_t getSize() const {return dim;}
      MlfDataAccessType putData(Type *d, int dim0[], MLF_ACCESSTYPE at) {
        dim = dim0[0];
        return Super::_putData(d, getSize(), at);
      }
      // WARNING:
      // Fortran geometry (from 1 to N)
      // Not the same as the convention used by the operator[]
      auto operator()(const size_t i) const {
#ifdef _MLF_BOUND_CHECK
        if(!Super::data)
          throw MlfRessourceError(MlfRscErrorType::NotAllocated);
        if(i <= 0 || i > dim)
          throw MlfOutOfBounds();
#endif
        return Super::operator[](i-1);
      }
  };

  template<typename Type, bool readOnly = true>
  ostream& operator<<(ostream &os, const MlfDataVector<Type,readOnly> &v) {
    for(size_t i=1; i <= v.getSize(); ++i)
      os << v(i) << ";";
    return os;
  }

  // WARNING:
  // THIS CLASS USES FORTRAN ARRAY CONVENTION!!!!!
  template<typename Type, bool readOnly = true>
  class MlfDataMatrix : public MlfData<Type,readOnly> {
    protected:
      size_t dims[2];
    public:
      typedef MlfData<Type,readOnly> Super;
      const static int rank = 2;
      size_t getSize() const {return dims[0]*dims[1];}
      MlfDataAccessType putData(Type *d, int dim0[], MLF_ACCESSTYPE at) {
        dims[0] = dim0[0]; dims[1] = dim0[1];
        return Super::_putData(d, getSize(), at);
      }
      auto operator()(const size_t i, const size_t j) const {
#ifdef _MLF_BOUND_CHECK
        if(!Super::data)
          throw MlfRessourceError(MlfRscErrorType::NotAllocated);
        if(i <= 0 || i >= dims[0])
          throw MlfOutOfBounds();
        if(j <= 0 || j >= dims[1])
          throw MlfOutOfBounds();
#endif
        return Super::operator[](i-1+(j-1)*dims[0]);
      }
  };

  // WARNING:
  // THIS CLASS USES FORTRAN ARRAY CONVENTION!!!!!
  template<typename Type, bool readOnly = true>
  class MlfData3DMatrix : public MlfData<Type,readOnly> {
    protected:
      size_t dims[3];
    public:
      typedef MlfData<Type,readOnly> Super;
      const static int rank = 3;
      size_t getSize() const {return dims[0]*dims[1]*dims[2];}
      MlfDataAccessType putData(Type *d, int dim0[], MLF_ACCESSTYPE at) {
        dims[0] = dim0[0]; dims[1] = dim0[1]; dims[2] = dim0[2];
        return Super::_putData(d, getSize(), at);
      }
      auto operator()(const size_t i, const size_t j, const size_t k) const {
#ifdef _MLF_BOUND_CHECK
        if(!Super::data)
          throw MlfRessourceError(MlfRscErrorType::NotAllocated);
        if(i <= 0 || i >= dims[0])
          throw MlfOutOfBounds();
        if(j <= 0 || j >= dims[1])
          throw MlfOutOfBounds();
        if(k <= 0 || k >= dims[2])
          throw MlfOutOfBounds();
#endif
        return Super::operator[](i-1+(j-1)*dims[0]);
      }
  };

  class MlfObject {
    private:
      // Make this kind of object non copiable (required by MLF_OBJ Fortran interface)
      MlfObject(const MlfObject&) = delete;
      MlfObject &operator=(const MlfObject&) = delete;
    protected:
      MLF_OBJ *obj = nullptr;
      unique_ptr<string_view[]> obj_names;
      int nrsc;
      bool is_copy = 0;
      bool updateIdMap();
    public:
      MlfObject() = default;
      MlfObject(MLF_OBJ *obj) :  obj(obj) {}
      int setObject(MLF_OBJ *object, bool isCopy = false) {
        if(!object)
          return -1;
        obj = object;
        is_copy = isCopy;
        return 0;
      }
      ~MlfObject() {
        if(!is_copy)
          finalize();
      }
      virtual void finalize() {
        if(is_copy)
          throw MlfFinalizeCopy();
        if(obj)
          mlf_dealloc(obj);
        obj = nullptr;
      }
      MLF_OBJ *get() const {
        return obj;
      }
      const char *getName(int id) const {
        return mlf_getinfo(obj, id, mlf_NAME);
      }
      const char *getDesc(int id) const {
        return mlf_getinfo(obj, id, mlf_DESC);
      }
      const char *getFields(int id) const {
        return mlf_getinfo(obj, id, mlf_FIELDS);
      }
      void *getRsc(int id, MLF_DT &dt, int &rank, int dims[], void *data=nullptr) const {
        return mlf_getrsc(obj, id, &dt, &rank, dims, data);
      }
      int updateRsc(int id, void *data) const {
        return mlf_updatersc(obj, id, data);
      }
      int getNumRsc() const {
        return mlf_getnumrsc(obj);
      }
      template<typename MData>
      MlfDataAccessType getRsc(int id, MData &md) const {
        typedef typename MData::Type T;
        int rank, dims[MLF_MAXRANK];
        MLF_DT dt;
        void *data = getRsc(id, dt, rank, dims);
        if(dt.dt != MData::dataType)
          throw MlfException(mlf_WRONGTYPE);
        if(rank != MData::rank)
          throw MlfException(mlf_WRONGRANK);
        MlfDataAccessType at = md.putData((T*) data, dims, (MLF_ACCESSTYPE) dt.access);
        
        if(at == MlfDataAccessType::NotLoaded) {
          data = getRsc(id, dt, rank, dims, (void*) md.getData());
        }
        return at;
      }
      int getIdName(const string &name, bool updated = false);
  };

  class MlfStepObject : public MlfObject {
    protected:
      MlfDataVector<int64_t> idata;
      MlfDataVector<double> rdata;
      void initOutput();
      bool is_initData = false;
    public:
      using MlfObject::MlfObject;
      int step() const {
        double dt = 0;
        int64_t niter = 1;
        return mlf_step(get(), &dt, &niter);
      }
      int step(double &dt) const {
        int64_t niter = 1;
        return mlf_step(get(), &dt, &niter);
      }
      int step(int64_t &nstep) const {
        double dt = 0;
        return mlf_step(get(), &dt, &nstep);
      }
      int step(double &dt, int64_t &nstep) const {
        return mlf_step(get(), &dt, &nstep);
      }
      void printLine(ostream& os);
      void printFields(ostream& os);
  };

  class MlfHdf5 : public MlfObject {
    protected:
      bool has_data = false, is_init = false;
      int overwrite = 0;
    public:
      int createFile(string &fileName, bool trunk = true) {
        int r = setObject(mlf_hdf5_createFile(fileName.c_str(), trunk ? 1 : 0));
        std::cout << "Creating " << fileName << r << std::endl;
        is_init = r >= 0;
        overwrite = 0;
        has_data = false;
        return r;
      }
      int openFile(string &fileName, bool rw = true) {
        int r = setObject(mlf_hdf5_openFile(fileName.c_str(), rw ? 1 : 0));
        std::cout << "Opening " << fileName << r << std::endl;
        is_init = r >= 0;
        overwrite = 1;
        has_data = is_init;
        return r;
      }
      int pushState(MlfObject &elt) {
        int r = mlf_pushState(get(), elt.get(), overwrite);
        std::cout << "PushState " << overwrite << std::endl;
        if(r >= 0)
          overwrite = 1;
        return r;
      }
      bool readWOrCreate(string &fileName);
      bool hasData() {return has_data;}
      bool isInit() {return is_init;}
  };

  class MlfOptimObject : public MlfStepObject {
    protected:
      MlfDataVector<double> minX;
    public:
      int setAlgo(const string &nalg, MlfObject &funobj, double target, int lambda, int mu, double sigma = 1.0) {
        return setObject(mlf_getoptimobj(nalg.c_str(), funobj.get(), nullptr, target, lambda, mu, sigma));
      }
      int setAlgo(const string &nalg, MlfObject &funobj, MlfHdf5 &handler, double target, int lambda, int mu, double sigma = 1.0) {
        return setObject(mlf_getoptimobj(nalg.c_str(), funobj.get(), handler.hasData() ? handler.get() : nullptr, target, lambda, mu, sigma));
      }
      void printMinX(ostream& os);
  };
}

