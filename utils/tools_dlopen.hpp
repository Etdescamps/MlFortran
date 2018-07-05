#pragma once
#include <cstdint>
#include <string>
#include <dlfcn.h>
#include "mlf.hpp"


namespace ToolsDlopen {

  using std::string;
  using MlFortran::MlfException;

  enum class DlErrorType {FilePathError, InvalidDll, MissingFunctions, DataError, InvalidFunctionType};
  class DlException : MlfException {
    public:
      DlErrorType errorDl;
      DlException(DlErrorType e) : MlfException(mlf_OTHERERROR), errorDl(e) {}
      const char* what() const noexcept override;

  };

  class DlLoader {
    protected:
      void *handle = nullptr;
    public:
      void init(string path);
      template<typename FType>
      FType getSymOrNull(string name) {
        return (FType) dlsym(handle, name.c_str());
      }
      template<typename FType>
      FType getSym(string name) {
        void *address = dlsym(handle, name.c_str());
        if(!address)
          throw DlException(DlErrorType::MissingFunctions);
        return (FType) address;
      }
      ~DlLoader();
  };

  enum class LibraryFunType {OptimFun, BasisFun};
  class LibraryFun : protected DlLoader{
    protected:
      void *data = nullptr;
      mlf_free_fun ffree = nullptr;
      char *description[mlf_FIELDS+1];
      MLF_OBJ *object = nullptr;
    public:
      void init(string path, string funPrefix, LibraryFunType typeFun, string fileName = "");
      ~LibraryFun();
  };
}

