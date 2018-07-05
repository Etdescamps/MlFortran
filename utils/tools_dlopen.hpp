#pragma once
#include <cstdint>
#include <string>
#include <dlfcn.h>
#include "mlf.hpp"


namespace ToolsDlopen {

  using std::string;
  using MlFortran::MlfException;

  enum class DlErrorType {FilePathError, InvalidDll, MissingFunctions, DataError};
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
      void init(string path, string funPrefix, LibraryFunType t);
      ~LibraryFun();
  };
}

