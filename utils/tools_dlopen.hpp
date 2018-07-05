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

  typedef enum {OptimFun, BasisFun} FunType;
  class LibraryFun {
    protected:
      void *handle = nullptr;
      void *data = nullptr;
      mlf_free_fun ffree = nullptr;
      char *description[mlf_FIELDS+1];
      MLF_OBJ *object = nullptr;
    public:
      int init(string path, string funPrefix, FunType t = OptimFun);
      ~LibraryFun();
  };
}

