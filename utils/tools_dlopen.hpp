#pragma once
#include <cstdint>
#include <string>
#include <exception>
#include <dlfcn.h>
#include "mlf_funintf.h"


namespace ToolsDlopen {

  using std::string;
  using std::exception;

  class  LibraryException : exception {

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

