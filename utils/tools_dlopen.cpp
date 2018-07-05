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
      default:
        return MlfException::what();
    }
  }
  int LibraryFun::init(string path, string funPrefix, FunType t) {
    if(handle) {
      dlclose(handle);
      handle = nullptr;
    }
    handle = dlopen(path.c_str(), RTLD_LAZY);
    if(!handle)
      throw DlException(DlErrorType::FilePathError);

  }
  LibraryFun::~LibraryFun() {
    if(data) {
      if(ffree)
        ffree(data);
      else
        free(data);
    }
    if(handle) {
      dlclose(handle);
    }
    if(object) {
      mlf_dealloc(object);
    }
  }
}

