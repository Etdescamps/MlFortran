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

  void DlLoader::init(string path) {
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
  void LibraryFun::init(string path, string funPrefix, LibraryFunType t) {
    DlLoader::init(path);
  }
  LibraryFun::~LibraryFun() {
    if(data) {
      if(ffree)
        ffree(data);
      else
        free(data);
    }
    if(object) {
      mlf_dealloc(object);
    }
  }
}

