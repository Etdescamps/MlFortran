#include <dlfcn.h>
#include <cstdlib>
#include "tools_dlopen.hpp"

namespace ToolsDlopen {
  int LibraryFun::init(string path, string funPrefix, FunType t) {
    if(handle) {
      dlclose(handle);
      handle = nullptr;
    }
    handle = dlopen(path.c_str(), RTLD_LAZY);
    if(!handle)
      return -1;
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

