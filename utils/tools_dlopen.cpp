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
      case DlErrorType::InvalidFunctionType:
        return "DlException: Invalid function type";
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

  void LibraryFun::init(string path, string funPrefix, LibraryFunType typeFun, string fileName) {
    DlLoader::init(path);
    mlf_init_fun finit = DlLoader::getSym<mlf_init_fun>(funPrefix+"_init");
    ffree = DlLoader::getSym<mlf_free_fun>(funPrefix+"_free");
    mlf_getinfo_fun finfo = DlLoader::getSym<mlf_getinfo_fun>(funPrefix+"_getinfo");
    data = finit(fileName.c_str());
    for(int i=mlf_NAME; i <= mlf_FIELDS; i++)
      description[i] = finfo(data, i);
    switch(typeFun) {
      case LibraryFunType::OptimFun:
        {
          mlf_objective_fun fobj = DlLoader::getSym<mlf_objective_fun>(funPrefix+"_objfun");
          mlf_objective_fun fcstr = DlLoader::getSymOrNull<mlf_objective_fun>(funPrefix+"_cstrfun");
          object = mlf_objfunction(fobj, data, fcstr);
        }
        break;
      case LibraryFunType::BasisFun:
        {
          mlf_basis_fun fbasis = DlLoader::getSym<mlf_basis_fun>(funPrefix+"_basisfun");
          object = mlf_basisfunction(fbasis, data);
        }
        break;
      default:
        throw DlException(DlErrorType::InvalidFunctionType);
    }
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

