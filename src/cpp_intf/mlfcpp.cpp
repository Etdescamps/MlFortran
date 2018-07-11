#include "mlf.hpp"

namespace MlFortran {
  const char *MlfException::what() const noexcept {
    switch(errorType) {
      case mlf_OK:
        return "MlfException: Error 0 (shall not happen)";
      case mlf_UNINIT:
        return "MlfException: Uninitialized element";
      case mlf_FUNERROR:
        return "MlfException: Error with function type";
      case mlf_WRONGTYPE:
        return "MlfException: Wrong type of object";
      case mlf_WRONGRANK:
        return "MlfException: Wrong rank";
      case mlf_FILENOTFOUND:
        return "MlfException: File not found";
      case mlf_FILEERROR:
        return "MlfException: Cannot access file";
      case mlf_OTHERERROR:
        return "MlfException: Other error type (shall use a different exception class?)";
      default:
        return "MlfException: Unknown error (shall not happen)";
    }
  }

  const char *MlfOutOfBounds::what() const noexcept {
    return "MlfOutOfBounds";
  }

  const char *MlfRessourceError::what() const noexcept {
    switch(rscError) {
      case MlfRscErrorType::NotWritable:
        return "MlfRessource not writable (shall use only a readonly structure)";
      case MlfRscErrorType::NotReadable:
        return "MlfRessource not readable (shall use a writable structure)";
      case MlfRscErrorType::InvalidAccessType:
        return "MlfRessource invalid MLF_ACCESSTYPE";
      case MlfRscErrorType::NotAllocated:
        return "MlfRessource not allocated";
      case MlfRscErrorType::NotFound:
        return "MlfRessource not found";
      default:
        return "MlfRessource error unknown";
    }
  }

  bool MlfObject::updateIdMap() {
    nrsc = getNumRsc();
    obj_names = std::make_unique<string_view[]>(nrsc+1);
    for(int i=1; i <= nrsc; ++i) {
      const char *name = getName(i);
      obj_names[i] = name ? name : "";
    }
    return true;
  }

  int MlfObject::getIdName(const string &name, bool updated) {
    if(!obj_names)
      updated = updateIdMap();
    for(int i=1; i <= nrsc; ++i) {
      if(obj_names[i] == name) {
        if(string_view(getName(i)) != name) {
          if(updated)
            throw MlfRscErrorType(MlfRscErrorType::NotFound);
          return getIdName(name, updateIdMap());
        }
        return i;
      }
    }
    if(updated)
      throw MlfRscErrorType(MlfRscErrorType::NotFound);
    return getIdName(name, updateIdMap());
  }

  void MlfStepObject::initOutput() {
    int idrpar = MlfObject::getIdName("rpar");
    int idipar = MlfObject::getIdName("ipar");
    MlfObject::getRsc(idrpar, rdata);
    MlfObject::getRsc(idipar, idata);
    is_initData = true;
  }

  void MlfStepObject::printLine(ostream& os) {
    if(!is_initData)
      initOutput();
    os << idata << rdata << std::endl;
  }
  void MlfStepObject::printFields(ostream& os) {
    int idrpar = MlfObject::getIdName("rpar");
    int idipar = MlfObject::getIdName("ipar");
    os << MlfObject::getFields(idipar) << MlfObject::getFields(idrpar) << std::endl;
  }

  bool MlfHdf5::readWOrCreate(string &fileName) {
    if(openFile(fileName) < 0) {
      if(createFile(fileName) < 0)
        throw MlfException(mlf_FILEERROR);
    }
    return hasData();
  }

}

