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
      default:
        return "MlfRessource error unknown";
    }
  }
}

