/*-
 * Copyright (c) 2017-2018 Etienne Descamps
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation and/or
 *    other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be
 *    used to endorse or promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "mlf.hpp"
#include <unistd.h>

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

  const char *MlfFinalizeCopy::what() const noexcept {
    return "MlfFinalizeCopy";
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
    int idRVar = MlfObject::getIdName("rVar");
    int idIVar = MlfObject::getIdName("iVar");
    MlfObject::getRsc(idRVar, rVar);
    MlfObject::getRsc(idIVar, iVar);
    is_initData = true;
  }

  void MlfStepObject::printLine(ostream& os) {
    if(!is_initData)
      initOutput();
    os << iVar << rVar << std::endl;
  }
  void MlfStepObject::printFields(ostream& os) {
    int idRVar = MlfObject::getIdName("rVar");
    int idIVar = MlfObject::getIdName("iVar");
    os << MlfObject::getFields(idIVar) << ';' <<MlfObject::getFields(idRVar) << std::endl;
  }

  bool MlfHdf5::readWOrCreate(string &fileName) {
    if(access(fileName.c_str(), F_OK) < 0) {
      if(createFile(fileName) < 0)
        throw MlfException(mlf_FILEERROR);
    }
    else if(openFile(fileName) < 0)
      throw MlfException(mlf_FILEERROR);
    return hasData();
  }

  void MlfOptimObject::printMinX(ostream& os) {
    if(!minX.associated()) {
      int id = MlfObject::getIdName("minX");
      MlfObject::getRsc(id, minX);
    }
    os << minX << std::endl;
  }

}

