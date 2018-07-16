! Copyright (c) 2017-2018 Etienne Descamps
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation and/or
!    other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


! Contains is an interface to all optimization methods

Module mlf_optim_intf
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_cfuns
  Use mlf_rsc_array
  Use mlf_step_algo
  Use mlf_optim
  Use mlf_fun_intf
  Use mlf_matrix
  Use mlf_rand
  Use mlf_utils
  Use mlf_errors
  Use mlf_cmaes
  IMPLICIT NONE
  PRIVATE

  Public :: c_getoptimobj
Contains
  type(c_ptr) Function c_getoptimobj(calgoname, cfunobj, chandler, ftarget, lambda, &
      mu, sigma) result(cptr) bind(C, name="mlf_getoptimobj")
    type(c_ptr), value :: calgoname, cfunobj, chandler
    type(mlf_optim_param) :: p
    integer(c_int), value :: lambda, mu
    real(c_double), value :: sigma, ftarget
    class(mlf_objective_fun), pointer :: funobj
    class(*), pointer :: obj 
    character(len=:, kind=c_char), allocatable, target :: algoname
    class(mlf_data_handler), pointer :: data_handler => NULL()
    call mlf_stringFromC(calgoname, algoname)
    obj => mlf_getobjfromc(cfunobj)
    if(associated(obj)) then
      select type (obj)
      class is (mlf_objective_fun)
        funobj => obj
      end select
    endif
    obj => mlf_getobjfromc(chandler)
    if(associated(obj)) then
      select type (obj)
      class is (mlf_data_handler)
        data_handler => obj
      end select
    endif
    obj => NULL()
    call p%init(funobj, lambda, mu, sigma = sigma, targetFun = ftarget)
    select case (algoname)
      case("maes") ! MA-ES
        obj => mlf_cmaes_objcreate(.TRUE., data_handler, p)
      case default ! CMA-ES
        obj => mlf_cmaes_objcreate(.FALSE., data_handler, p)
    end select
    cptr = c_allocate(obj)
  End Function c_getoptimobj


End Module mlf_optim_intf


