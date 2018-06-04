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

Program test_cmaes
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_utils
  Use mlf_matrix
  Use mlf_fun_intf
  Use mlf_cmaes
  Use test_functions
  integer :: lambda=120, ND = 30, Nmax = 2000, Niter = 30, NFirst
  integer :: larg = 32, info
  real(c_double) :: ftarget = 1.0e-10
  real(c_double), allocatable :: X0(:)
  character(len=32) :: arg
  type (mlf_maes_obj), target :: maesOpt
  type (mlf_cmaes_obj), target :: cmaesOpt
  class (mlf_matrix_es_obj), pointer :: optObj => maesOpt
  type (mlf_objective_test) :: fTest
  procedure (mlf_obj_test_fun), pointer :: optFun
  optFun => OFrosenbrock
  if(COMMAND_ARGUMENT_COUNT()>0) then
    call GET_COMMAND_ARGUMENT(1, arg, larg)
    select case (arg)
    case('sphere')
      optFun => OFsphere
    case('cigar')
      optFun => OFcigar
    case('tablet')
      optFun => OFtablet
    case('ellipsoid')
      optFun => OFellipsoid
    case('rosenbrock')
      optFun => OFrosenbrock
    case default
      print *, "Error invalid function name"
      optFun => OFrosenbrock
    end select
  endif 

  if(COMMAND_ARGUMENT_COUNT()>1) then
    larg = 32
    call GET_COMMAND_ARGUMENT(2, arg, larg)
    select case (arg)
    case('cmaes')
      optObj => cmaesOpt
    case('maes')
      optObj => maesOpt
    case default
      print *, "Error invalid method name"
      optObj => cmaesOpt
    end select
  else
    optObj => cmaesOpt
  endif
   if(COMMAND_ARGUMENT_COUNT()>2) then
    larg = 32
    call GET_COMMAND_ARGUMENT(3, arg, larg)
    read(arg,*,iostat=info) ND
  endif
  if(COMMAND_ARGUMENT_COUNT()>3) then
    larg = 32
    call GET_COMMAND_ARGUMENT(4, arg, larg)
    read(arg,*,iostat=info) lambda
  endif
  if(COMMAND_ARGUMENT_COUNT()>4) then
    larg = 32
    call GET_COMMAND_ARGUMENT(5, arg, larg)
    read(arg,*,iostat=info) Nmax
  endif
  if(COMMAND_ARGUMENT_COUNT()>5) then
    larg = 32
    call GET_COMMAND_ARGUMENT(6, arg, larg)
    read(arg,*,iostat=info) ftarget
  endif
  if(COMMAND_ARGUMENT_COUNT()>6) then
    larg = 32
    call GET_COMMAND_ARGUMENT(7, arg, larg)
    read(arg,*,iostat=info) Niter
  endif
  call fTest%init(optFun, ND)
  ALLOCATE(X0(ND))
  X0 = 1d0 
  Select Type(optObj)
  class is (mlf_maes_obj)
    info = optObj%init(fTest, X0 = X0, lambdaIn = lambda)
  class is (mlf_cmaes_obj)
    info = optObj%init(fTest, X0 = X0, lambdaIn = lambda, covEvery = 1)
  End Select
  DEALLOCATE(X0)
  call optObj%constraints_step(niterMax = int(nMax, kind=8))
  call optObj%constraints_optim(targetFun = ftarget)
  Do While (optObj%step() == 0)
    print *, optObj%niter, optObj%cputime, optObj%realtime, optObj%minFun, optObj%sigma
  End Do
  print *, optObj%minX
  Select Type(optObj)
  class is (mlf_maes_obj)
    print *, "Algorithm MA-ES"
  class is (mlf_cmaes_obj)
    print *, "Algorithm CMA-ES"
  End Select
  if(COMMAND_ARGUMENT_COUNT()>0) then
    call GET_COMMAND_ARGUMENT(1, arg, larg)
    print *, "Function: ", arg
  else
    print *, "Function: rosenbrock"
  endif
  print *, "NDim: ", ND, "lambda: ", lambda, "mu: ", lambda/2
End

