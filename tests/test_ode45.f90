Program test_ode45
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  Use mlf_ode_class
  Use mlf_ode45
  Use test_functions
  call test()
Contains
  Subroutine test()
    type(mlf_odeTest) :: fun
    type(mlf_ode45_obj) :: ode
    integer :: i, N=100000000, info
    integer(c_int64_t) :: nstep
    call fun%init(FArenstorf)
    info = ode%init(nX = 4_8, atoli = 1d-6, rtoli = 1d-6)
    If(info < 0) RETURN
    CALL ode%setFun(fun)
    info = ode%initODE(X0Arenstorf, tMax = TEndArenstorf)
    If(info < 0) RETURN
    PRINT *, 0, ode%t0, ode%X0
    Do i=1,N
      nstep = 2 
      info = ode%step(niter = nstep)
      PRINT *, ode%nFun, ode%t, ode%X, int(nstep,4)
      If(info == mlf_ODE_StopTime) EXIT
    End Do
    
  End Subroutine test

End

