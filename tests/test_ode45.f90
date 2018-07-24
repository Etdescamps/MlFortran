Program test_ode45
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  Use mlf_ode45
  Use test_functions
  call test()
Contains
  Subroutine test()
    type(mlf_odeTest) :: fun
    type(mlf_ode45_obj) :: ode
    integer :: i, N=100000000, info
    integer(c_int64_t) :: nstep =20
    call fun%init(FArenstorf)
    info = ode%init(fun, X0Arenstorf, tMax = TEndArenstorf, atoli = 1d-8, rtoli = 1d-8)
    if(info < 0) RETURN
    print *, 0, ode%t0, ode%X0
    Do i=1,N
      info = ode%step(niter = nstep)
      print *, i, ode%t, ode%X
      if(info /= 0) EXIT
    End Do
    
  End Subroutine test

End

