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
  Use mlf_hdf5
  Use test_functions
  integer :: infx
  integer, parameter :: Npoints = 1000000
  infx = mlf_init()
  call test()
  infx = mlf_quit()
Contains
  Subroutine test()
    type(mlf_odeTestCstr) :: fun
    !type(mlf_odeTestCstrIds) :: fun
    type(mlf_ode45_obj) :: ode
    type(mlf_hdf5_file) :: h5f
    integer :: i, N=1000000, info, j, k
    real(c_double), allocatable :: trajectory(:,:), steps(:,:)
    real(c_double) :: t
    ALLOCATE(trajectory(4, Npoints), steps(4, N))
    info = h5f%createFile("arenstorf.h5")
    CALL fun%init(FArenstorf, RESHAPE([1d0,0d0,0d0,0d0], [4,1]))
    !call fun%init(FArenstorf, [1])
    info = ode%init(fun, X0Arenstorf, tMax = TEndArenstorf, atoli = 1d-6, rtoli = 1d-6)
    If(info < 0) RETURN
    j = 1
    PRINT *, ode%t0, ode%X0
    Do i=1,N
      info = ode%step()
      If(info == mlf_ODE_SoftCstr) Then
        PRINT *, ode%t, ode%X
        !ode%X(1) = 0
      Endif
      Do k = j,Npoints
        t = real(k-1,8)/real(Npoints-1,8)*ode%tMax
        If(t>ode%t) EXIT
        call ode%denseEvaluation(t, trajectory(:,k))
      End Do
      steps(:,i) = ode%X
      j = k
      If(info == mlf_ODE_StopTime) EXIT
    End Do
    PRINT *, ode%t, ode%X
    Do k = j,Npoints
      t = real(k-1,8)/real(Npoints-1,8)*ode%tMax
      call ode%denseEvaluation(t, trajectory(:,k))
    End Do
    PRINT *, "NSTEPS = ", i
    info = h5f%pushData(trajectory, 'trajectory')
    info = h5f%pushData(steps(:, 1:i), 'steps')
    call h5f%finalize()
  End Subroutine test

End

