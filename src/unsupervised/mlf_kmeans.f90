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

Module mlf_kmeans
  Use ieee_arithmetic
  Use iso_c_binding
  Use iso_fortran_env
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_step_algo
  Use mlf_models
  Use mlf_utils
  Use mlf_rand
  Use mlf_errors
  IMPLICIT NONE
  
  PRIVATE
  public :: mlf_kmeans_c, mlf_match_points

  ! Kmean object type
  Type, Public, extends(mlf_class_model) :: mlf_algo_kmeans
    real(c_double), pointer :: X(:,:), Mu(:,:), MDist
    real(c_double), allocatable :: minDist(:)
    integer(c_int32_t), allocatable :: cl(:)
    logical :: initialized
  Contains
    procedure :: init => mlf_algo_kmeans_init
    procedure :: reinit => mlf_algo_kmeans_reinit
    procedure :: stepF => mlf_algo_kmeans_stepF
    procedure :: getClass => mlf_algo_kmean_getClass
  End Type mlf_algo_kmeans

Contains
  ! Simple K-means algorithm implementation

  ! Init function for the step object
  integer Function mlf_algo_kmeans_init(this, X, nC, Mu, data_handler) result(info)
    class(mlf_algo_kmeans), intent(out), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), target :: X(:,:)
    real(c_double),  optional :: Mu(:,:)
    integer, intent(in), optional :: nC
    logical :: fixed_dims(2) = [.TRUE., .FALSE.]
    integer(c_int64_t) :: nipar, nrpar, nY, nX, nd(2)
    integer :: nrsc
    nipar = 0; nrpar = 1; nrsc = 1
    info = mlf_step_obj_init(this, nipar, nrpar, nrsc, C_CHAR_"", C_CHAR_"Mdist;", data_handler = data_handler)
    if(info < 0) RETURN
    nX = size(X,2); nY = size(X,1)
    nd(1) = nY
    if(present(Mu)) then
      nd = shape(Mu)
      if(nd(1) /= nY) then
        write (error_unit, *) "Dimensions of X and Mu don't match: ", nY, " =/ ", nd(1)
        info = -1
        return
      endif
      fixed_dims = .TRUE.
    else if(present(nC)) then
      nd = [nY, int(nC, kind=8)]
      fixed_dims = .TRUE.
    else if(.NOT. present(data_handler)) then
      info = -1; RETURN
    endif
    this%X => X
    info = this%add_rmatrix(nrsc, nd, this%Mu, C_CHAR_"Mu", data_handler = data_handler, &
      fixed_dims = fixed_dims)
    if(CheckF(info, "Error adding Mu")) RETURN
    if(present(data_handler)) then
      this%initialized = .TRUE.
    else if(present(Mu)) then
      info = this%reinit()
      this%Mu = Mu
      this%initialized = .TRUE.
    else
      info = this%reinit()
    endif
    if(CheckF(info, "Error reinit")) RETURN
    ALLOCATE(this%minDist(nX), this%cl(nX))
    this%MDist => this%rpar(nrpar)
  End Function mlf_algo_kmeans_init

  ! C interface to k-means algorithm
  type(c_ptr) Function mlf_kmeans_c(pX, nX, nY, nC, pMu) result(cptr) bind(C, name="mlf_kmeans_c")
    type(c_ptr), value :: pX, pMu
    real(c_double), pointer :: X(:,:), Mu(:,:)
    integer(c_int), value :: nX, nY, nC
    integer :: info
    type(mlf_algo_kmeans), pointer :: this
    class (mlf_obj), pointer :: obj
    ALLOCATE(this)
    call C_F_POINTER(pX, X, [nY, nX])
    if(C_ASSOCIATED(pMu)) then
      call C_F_POINTER(pMu, Mu, [nC, nX])
      info = this%init(X, Mu = Mu)
    else
      info = this%init(X, nC = nC)
    endif
    obj => this
    cptr = c_allocate(obj)
  End Function mlf_kmeans_c

  ! Reinit algorithm parameters
  integer Function mlf_algo_kmeans_reinit(this) result(info)
    class(mlf_algo_kmeans), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    if(associated(this%Mu)) this%Mu = 0
    this%initialized = .FALSE.
  End Function mlf_algo_kmeans_reinit

  ! Algorithm step function
  integer Function mlf_algo_kmeans_stepF(this, niter) result(info)
    class(mlf_algo_kmeans), intent(inout), target :: this
    integer(kind=8), intent(inout), optional :: niter
    integer(kind=8) :: i, N
    N = 1
    if(present(niter)) N = niter
    do i=1,N
      if(.NOT. this%initialized) then
        call InitKMeans(this%X, this%Mu)
        this%initialized = .TRUE.
      else
        call EvaluateClass(this%X, this%Mu, this%cl, this%minDist)
        this%MDist = mean(this%minDist)
        call EvaluateCentres(this%X, this%Mu, this%cl)
      endif
    end do
    info = 0
  End Function mlf_algo_kmeans_stepF

  integer Function mlf_algo_kmean_getClass(this, X, Cl) result(info)
    class(mlf_algo_kmeans), intent(in), target :: this
    real(c_double), intent(in) :: X(:,:)
    integer(c_int), intent(out) :: Cl(:)
    info = 0
    call EvaluateClass(X, this%Mu, Cl)
  End Function mlf_algo_kmean_getClass

  Subroutine EvaluateClass(X, Mu, cl, minDist)
    ! Naive algorithm (TODO use a better structure when the number of points is huge)
    real(c_double), intent(in) :: X(:,:), Mu(:,:)
    real(c_double), intent(out), optional :: minDist(:)
    integer, intent(out), optional :: cl(:)
    real(c_double) :: Dist(size(X,2),size(Mu,2))
    integer :: ND, NX, NC, i, j
    NX = size(X,2)
    ND = size(X,1)
    NC = size(Mu,2)
    Do i=1,NC
      Forall (j=1:NX) dist(j,i) = norm2(X(:,j)-Mu(:,i))
    End Do
    if(present(cl)) cl = minloc(dist, 2) 
    if(present(minDist)) minDist = minval(dist, 2)
  End Subroutine EvaluateClass

  subroutine EvaluateCentres(X, Mu, cl)
    real(c_double), intent(in) :: X(:,:)
    real(c_double), intent(out) :: Mu(:,:)
    integer, intent(in) :: cl(:)
    integer :: NX, NC, i, k, Ncl(size(cl,1))
    NX = size(X,2)
    NC = size(Mu,2)
    Mu = 0
    Ncl = 0
    Do i=1, NX
      k = cl(i)
      Ncl(k) = Ncl(k)+1
      Mu(:,k) = Mu(:,k)+X(:,i)
    End Do
    Do i=1,NC
      if(Ncl(i) == 0) then
        ! Set to the farest position
        Mu(:,i) = GetFarPoint(X,Mu)
      else
        Mu(:,i) = Mu(:,i)/real(Ncl(i), kind=c_double)
      endif
    End Do
  End Subroutine EvaluateCentres

  Function GetFarPoint(X,Mu) result(V)
    real(c_double), intent(in) :: X(:,:), Mu(:,:)
    real(c_double) :: V(size(X,1)), minDist(size(X,2))
    integer :: k(1)
    ! Use in this case to get the farest point from the selected centres Mu
    call EvaluateClass(X,Mu, minDist = minDist)
    call mlf_rand_class(mlf_cumulativefromvect(minDist*minDist), k)
    V = X(:,k(1))
  End Function GetFarPoint

  Subroutine InitKMeans(X, Mu)
    real(c_double), intent(in) :: X(:,:)
    real(c_double), intent(out) :: Mu(:,:)
    integer :: NX, NC, i
    NC = size(Mu,2)
    NX = size(X,2)
    ! Determine mean of the points
    Mu(:,1) = mean(X,2)
    ! Select the farthest point
    Mu(:,1) = GetFarPoint(X, Mu(:,1:1))
    Do i=2,NC
      ! Reiter the farthest point selection method
      Mu(:,i) = GetFarPoint(X, Mu(:,1:i-1))
    End Do
  End Subroutine InitKMeans

  Subroutine mlf_match_points(X, Y, idx)
    real(c_double), intent(in) :: X(:,:), Y(:,:)
    integer(c_int), intent(out) :: idx(:)
    real(c_double) :: minDist, d
    integer :: i, j, N, k
    N = size(X,2)
    idx = [(i, i=1,N)]
    Do i=1,N
      minDist = norm2(X(:,i)-Y(:,idx(i)))
      k = i
      Do j=i+1,N
        d = norm2(X(:,i)-Y(:,idx(j)))
        if(d < minDist) then
          d = minDist
          k = j
        endif
      End Do
      call XChange(idx(i), idx(k))
    End Do
  End Subroutine mlf_match_points
End Module mlf_kmeans

