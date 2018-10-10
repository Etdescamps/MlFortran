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
  Public :: mlf_match_points, mlf_GetFarPoint, mlf_EvaluateClass, mlf_InitKMeans

  ! Kmean object type
  Type, Public, Abstract, extends(mlf_step_obj) :: mlf_algo_kmeans
    real(c_double), pointer :: X(:,:), Mu(:,:), meanDist
    real(c_double), allocatable :: minDist(:)
    integer(c_int32_t), allocatable :: cl(:)
    logical :: initialized
  Contains
    procedure :: init => mlf_algo_kmeans_init
    procedure :: reinit => mlf_algo_kmeans_reinit
  End Type mlf_algo_kmeans

  Type, Public, extends(mlf_class_model) :: mlf_model_kmeans
    class(mlf_algo_kmeans), pointer :: top
  Contains
    procedure :: getClass => mlf_model_kmean_getClass
  End Type mlf_model_kmeans


Contains
  ! Model initilisator
  Subroutine mlf_model_kmean_init(this, top)
    class(mlf_model), intent(out), allocatable :: this
    class(mlf_algo_kmeans), intent(in), target :: top
    ALLOCATE(mlf_model_kmeans :: this)
    select type(this)
      class is (mlf_model_kmeans)
        this%top => top
    end select
  End Subroutine mlf_model_kmean_init

  integer Function mlf_model_kmean_getClass(this, X, Cl) result(info)
    class(mlf_model_kmeans), intent(in), target :: this
    real(c_double), intent(in) :: X(:,:)
    integer(c_int), intent(out) :: Cl(:)
    info = 0
    call mlf_EvaluateClass(X, this%top%Mu, Cl)
  End Function mlf_model_kmean_getClass

  ! Init function for the step object
  integer Function mlf_algo_kmeans_init(this, X, nC, Mu, data_handler) result(info)
    class(mlf_algo_kmeans), intent(out), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), target :: X(:,:)
    real(c_double),  optional :: Mu(:,:)
    integer, intent(in), optional :: nC
    type(mlf_step_numFields) :: numFields
    logical :: fixed_dims(2)
    integer(c_int64_t) :: nY, nX, nd(2)
    fixed_dims = [.TRUE., .FALSE.]
    call numFields%initFields(nRVar = 1, nRsc = 1)
    info = mlf_step_obj_init(this, numFields, data_handler = data_handler)
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
    info = this%add_rmatrix(numFields, nd, this%Mu, C_CHAR_"Mu", data_handler = data_handler, &
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
    call this%addRVar(numFields, this%meanDist, "meanDist")
    if(.NOT. ALLOCATED(this%model)) call mlf_model_kmean_init(this%model, this)
  End Function mlf_algo_kmeans_init

  ! Reinit algorithm parameters
  integer Function mlf_algo_kmeans_reinit(this) result(info)
    class(mlf_algo_kmeans), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    if(associated(this%Mu)) this%Mu = 0
    this%initialized = .FALSE.
  End Function mlf_algo_kmeans_reinit

  Function mlf_GetFarPoint(X,Mu) result(V)
    real(c_double), intent(in) :: X(:,:), Mu(:,:)
    real(c_double) :: V(size(X,1)) ! Output of the function
    real(c_double) :: minDist(size(X,2))
    integer :: k(1)
    ! Use in this case to get the farest point from the selected centres Mu
    call mlf_EvaluateClass(X,Mu, minDist = minDist)
    call mlf_rand_class(mlf_cumulativefromvect(minDist*minDist), k)
    V = X(:,k(1))
  End Function mlf_GetFarPoint

  Subroutine mlf_EvaluateClass(X, Mu, cl, minDist)
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
  End Subroutine mlf_EvaluateClass


  Subroutine mlf_InitKMeans(X, Mu)
    real(c_double), intent(in) :: X(:,:)
    real(c_double), intent(out) :: Mu(:,:)
    integer :: NX, NC, i
    NC = size(Mu,2)
    NX = size(X,2)
    ! Determine mean of the points
    Mu(:,1) = mean(X,2)
    ! Select the farthest point
    Mu(:,1) = mlf_GetFarPoint(X, Mu(:,1:1))
    Do i=2,NC
      ! Reiter the farthest point selection method
      Mu(:,i) = mlf_GetFarPoint(X, Mu(:,1:i-1))
    End Do
  End Subroutine mlf_InitKMeans

  Subroutine mlf_match_points(X, Y, idx)
    real(c_double), intent(in) :: X(:,:), Y(:,:)
    integer(c_int), intent(out) :: idx(:)
    real(c_double) :: dist(size(X,2), size(X,2))
    integer :: i, j, N, k, pos(2)
    N = size(X,2)
    ! Compute all pair of distance between the points of X and Y
    forall(i=1:N, j=1:N) dist(i,j) = norm2(X(:,i)-Y(:,j))
    Do k=1,N
      ! Find the best match
      pos = MINLOC(dist)
      i = pos(1); j = pos(2)
      ! Remove i and j from possible solutions
      dist(i,:) = HUGE(1d0)
      dist(:,j) = HUGE(1d0)
      idx(i) = j ! Associate j with i
    End Do
  End Subroutine mlf_match_points
End Module mlf_kmeans

