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

Module mlf_pareto
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_cfuns
  Use mlf_utils
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_pareto_dominate, mlf_getPareto, mlf_crowndingDistance, mlf_paretoSelect

  Enum, bind(C)
    enumerator :: mlf_PARETO_RM_DUPLICATES = 0, mlf_PARETO_KEEP_DUPLICATES = 1,  mlf_PARETO_MOVE_PF_DUPLICATES = 2
  End Enum
  Public :: mlf_PARETO_RM_DUPLICATES, mlf_PARETO_KEEP_DUPLICATES, mlf_PARETO_MOVE_PF_DUPLICATES

Contains

  integer Function clear_ids(idx) result(i)
    integer(c_int), intent(inout) :: idx(:)
    integer :: k
    i = 0
    Do k=1,size(idx)
      if(idx(k)==0) CYCLE
      i = i+1
      idx(i) = idx(k)
    End Do
  End Function clear_ids

  integer Function clear_ids_Y(idx, Y) result(i)
    integer(c_int), intent(inout) :: idx(:)
    real(c_double), intent(inout) :: Y(:)
    integer :: k
    i = 0
    Do k=1,size(idx)
      if(idx(k)==0) CYCLE
      i = i+1
      Y(i) = Y(k)
      idx(i) = idx(k)
    End Do
  End Function clear_ids_Y


  ! Paretor dominance function for an indexed case with X before Y
  ! The objective are minimized in this case (use -X and -Y for maximization)
  integer Function mlf_pareto_dominate(X, Y) result(r)
    real(c_double), intent(in) :: X(:), Y(:)
    integer :: i, k, N
    N = size(X)
    r = 0
    do i=1,N
      ! Look at the first indice that unmatches
      if(X(i)<Y(i)) then
        r = 1; k=i
        EXIT; ! X dominate Y on i 
      else if(X(i)>Y(i)) then
        r = -1; k=i
        EXIT; ! X is dominated by Y on i 
      endif
    end do
    if(r > 0) then
      do i=k+1,N
        if(X(i)>Y(i)) then
          ! No dominance at all
          r = 0; RETURN
        endif
      end do
      ! X is dominated by Y
    else if(r < 0) then
      do i=k+1,N
        if(X(i)<Y(i)) then
          ! No dominance at all
          r = 0; RETURN
        endif
      end do
      ! X dominates Y
    else
      ! In this case: X == Y
      r = 2 ! X dominate Y in this indexed method
      ! (index of X before those of Y and X==Y => X dominate Y)
      ! In other case r == 2 shall be interpreted as a non dominance
    end if
  End Function mlf_pareto_dominate

  integer Function ParetoReduce(Y, idList, idPareto) result(fDom)
    real(c_double), intent(inout) :: Y(:,:)
    integer(c_int), intent(inout) :: idList(:)
    integer(c_int), intent(out) :: idPareto(:)
    integer :: lastDom, nDom, i, j, K
    nDom = 1; fDom = 1; lastDom =1; K = size(Y,2)
    idPareto(fDom) = idList(1)
    idList(1) = 0
    MainLoop: Do i=2,K
      if(idList(i)<0) then
        if(idList(i-1) == 0) then
          fDom = fDom+1
          idPareto(fDom) = idList(i)
          idList(i) = 0
        endif
        CYCLE MainLoop
      endif
      SecondLoop: Do j=nDom,1,-1
        if(j<lastDom) EXIT SecondLoop
        Select case(mlf_pareto_dominate(Y(:,j), Y(:,i)))
        case(0)
          CYCLE SecondLoop ! Y(:,i) is not dominated
        case(-1)
          ! Y(:,i) dominate, on index /= 1, the element j
          if(j>lastDom) Y(:,lastDom+1:j) = Y(:,lastDom:j-1)
          ! Faster than the evaluation of the rest of the Pareto front
          ! and remove one element to evaluate
          lastDom = lastDom+1
          ! If an element is not dominated by i, it won't be dominated by j,
          ! So j is removed from the tested Pareto set
          EXIT SecondLoop ! i is not Pareto dominated by the others elements
        ! case(0) -> continue evaluation of Pareto front
        case default
          CYCLE MainLoop ! Y(:,i) is dominated
        End Select
      End Do SecondLoop
      ! i is not Pareto dominated by previous elements
      ! so it is in Pareto front
      nDom = nDom + 1
      Y(:,nDom) = Y(:,i)
      fDom = fDom+1
      idPareto(fDom) = idList(i)
      idList(i) = 0
    End Do MainLoop
    ! return nDom
  End Function ParetoReduce

  Function ParetoFront2(Y, idList, idPareto) result(numP)
    integer(c_int), allocatable :: numP(:), lP(:)
    real(c_double), intent(inout) :: Y(:)
    integer(c_int), intent(inout) :: idList(:)
    integer(c_int), intent(out) :: idPareto(:)
    integer :: i, j, K, nRemain, nPareto, M
    real(c_double) :: lastMin
    K = size(Y)
    nRemain = K; nPareto = 0; M = K
    ALLOCATE(lP(K))
    D1: Do i=1,K
      lastMin = HUGE(lastMin)
      D2: Do j=1,M
        if(idList(j)==0) CYCLE D2
        if(Y(j) > lastMin) CYCLE D2
        if(Y(j) == lastMin) then
          if(idList(j) < 0 .AND. idList(j-1) /= 0) CYCLE D2
          ! Equals elements are always in the same Pareto front and placed consecutively
        endif
        nPareto = nPareto+1
        idPareto(nPareto) = idList(j)
        idList(j) = 0
        nRemain = nRemain-1
        if(nRemain == 0) EXIT D2
        lastMin = Y(j)
      End Do D2
      lP(i) = nPareto
      if(nRemain == 0) then
        numP = lP(:i)
        RETURN
      endif
      if(nRemain > 16 .AND. nRemain*2 < M) then
        M = clear_ids_Y(idList, Y)
      endif
    End Do D1
    numP = lP ! Corner case when each point has its paretor front
  End Function ParetoFront2

  Function mlf_getPareto(V, idPareto, mu0, equalcase0) result(numP)
    real(c_double), intent(in), contiguous :: V(:,:)
    integer(c_int), intent(out) :: idPareto(:)
    real(c_double), allocatable :: Y(:,:)
    integer(c_int), allocatable :: idSorted(:), numP(:)
    integer, intent(in), optional :: mu0, equalcase0
    integer(c_int) :: N, ND, K, i, j, M, NumPareto, mu, equalcase
    ND = size(V,1); N = size(V,2)
    ALLOCATE(idSorted(N))
    equalcase = mlf_PARETO_RM_DUPLICATES
    if(present(equalcase0)) equalcase = equalcase0
    K = N
    select case (equalcase)
    ! Use QuickSort that uses lexicographical order on the objectives
    case(mlf_PARETO_RM_DUPLICATES)
      ! Variant that removes the duplicates entries
      K = c_qsort_unify(V, idSorted, N, ND, ND, N)
    case(mlf_PARETO_KEEP_DUPLICATES)
      ! Variant that keep duplicates entries on the same pareto front
      CALL c_qsort_neg(V, idSorted, N, ND, ND, N)
    case(mlf_PARETO_MOVE_PF_DUPLICATES)
      ! Variant that keep one duplicate on the main Pareto front (those with the lowest id)
      ! and move to other Pareto fronts the other ones
      CALL c_qsort(V, idSorted, N, ND, ND, N)
    end select
    mu = K
    if(present(mu0)) then
      if(mu0<K) mu = mu0
    endif
    ! Lexicographic order so: if i<j, j cannot Pareto dominate i
    ALLOCATE(Y(ND-1,K))
    j = 1; M = K
    numPareto = -1 ! Remove warning
    if(ND == 2) then ! Easier case
      Y(:,:M) = V(2:, idSorted(:M))
      numP = ParetoFront2(Y(1,:M), idSorted(:M), idPareto)
      RETURN
    endif
    Do i=1,K
      if(M <= 1) then
        idPareto(j) = idSorted(1)
        NumPareto = i ! Don't evaluate the Pareto front if the number of element id 1
        idSorted(N-i+1) = j
        EXIT
      endif
      Y(:,:M) = V(2:, idSorted(:M))
      j = j + ParetoReduce(Y(:,:M), idSorted(:M), idPareto(j:))
      M = clear_ids(idSorted(:M))
      ! clear_ids remove at least one element (so idSorted(N-i+1) is free)
      idSorted(N-i+1) = j-1
      if(j >= mu) then
        NumPareto = i ! Normal case
        EXIT
      endif
    End Do
    ALLOCATE(numP(numPareto))
    forall(i=1:numPareto) numP(i) = idSorted(N-i+1)
  End Function mlf_getPareto

  Subroutine mlf_crowndingDistance(V, cDist)
    real(c_double), intent(in) :: V(:,:)
    real(c_double), intent(out) :: cDist(:)
    real(c_double), allocatable :: Y(:)
    integer(c_int), allocatable :: idSorted(:)
    integer(c_int) :: N, ND, i, j 
    real(c_double) :: dx
    ND = size(V,1); N = size(V,2)
    ALLOCATE(idSorted(N), Y(N))
    cDist = 0
    Do i=1,ND
      Y = V(i,:)
      ! QuickSort using lexicographical order on the objectives
      CALL QSortIdx(Y, idSorted)
      cDist(idSorted(1))=HUGE(dx)
      cDist(idSorted(N))=HUGE(dx)
      dx = 1d0/(Y(idSorted(1))+Y(idSorted(N)))
      Forall(j=2:N-1) cDist(idSorted(j)) = cDist(idSorted(j))+dx*(Y(idSorted(j+1))-Y(idSorted(j-1)))
    End Do
  End Subroutine mlf_crowndingDistance

  Subroutine mlf_paretoSelect(V, mu, idx, idP)
    real(c_double), intent(in), contiguous :: V(:,:)
    integer(c_int), intent(inout) :: idx(:), idP(:)
    integer, intent(in) :: mu
    integer(c_int), allocatable :: idC(:), idD(:)
    real(c_double), allocatable :: cDist(:)
    integer :: i, nP, k0, k1
    nP = size(idP)
    Do i=1,nP
      if(idP(i) == mu) RETURN
      if(idP(i) > mu) then
        k0 = 1
        if(i>1) k0 = idP(i-1)+1
        k1 = idP(i)
        ALLOCATE(cDist(k1-k0+1), idD(mu-k0+1), idC(k1-k0+1))
        idC = idx(k0:k1)
        call mlf_crowndingDistance(V(:,idC), cDist)
        call QSortIdx(cDist, idD, mu-k0+1)
        idx(k0:mu) = idC(idD)
        RETURN
      endif
    End Do
  End Subroutine mlf_paretoSelect

End Module mlf_pareto

