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

! TODO: find a better way for handling multiple dimension and data type
Module mlf_rsc_array
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_errors
  Use mlf_cfuns
  IMPLICIT NONE
  PRIVATE
  Public :: mlf_arr_init

  Type, Public :: mlf_rsc_numFields
    integer(c_int64_t) :: nIPar = 0
    integer(c_int64_t) :: nRPar = 0
    integer(c_int) :: nRsc = 0
  Contains
    procedure :: incNRsc => mlf_rsc_incNRsc
  End Type mlf_rsc_numFields

  ! Integer (64) ressource handler
  Type, Public, extends(mlf_rsc_intf) :: mlf_rsc_int64_1d
    integer(c_int64_t), allocatable :: V(:)
  Contains
    procedure :: getd => mlf_int64_1d_getd
    procedure :: updated => mlf_int64_1d_updated
  End Type mlf_rsc_int64_1d

  ! Integer (32) ressource handler
  Type, Public, extends(mlf_rsc_intf) :: mlf_rsc_int32_1d
    integer(c_int32_t), allocatable :: V(:)
  Contains
    procedure :: getd => mlf_int32_1d_getd
    procedure :: updated => mlf_int32_1d_updated
  End Type mlf_rsc_int32_1d

  ! Double ressource handler
  Type, Public, extends(mlf_rsc_intf) :: mlf_rsc_double1d
    real(c_double), allocatable :: V(:)
  Contains
    procedure :: getd => mlf_double1d_getd
    procedure :: updated => mlf_double1d_updated
  End Type mlf_rsc_double1d

  ! Double ressource handler
  Type, Public, extends(mlf_rsc_intf) :: mlf_rsc_double2d
    real(c_double), allocatable :: V(:,:)
  Contains
    procedure :: getd => mlf_double2d_getd
    procedure :: updated => mlf_double2d_updated
  End Type mlf_rsc_double2d

  ! Double ressource handler
  Type, Public, extends(mlf_rsc_intf) :: mlf_rsc_double3d
    real(c_double), allocatable :: V(:,:,:)
  Contains
    procedure :: getd => mlf_double3d_getd
    procedure :: updated => mlf_double3d_updated
  End Type mlf_rsc_double3d

  ! Double ressource handler
  Type, Public, extends(mlf_rsc_intf) :: mlf_rsc_float1d
    real(c_float), allocatable :: V(:)
  Contains
    procedure :: getd => mlf_float1d_getd
    procedure :: updated => mlf_float1d_updated
  End Type mlf_rsc_float1d

  ! Double ressource handler
  Type, Public, extends(mlf_rsc_intf) :: mlf_rsc_float2d
    real(c_float), allocatable :: V(:,:)
  Contains
    procedure :: getd => mlf_float2d_getd
    procedure :: updated => mlf_float2d_updated
  End Type mlf_rsc_float2d

  ! Double ressource handler
  Type, Public, extends(mlf_rsc_intf) :: mlf_rsc_float3d
    real(c_float), allocatable :: V(:,:,:)
  Contains
    procedure :: getd => mlf_float3d_getd
    procedure :: updated => mlf_float3d_updated
  End Type mlf_rsc_float3d

  
  ! Generic object handler containing main parameters
  Type, Public, extends(mlf_obj) :: mlf_arr_obj
    integer(c_int64_t), pointer :: iPar(:)
    real(c_double), pointer :: rPar(:)
  Contains
    procedure :: print_rpar => mlf_print_rpar
    procedure :: print_ipar => mlf_print_ipar
    procedure :: add_i64array => mlf_obj_addI64array
    procedure :: add_i32array => mlf_obj_addI32array
    procedure :: add_rarray => mlf_obj_addRArray
    procedure :: add_rmatrix => mlf_obj_addRMatrix
    procedure :: add_rmatrix3D => mlf_obj_addRMatrix3D
    procedure :: add_farray => mlf_obj_addFArray
    procedure :: add_fmatrix => mlf_obj_addFMatrix
    procedure :: add_fmatrix3D => mlf_obj_addFMatrix3D
    procedure :: addRPar => mlf_obj_addRPar
    procedure :: addIPar => mlf_obj_addIPar
  End Type mlf_arr_obj

  ! Data handler for the HDF5 data import
  Type, Public, extends(mlf_rsc_handler), abstract :: mlf_data_handler
  Contains
    procedure (mlf_data_handler_getdata_double1d), deferred :: getdata_double1d
    procedure (mlf_data_handler_getdata_double2d), deferred :: getdata_double2d
    procedure (mlf_data_handler_getdata_double3d), deferred :: getdata_double3d
    procedure (mlf_data_handler_getdata_float1d), deferred :: getdata_float1d
    procedure (mlf_data_handler_getdata_float2d), deferred :: getdata_float2d
    procedure (mlf_data_handler_getdata_float3d), deferred :: getdata_float3d
    procedure (mlf_data_handler_getdata_int64), deferred :: getdata_int64
    procedure (mlf_data_handler_getdata_int32), deferred :: getdata_int32
    generic :: getData => getdata_double1d, getdata_double2d, getdata_double3d, &
      getdata_int64, getdata_int32, getdata_float1d, getdata_float2d, getdata_float3d
    procedure (mlf_date_handler_getSubobject), deferred :: getSubObject
  End Type mlf_data_handler

  Abstract Interface
    Function mlf_date_handler_getSubobject(this, obj_name) Result(ptr)
      Use iso_c_binding
      import :: mlf_data_handler
      class(mlf_data_handler), intent(inout), target :: this
      character(len=*,kind=c_char), intent(in) :: obj_name
      class(mlf_data_handler), pointer :: ptr
    End Function mlf_date_handler_getSubobject

    ! Get data functions
    Integer Function mlf_data_handler_getdata_double1d(this, rsc_data, rsc_name, dims, fixed_dims)
      Use iso_c_binding
      import :: mlf_data_handler
      class(mlf_data_handler), intent(inout), target :: this
      real(c_double), intent(out), target, allocatable :: rsc_data(:)
      character(len=*,kind=c_char), intent(in) :: rsc_name
      integer(c_int64_t), intent(in), optional :: dims(:)
      logical, intent(in), optional :: fixed_dims(:)
    End Function mlf_data_handler_getdata_double1d

    Integer Function mlf_data_handler_getdata_double2d(this, rsc_data, rsc_name, dims, fixed_dims)
      Use iso_c_binding
      import :: mlf_data_handler
      class(mlf_data_handler), intent(inout), target :: this
      real(c_double), intent(out), target, allocatable :: rsc_data(:,:)
      character(len=*,kind=c_char), intent(in) :: rsc_name
      integer(c_int64_t), intent(in), optional :: dims(:)
      logical, intent(in), optional :: fixed_dims(:)
    End Function mlf_data_handler_getdata_double2d

    Integer Function mlf_data_handler_getdata_double3d(this, rsc_data, rsc_name, dims, fixed_dims)
      Use iso_c_binding
      import :: mlf_data_handler
      class(mlf_data_handler), intent(inout), target :: this
      real(c_double), intent(out), target, allocatable :: rsc_data(:,:,:)
      character(len=*,kind=c_char), intent(in) :: rsc_name
      integer(c_int64_t), intent(in), optional :: dims(:)
      logical, intent(in), optional :: fixed_dims(:)
    End Function mlf_data_handler_getdata_double3d

    Integer Function mlf_data_handler_getdata_float1d(this, rsc_data, rsc_name, dims, fixed_dims)
      Use iso_c_binding
      import :: mlf_data_handler
      class(mlf_data_handler), intent(inout), target :: this
      real(c_float), intent(out), target, allocatable :: rsc_data(:)
      character(len=*,kind=c_char), intent(in) :: rsc_name
      integer(c_int64_t), intent(in), optional :: dims(:)
      logical, intent(in), optional :: fixed_dims(:)
    End Function mlf_data_handler_getdata_float1d

    Integer Function mlf_data_handler_getdata_float2d(this, rsc_data, rsc_name, dims, fixed_dims)
      Use iso_c_binding
      import :: mlf_data_handler
      class(mlf_data_handler), intent(inout), target :: this
      real(c_float), intent(out), target, allocatable :: rsc_data(:,:)
      character(len=*,kind=c_char), intent(in) :: rsc_name
      integer(c_int64_t), intent(in), optional :: dims(:)
      logical, intent(in), optional :: fixed_dims(:)
    End Function mlf_data_handler_getdata_float2d

    Integer Function mlf_data_handler_getdata_float3d(this, rsc_data, rsc_name, dims, fixed_dims)
      Use iso_c_binding
      import :: mlf_data_handler
      class(mlf_data_handler), intent(inout), target :: this
      real(c_float), intent(out), target, allocatable :: rsc_data(:,:,:)
      character(len=*,kind=c_char), intent(in) :: rsc_name
      integer(c_int64_t), intent(in), optional :: dims(:)
      logical, intent(in), optional :: fixed_dims(:)
    End Function mlf_data_handler_getdata_float3d

    Integer Function mlf_data_handler_getdata_int64(this, rsc_data, rsc_name, dims, fixed_dims)
      Use iso_c_binding
      import :: mlf_data_handler
      class(mlf_data_handler), intent(inout), target :: this
      integer(c_int64_t), intent(out), target, allocatable :: rsc_data(:)
      character(len=*,kind=c_char), intent(in) :: rsc_name
      integer(c_int64_t), intent(in), optional :: dims(:)
      logical, intent(in), optional :: fixed_dims(:)
    End Function mlf_data_handler_getdata_int64
    
    Integer Function mlf_data_handler_getdata_int32(this, rsc_data, rsc_name, dims, fixed_dims)
      Use iso_c_binding
      import :: mlf_data_handler
      class(mlf_data_handler), intent(inout), target :: this
      integer(c_int32_t), intent(out), target, allocatable :: rsc_data(:)
      character(len=*,kind=c_char), intent(in) :: rsc_name
      integer(c_int64_t), intent(in), optional :: dims(:)
      logical, intent(in), optional :: fixed_dims(:)
    End Function mlf_data_handler_getdata_int32
  End Interface

Contains

  Subroutine mlf_print_rpar(this)
    class(mlf_arr_obj), intent(in), target :: this
    print *, this%v(2)%r_fields
    print *, this%rpar
  End Subroutine mlf_print_rpar
  Subroutine mlf_print_ipar(this)
    class(mlf_arr_obj), intent(in), target :: this
    print *, this%v(1)%r_fields
    print *, this%ipar
  End Subroutine mlf_print_ipar

  Integer Function mlf_rsc_incNRsc(this) result(id)
    class(mlf_rsc_numFields), intent(inout) :: this
    id = this%nRsc+1
    this%nRsc = id
  End Function mlf_rsc_incNRsc

  Integer Function mlf_obj_addI64array(this, numFields, nd, pnt, rsc_name, rsc_fields, &
      data_handler, fixed_dims) result(info)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    class(mlf_rsc_numFields), intent(inout) :: numFields
    integer(c_int64_t), intent(out), pointer, optional :: pnt(:)
    integer(c_int64_t), intent (inout) :: nd
    character(len=*,kind=c_char), optional :: rsc_name, rsc_fields
    logical, intent(in), optional :: fixed_dims(:)
    integer :: id
    type(mlf_rsc_int64_1d) :: rsc_int
    id = numFields%incNRsc()
    this%v(id)%r = rsc_int
    if(present(rsc_name)) call this%v(id)%set_str(mlf_NAME, rsc_name)
    if(present(rsc_fields)) call this%v(id)%set_str(mlf_FIELDS, rsc_fields)
    Associate(ip => this%v(id)%r)
      select type(ip)
      class is (mlf_rsc_int64_1d)
        if(present(data_handler)) then
          info = data_handler%getData(ip%V, rsc_name, [nd], fixed_dims)
          if(info<0) RETURN
        endif
        if(.NOT. allocated(ip%V)) ALLOCATE(ip%V(nd), stat=info)
        if(CheckNZ(info, "Error array allocation", [nd])) RETURN
        nd = size(ip%V)
        if(present(pnt)) pnt => ip%V
      end select
    End Associate
  End Function mlf_obj_addI64array

  Integer Function mlf_obj_addI32array(this, numFields, nd, pnt, rsc_name, rsc_fields, &
      data_handler, fixed_dims) result(info)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    integer(c_int32_t), intent(out), pointer, optional :: pnt(:)
    class(mlf_rsc_numFields), intent(inout) :: numFields
    integer :: id
    integer(c_int64_t), intent (inout) :: nd
    character(len=*,kind=c_char), optional :: rsc_name, rsc_fields
    logical, intent(in), optional :: fixed_dims(:)
    type(mlf_rsc_int32_1d) :: rsc_int
    id = numFields%incNRsc()
    this%v(id)%r = rsc_int
    if(present(rsc_name)) call this%v(id)%set_str(mlf_NAME, rsc_name)
    if(present(rsc_fields)) call this%v(id)%set_str(mlf_FIELDS, rsc_fields)
    Associate(ip => this%v(id)%r)
      select type(ip)
      class is (mlf_rsc_int32_1d)
        if(present(data_handler)) then
          info = data_handler%getData(ip%V, rsc_name, [nd], fixed_dims)
          if(info<0) RETURN
        endif
        if(.NOT. allocated(ip%V)) ALLOCATE(ip%V(nd), stat=info)
        if(CheckNZ(info, "Error array allocation", [nd])) RETURN
        nd = size(ip%V)
        if(present(pnt)) pnt => ip%V
      end select
    End Associate
  End Function mlf_obj_addI32array

  Integer Function mlf_obj_addRArray(this, numFields, nd, pnt, rsc_name, rsc_fields, &
      data_handler, fixed_dims) result(info)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(out), pointer, optional :: pnt(:)
    class(mlf_rsc_numFields), intent(inout) :: numFields
    integer :: id
    integer(c_int64_t), intent (inout) :: nd
    character(len=*,kind=c_char), optional :: rsc_name, rsc_fields
    logical, intent(in), optional :: fixed_dims(:)
    type(mlf_rsc_double1d) :: rsc_double
    id = numFields%incNRsc()
    this%v(id)%r = rsc_double
    if(present(rsc_name)) call this%v(id)%set_str(mlf_NAME, rsc_name)
    if(present(rsc_fields)) call this%v(id)%set_str(mlf_FIELDS, rsc_fields)
    Associate(ip => this%v(id)%r)
      select type(ip)
      class is (mlf_rsc_double1d)
        if(present(data_handler)) then
          info = data_handler%getData(ip%V, rsc_name, [nd], fixed_dims)
          if(info<0) RETURN
        endif
        if(.NOT. allocated(ip%V)) ALLOCATE(ip%V(nd), stat=info)
        if(CheckNZ(info, "Error array allocation", [nd])) RETURN
        nd = size(ip%V)
        if(present(pnt)) pnt => ip%V
      end select
    End Associate
  End Function mlf_obj_addRArray

  Integer Function mlf_obj_addRMatrix(this, numFields, nd, pnt, rsc_name, rsc_fields, &
      data_handler, fixed_dims) result(info)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(out), pointer, optional :: pnt(:,:)
    class(mlf_rsc_numFields), intent(inout) :: numFields
    integer :: id
    integer(c_int64_t), intent (inout) :: nd(2)
    character(len=*,kind=c_char), optional :: rsc_name, rsc_fields
    logical, intent(in), optional :: fixed_dims(:)
    type(mlf_rsc_double2d) :: rsc_double
    id = numFields%incNRsc()
    this%v(id)%r = rsc_double
    if(present(rsc_name)) call this%v(id)%set_str(mlf_NAME, rsc_name)
    if(present(rsc_fields)) call this%v(id)%set_str(mlf_FIELDS, rsc_fields)
    Associate(ip => this%v(id)%r)
      select type(ip)
      class is (mlf_rsc_double2d)
        if(present(data_handler)) then
          info = data_handler%getData(ip%V, rsc_name, nd, fixed_dims)
          if(info<0) RETURN
        endif
        if(.NOT. allocated(ip%V)) ALLOCATE(ip%V(nd(1), nd(2)), stat=info)
        if(CheckNZ(info, "Error matrix allocation", nd)) RETURN
        nd = shape(ip%V)
        if(present(pnt)) pnt => ip%V
      end select
    End Associate
  End Function mlf_obj_addRMatrix

  Integer Function mlf_obj_addRMatrix3D(this, numFields, nd, pnt, rsc_name, rsc_fields, &
      data_handler, fixed_dims) result(info)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(out), pointer, optional :: pnt(:,:, :)
    class(mlf_rsc_numFields), intent(inout) :: numFields
    integer :: id
    integer(c_int64_t), intent (inout) :: nd(3)
    character(len=*,kind=c_char), optional :: rsc_name, rsc_fields
    logical, intent(in), optional :: fixed_dims(:)
    type(mlf_rsc_double3d) :: rsc_double
    id = numFields%incNRsc()
    this%v(id)%r = rsc_double
    if(present(rsc_name)) call this%v(id)%set_str(mlf_NAME, rsc_name)
    if(present(rsc_fields)) call this%v(id)%set_str(mlf_FIELDS, rsc_fields)
    Associate(ip => this%v(id)%r)
      select type(ip)
      class is (mlf_rsc_double3d)
        if(present(data_handler)) then
          info = data_handler%getData(ip%V, rsc_name, nd, fixed_dims) 
          if(info<0) RETURN
        endif
        if(.NOT. allocated(ip%V)) ALLOCATE(ip%V(nd(1), nd(2), nd(3)), stat=info)
        if(CheckNZ(info, "Error matrix allocation", nd)) RETURN
        nd = shape(ip%V)
        if(present(pnt)) pnt => ip%V
      end select
    End Associate
  End Function mlf_obj_addRMatrix3D

  Integer Function mlf_obj_addFArray(this, numFields, nd, pnt, rsc_name, rsc_fields, &
      data_handler, fixed_dims) result(info)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_float), intent(out), pointer, optional :: pnt(:)
    class(mlf_rsc_numFields), intent(inout) :: numFields
    integer :: id
    integer(c_int64_t), intent (inout) :: nd
    character(len=*,kind=c_char), optional :: rsc_name, rsc_fields
    logical, intent(in), optional :: fixed_dims(:)
    type(mlf_rsc_float1d) :: rsc_float
    id = numFields%incNRsc()
    this%v(id)%r = rsc_float
    if(present(rsc_name)) call this%v(id)%set_str(mlf_NAME, rsc_name)
    if(present(rsc_fields)) call this%v(id)%set_str(mlf_FIELDS, rsc_fields)
    Associate(ip => this%v(id)%r)
      select type(ip)
      class is (mlf_rsc_float1d)
        if(present(data_handler)) then
          info = data_handler%getData(ip%V, rsc_name, [nd], fixed_dims)
          if(info<0) RETURN
        endif
        if(.NOT. allocated(ip%V)) ALLOCATE(ip%V(nd), stat=info)
        if(CheckNZ(info, "Error array allocation", [nd])) RETURN
        nd = size(ip%V)
        if(present(pnt)) pnt => ip%V
      end select
    End Associate
  End Function mlf_obj_addFArray

  Integer Function mlf_obj_addFMatrix(this, numFields, nd, pnt, rsc_name, rsc_fields, &
      data_handler, fixed_dims) result(info)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_float), intent(out), pointer, optional :: pnt(:,:)
    class(mlf_rsc_numFields), intent(inout) :: numFields
    integer :: id
    integer(c_int64_t), intent (inout) :: nd(2)
    character(len=*,kind=c_char), optional :: rsc_name, rsc_fields
    logical, intent(in), optional :: fixed_dims(:)
    type(mlf_rsc_float2d) :: rsc_float
    id = numFields%incNRsc()
    this%v(id)%r = rsc_float
    if(present(rsc_name)) call this%v(id)%set_str(mlf_NAME, rsc_name)
    if(present(rsc_fields)) call this%v(id)%set_str(mlf_FIELDS, rsc_fields)
    Associate(ip => this%v(id)%r)
      select type(ip)
      class is (mlf_rsc_float2d)
        if(present(data_handler)) then
          info = data_handler%getData(ip%V, rsc_name, nd, fixed_dims)
          if(info<0) RETURN
        endif
        if(.NOT. allocated(ip%V)) ALLOCATE(ip%V(nd(1), nd(2)), stat=info)
        if(CheckNZ(info, "Error matrix allocation", nd)) RETURN
        nd = shape(ip%V)
        if(present(pnt)) pnt => ip%V
      end select
    End Associate
  End Function mlf_obj_addFMatrix

  Integer Function mlf_obj_addFMatrix3D(this, numFields, nd, pnt, rsc_name, rsc_fields, &
      data_handler, fixed_dims) result(info)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_float), intent(out), pointer, optional :: pnt(:,:, :)
    class(mlf_rsc_numFields), intent(inout) :: numFields
    integer :: id
    integer(c_int64_t), intent (inout) :: nd(3)
    character(len=*,kind=c_char), optional :: rsc_name, rsc_fields
    logical, intent(in), optional :: fixed_dims(:)
    type(mlf_rsc_float3d) :: rsc_float
    id = numFields%incNRsc()
    this%v(id)%r = rsc_float
    if(present(rsc_name)) call this%v(id)%set_str(mlf_NAME, rsc_name)
    if(present(rsc_fields)) call this%v(id)%set_str(mlf_FIELDS, rsc_fields)
    Associate(ip => this%v(id)%r)
      select type(ip)
      class is (mlf_rsc_float3d)
        if(present(data_handler)) then
          info = data_handler%getData(ip%V, rsc_name, nd, fixed_dims) 
          if(info<0) RETURN
        endif
        if(.NOT. allocated(ip%V)) ALLOCATE(ip%V(nd(1), nd(2), nd(3)), stat=info)
        if(CheckNZ(info, "Error matrix allocation", nd)) RETURN
        nd = shape(ip%V)
        if(present(pnt)) pnt => ip%V
      end select
    End Associate
  End Function mlf_obj_addFMatrix3D
      
  Integer Function mlf_arr_init(this, numFields, data_handler) result(info)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    class(mlf_rsc_numFields), intent(inout) :: numFields
    integer(c_int64_t) :: nIPar, nRPar, nRsc
    nIPar = numFields%nIPar; nRPar = numFields%nRPar; nRsc = numFields%nRsc+2
    ALLOCATE(this%v(nRsc))
    ! Reinit numFields
    numFields%nIPar = 0; numFields%nRPar = 0; numFields%nRsc = 0
    info = this%add_i64array(numFields, nIPar, this%iPar, C_CHAR_"iPar", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    if(CheckF(info, "Error creating ipar")) RETURN
    info = this%add_rarray(numFields, nRPar, this%rPar, C_CHAR_"rPar", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    if(CheckF(info, "Error creating rpar")) RETURN 
  End Function mlf_arr_init

  Subroutine mlf_obj_addRPar(this, numFields, pnt, field)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_rsc_numFields), intent(inout) :: numFields
    real(c_double), intent(out), pointer :: pnt
    character(len=*,kind=c_char), optional :: field
    numFields%nRPar = numFields%nRPar+1
    pnt => this%rPar(numFields%nRPar)
    if(present(field)) call this%v(2)%addField(field)
  End Subroutine mlf_obj_addRPar

  Subroutine mlf_obj_addIPar(this, numFields, pnt, field)
    class(mlf_arr_obj), intent(inout), target :: this
    class(mlf_rsc_numFields), intent(inout) :: numFields
    integer(c_int64_t), intent(out), pointer :: pnt
    character(len=*,kind=c_char), optional :: field
    numFields%nIPar = numFields%nIPar+1
    pnt => this%iPar(numFields%nIPar)
    if(present(field)) call this%v(1)%addField(field)
  End Subroutine mlf_obj_addIPar

  ! Generic interface for int arrays
  Function mlf_int32_getdata(array, nD, D, ptr) result(rptr)
    integer(c_int32_t), intent(in), target :: array(..)
    integer(c_int32_t) :: r
    integer(c_int), intent(out) :: nD, D(*)
    type(c_ptr) :: rptr, ptr
    integer :: i
    integer(c_size_t) :: N
    nD = rank(array)
    N = 1
    Do i=1,nD
      D(i) = size(array,i)
      N = N*D(i)
    End do
    rptr = C_LOC(array)
    if(C_ASSOCIATED(ptr)) then
      N = c_sizeof(r)*N
      rptr = c_memcpy(ptr, rptr, N)
    endif
  End Function mlf_int32_getdata
  Function mlf_int32_updatedata(array, ptr) result(info)
    integer(c_int32_t), intent(inout), target :: array(..)
    integer(c_int32_t) :: r
    type(c_ptr) :: rptr, ptr
    integer :: i, info, nD
    integer(c_size_t) :: N
    nD = rank(array)
    N = 1
    Do i=1,nD
      N = N*size(array,i)
    End do
    rptr = C_LOC(array)
    info = 0
    if(C_ASSOCIATED(ptr)) then
      N = c_sizeof(r)*N
      rptr = c_memcpy(rptr, ptr, N)
      if(.NOT. C_ASSOCIATED(rptr)) info = -1
    endif
  End Function mlf_int32_updatedata


  ! Generic interface for int(64) arrays
  Function mlf_int64_getdata(array, nD, D, ptr) result(rptr)
    integer(c_int64_t), intent(in), target :: array(..)
    integer(c_int64_t) :: r
    integer(c_int), intent(out) :: nD, D(*)
    type(c_ptr) :: rptr, ptr
    integer :: i
    integer(c_size_t) :: N
    nD = rank(array)
    N = 1
    Do i=1,nD
      D(i) = size(array,i)
      N = N*D(i)
    End do
    rptr = C_LOC(array)
    if(C_ASSOCIATED(ptr)) then
      N = c_sizeof(r)*N
      rptr = c_memcpy(ptr, rptr, N)
    endif
  End Function mlf_int64_getdata
  Function mlf_int64_updatedata(array, ptr) result(info)
    integer(c_int64_t), intent(inout), target :: array(..)
    integer(c_int64_t) :: r
    type(c_ptr) :: rptr, ptr
    integer :: i, info, nD
    integer(c_size_t) :: N
    nD = rank(array)
    N = 1
    Do i=1,nD
      N = N*size(array,i)
    End do
    rptr = C_LOC(array)
    info = 0
    if(C_ASSOCIATED(ptr)) then
      N = c_sizeof(r)*N
      rptr = c_memcpy(rptr, ptr, N)
      if(.NOT. C_ASSOCIATED(rptr)) info = -1
    endif
  End Function mlf_int64_updatedata

  ! Generic interface for double arrays
  Function mlf_double_getdata(array, nD, D, ptr) result(rptr)
    real(c_double), intent(in), target :: array(..)
    real(c_double) :: r
    integer(c_int), intent(out) :: nD, D(*)
    type(c_ptr) :: rptr, ptr
    integer :: i
    integer(c_size_t) :: N
    nD = rank(array)
    N = 1
    Do i=1,nD
      D(i) = size(array,i)
      N = N*D(i)
    End do
    rptr = C_LOC(array)
    if(C_ASSOCIATED(ptr)) then
      N = c_sizeof(r)*N
      rptr = c_memcpy(ptr, rptr, N)
    endif
  End Function mlf_double_getdata

  Function mlf_double_updatedata(array, ptr) result(info)
    real(c_double), intent(inout), target :: array(..)
    real(c_double) :: r
    type(c_ptr) :: rptr, ptr
    integer :: i, info, nD
    integer(c_size_t) :: N
    nD = rank(array)
    N = 1
    Do i=1,nD
      N = N*size(array,i)
    End do
    rptr = C_LOC(array)
    info = 0
    if(C_ASSOCIATED(ptr)) then
      N = c_sizeof(r)*N
      rptr = c_memcpy(rptr, ptr, N)
      if(.NOT. C_ASSOCIATED(rptr)) info = -1
    endif
  End Function mlf_double_updatedata
  
  ! Generic interface for float arrays
  Function mlf_float_getdata(array, nD, D, ptr) result(rptr)
    real(c_float), intent(in), target :: array(..)
    real(c_float) :: r
    integer(c_int), intent(out) :: nD, D(*)
    type(c_ptr) :: rptr, ptr
    integer :: i
    integer(c_size_t) :: N
    nD = rank(array)
    N = 1
    Do i=1,nD
      D(i) = size(array,i)
      N = N*D(i)
    End do
    rptr = C_LOC(array)
    if(C_ASSOCIATED(ptr)) then
      N = c_sizeof(r)*N
      rptr = c_memcpy(ptr, rptr, N)
    endif
  End Function mlf_float_getdata

  Function mlf_float_updatedata(array, ptr) result(info)
    real(c_float), intent(inout), target :: array(..)
    real(c_float) :: r
    type(c_ptr) :: rptr, ptr
    integer :: i, info, nD
    integer(c_size_t) :: N
    nD = rank(array)
    N = 1
    Do i=1,nD
      N = N*size(array,i)
    End do
    rptr = C_LOC(array)
    info = 0
    if(C_ASSOCIATED(ptr)) then
      N = c_sizeof(r)*N
      rptr = c_memcpy(rptr, ptr, N)
      if(.NOT. C_ASSOCIATED(rptr)) info = -1
    endif
  End Function mlf_float_updatedata

  ! Wrapper for the object int1d
  Function mlf_int32_1d_getd(this, nD, D, dt, ptr) result(cptr)
    class(mlf_rsc_int32_1d), intent(in), target :: this
    integer(c_int), intent(out) :: nD, D(*)
    type(mlf_dt), intent(out), pointer :: dt
    type(c_ptr) :: cptr, ptr
    if(associated(dt)) then
      dt%dt = mlf_INT; dt%acc = mlf_DIRECT
    endif
    cptr = mlf_int32_getdata(this%V, nD, D, ptr)
  End Function mlf_int32_1d_getd
  Function mlf_int32_1d_updated(this, ptr) result(info)
    class(mlf_rsc_int32_1d), intent(inout), target :: this
    type(c_ptr) :: ptr
    integer(c_int) :: info
    info = -1
    if(allocated(this%v)) info = mlf_int32_updatedata(this%V, ptr)
  End Function mlf_int32_1d_updated

  ! Wrapper for the object int1d(64)
  Function mlf_int64_1d_getd(this, nD, D, dt, ptr) result(cptr)
    class(mlf_rsc_int64_1d), intent(in), target :: this
    integer(c_int), intent(out) :: nD, D(*)
    type(mlf_dt), intent(out), pointer :: dt
    type(c_ptr) :: cptr, ptr
    if(associated(dt)) then
      dt%dt = mlf_INT64; dt%acc = mlf_DIRECT
    endif
    cptr = mlf_int64_getdata(this%V, nD, D, ptr)
  End Function mlf_int64_1d_getd

  Function mlf_int64_1d_updated(this, ptr) result(info)
    class(mlf_rsc_int64_1d), intent(inout), target :: this
    type(c_ptr) :: ptr
    integer(c_int) :: info
    info = -1
    if(allocated(this%v)) info = mlf_int64_updatedata(this%V, ptr)
  End Function mlf_int64_1d_updated

  ! Wrapper for the object double1d
  Function mlf_double1d_getd(this, nD, D, dt, ptr) result(cptr)
    class(mlf_rsc_double1d), intent(in), target :: this
    integer(c_int), intent(out) :: nD, D(*)
    type(mlf_dt), intent(out), pointer :: dt
    type(c_ptr) :: cptr, ptr
    if(associated(dt)) then
      dt%dt = mlf_DOUBLE; dt%acc = mlf_DIRECT
    endif
    cptr = mlf_double_getdata(this%V, nD, D, ptr)
  End Function mlf_double1d_getd

  Function mlf_double1d_updated(this, ptr) result(info)
    class(mlf_rsc_double1d), intent(inout), target :: this
    type(c_ptr) :: ptr
    integer(c_int) :: info
    info = -1
    if(allocated(this%v)) info = mlf_double_updatedata(this%V, ptr)
  End Function mlf_double1d_updated

  ! Wrappers for the object double2d
  Function mlf_double2d_getd(this, nD, D, dt, ptr) result(cptr)
    class(mlf_rsc_double2d), intent(in), target :: this
    integer(c_int), intent(out) :: nD, D(*)
    type(mlf_dt), intent(out), pointer :: dt
    type(c_ptr) :: cptr, ptr
    if(associated(dt)) then
      dt%dt = mlf_DOUBLE; dt%acc = mlf_DIRECT
    endif
    cptr = mlf_double_getdata(this%V, nD, D, ptr)
  End Function mlf_double2d_getd
  Function mlf_double2d_updated(this, ptr) result(info)
    class(mlf_rsc_double2d), intent(inout), target :: this
    type(c_ptr) :: ptr
    integer(c_int) :: info
    info = -1
    if(allocated(this%v)) info = mlf_double_updatedata(this%V, ptr)
  End Function mlf_double2d_updated
  ! Wrappers for the object double3d
  Function mlf_double3d_getd(this, nD, D, dt, ptr) result(cptr)
    class(mlf_rsc_double3d), intent(in), target :: this
    integer(c_int), intent(out) :: nD, D(*)
    type(mlf_dt), intent(out), pointer :: dt
    type(c_ptr) :: cptr, ptr
    if(associated(dt)) then
      dt%dt = mlf_DOUBLE; dt%acc = mlf_DIRECT
    endif
    cptr = mlf_double_getdata(this%V, nD, D, ptr)
  End Function mlf_double3d_getd
  Function mlf_double3d_updated(this, ptr) result(info)
    class(mlf_rsc_double3d), intent(inout), target :: this
    type(c_ptr) :: ptr
    integer(c_int) :: info
    info = -1
    if(allocated(this%v)) info = mlf_double_updatedata(this%V, ptr)
  End Function mlf_double3d_updated

  ! Wrapper for the object float1d
  Function mlf_float1d_getd(this, nD, D, dt, ptr) result(cptr)
    class(mlf_rsc_float1d), intent(in), target :: this
    integer(c_int), intent(out) :: nD, D(*)
    type(mlf_dt), intent(out), pointer :: dt
    type(c_ptr) :: cptr, ptr
    if(associated(dt)) then
      dt%dt = mlf_FLOAT; dt%acc = mlf_DIRECT
    endif
    cptr = mlf_float_getdata(this%V, nD, D, ptr)
  End Function mlf_float1d_getd

  Function mlf_float1d_updated(this, ptr) result(info)
    class(mlf_rsc_float1d), intent(inout), target :: this
    type(c_ptr) :: ptr
    integer(c_int) :: info
    info = -1
    if(allocated(this%v)) info = mlf_float_updatedata(this%V, ptr)
  End Function mlf_float1d_updated

  ! Wrappers for the object float2d
  Function mlf_float2d_getd(this, nD, D, dt, ptr) result(cptr)
    class(mlf_rsc_float2d), intent(in), target :: this
    integer(c_int), intent(out) :: nD, D(*)
    type(mlf_dt), intent(out), pointer :: dt
    type(c_ptr) :: cptr, ptr
    if(associated(dt)) then
      dt%dt = mlf_FLOAT; dt%acc = mlf_DIRECT
    endif
    cptr = mlf_float_getdata(this%V, nD, D, ptr)
  End Function mlf_float2d_getd
  Function mlf_float2d_updated(this, ptr) result(info)
    class(mlf_rsc_float2d), intent(inout), target :: this
    type(c_ptr) :: ptr
    integer(c_int) :: info
    info = -1
    if(allocated(this%v)) info = mlf_float_updatedata(this%V, ptr)
  End Function mlf_float2d_updated
  ! Wrappers for the object float3d
  Function mlf_float3d_getd(this, nD, D, dt, ptr) result(cptr)
    class(mlf_rsc_float3d), intent(in), target :: this
    integer(c_int), intent(out) :: nD, D(*)
    type(mlf_dt), intent(out), pointer :: dt
    type(c_ptr) :: cptr, ptr
    if(associated(dt)) then
      dt%dt = mlf_FLOAT; dt%acc = mlf_DIRECT
    endif
    cptr = mlf_float_getdata(this%V, nD, D, ptr)
  End Function mlf_float3d_getd
  Function mlf_float3d_updated(this, ptr) result(info)
    class(mlf_rsc_float3d), intent(inout), target :: this
    type(c_ptr) :: ptr
    integer(c_int) :: info
    info = -1
    if(allocated(this%v)) info = mlf_float_updatedata(this%V, ptr)
  End Function mlf_float3d_updated

End Module mlf_rsc_array

