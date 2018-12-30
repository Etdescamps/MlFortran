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

Module mlf_hdf5
  Use ieee_arithmetic
  Use iso_c_binding
  Use iso_fortran_env
  Use mlf_cfuns
  Use mlf_intf
  Use mlf_errors
  Use mlf_rsc_array
  Use mlf_utils
  Use hdf5
  IMPLICIT NONE
  PRIVATE

  integer, parameter :: h5_maxrank = 32

  Public :: mlf_hdf5_closeHDFIds, mlf_hdf5_createOrWrite, mlf_hdf5_openSpaceResource

  Type, Public, abstract, extends(mlf_data_handler) :: mlf_hdf5_handler
  Contains
    procedure(mlf_hdf5_getId), deferred :: getId ! Get file or group identifier
    procedure :: open_group => mlf_hdf5_openGroup
    procedure :: getSubObject => mlf_hdf5_getSubGroup
    procedure :: pushState => mlf_hdf5_pushState
    procedure :: getdata_double1d => mlf_hdf5_getDouble1d
    procedure :: getdata_double2d => mlf_hdf5_getDouble2d
    procedure :: getdata_double3d => mlf_hdf5_getDouble3d
    procedure :: getdata_float1d => mlf_hdf5_getFloat1d
    procedure :: getdata_float2d => mlf_hdf5_getFloat2d
    procedure :: getdata_float3d => mlf_hdf5_getFloat3d
    procedure :: getdata_int64 => mlf_hdf5_getInt64
    procedure :: getdata_int32 => mlf_hdf5_getInt32
    procedure :: mlf_pushdata_double2d, mlf_pushdata_double1d, mlf_pushdata_double3d
    procedure :: mlf_pushdata_float2d, mlf_pushdata_float1d, mlf_pushdata_float3d
    procedure :: mlf_pushdata_int32_1d, mlf_pushdata_int64_1d
    generic ::pushData => mlf_pushdata_double2d, mlf_pushdata_double1d, &
      mlf_pushdata_float2d, mlf_pushdata_float1d, mlf_pushdata_float3d, &
      mlf_pushdata_double3d, mlf_pushdata_int32_1d, mlf_pushdata_int64_1d
  End Type mlf_hdf5_handler

  Type, Public, extends(mlf_hdf5_handler) :: mlf_hdf5_group
    integer(HID_T) :: group_id = -1
  Contains
    procedure :: getId => mlf_hdf5_group_getId
    procedure :: finalize => mlf_hdf5_group_finalize
  End Type mlf_hdf5_group

  Type, Public, extends(mlf_hdf5_handler) :: mlf_hdf5_file
    integer(HID_T) :: file_id = -1
  Contains
    procedure :: getId => mlf_hdf5_file_getId
    procedure :: createFile => mlf_hdf5_createFile
    procedure :: openFile => mlf_hdf5_openFile
    procedure :: finalize => mlf_hdf5_file_finalize
  End Type mlf_hdf5_file

  Abstract Interface
    Function mlf_hdf5_getId(this) result(id)
      Use hdf5
      import :: mlf_hdf5_handler
      class(mlf_hdf5_handler), intent(in) :: this
      integer(HID_T) :: id
    End Function mlf_hdf5_getId
  End Interface

Contains
  Function mlf_hdf5_group_getId(this) result(id)
    class(mlf_hdf5_group), intent(in) :: this
    integer(HID_T) :: id
    id = this%group_id
  End Function mlf_hdf5_group_getId

  Function mlf_hdf5_file_getId(this) result(id)
    class(mlf_hdf5_file), intent(in) :: this
    integer(HID_T) :: id
    id = this%file_id
  End Function mlf_hdf5_file_getId

  Function mlf_hdf5_openGroup(this, groupName, create) Result(handler)
    class(mlf_hdf5_handler), intent(inout), target :: this
    type(mlf_hdf5_group), pointer :: handler
    class(mlf_obj), pointer :: np
    character(LEN=*), intent(in) :: groupName
    logical, optional :: create
    logical :: do_create
    integer(HID_T) :: id
    integer :: hdferr
    CALL InitOrDefault(do_create, .FALSE., create)
    handler => NULL()
    id = this%getId()
    If(id<0) RETURN
    ALLOCATE(handler)
    If(PresentAndTrue(create)) Then
      CALL H5GCreate_f(id, groupName, handler%group_id, hdferr)
    Else
      CALL H5GOpen_f(id, groupName, handler%group_id, hdferr)
    Endif
    If(hdferr<0) GOTO 10
    np => handler
    CALL this%add_subobject(groupName, np)
    RETURN
10  DEALLOCATE(handler)
    handler => NULL()
  EndFunction mlf_hdf5_openGroup
  
  Function mlf_hdf5_getSubGroup(this, obj_name) Result(ptr)
    class(mlf_hdf5_handler), intent(inout), target :: this
    class(mlf_data_handler), pointer :: ptr
    character(LEN=*,kind=c_char), intent(in) :: obj_name
    class(mlf_obj), pointer :: np
    ptr => NULL()
    np => this%get_subobject(obj_name)
    If(ASSOCIATED(np)) then
      sELECT TYPE(np)
      Class is (mlf_hdf5_group)
        ptr => np
      Class default
        write (error_unit, *) 'Hdf5: subobject not from mlf_hdf5_group class', obj_name
      END SELECT
    Endif
    ptr => this%open_group(obj_name) ! Open existing group
    If(ASSOCIATED(ptr)) RETURN
    ptr => this%open_group(obj_name, create= .TRUE.) ! Create new group
  End Function mlf_hdf5_getSubGroup

  integer Function mlf_hdf5_createFile(this, fname, access_flag) result(hdferr)
    class(mlf_hdf5_file), intent(inout) :: this
    character(LEN=*), intent(in) :: fname
    integer, intent(in), optional :: access_flag
    ! Close file if the handler as been previously used
    CALL this%finalize()
    If(present(access_flag)) Then
      CALL H5Fcreate_f(fname, access_flag, this%file_id, hdferr)
    Else
      ! By default, delete the file
      CALL H5Fcreate_f(fname, H5F_ACC_TRUNC_F, this%file_id, hdferr)
    Endif
    If(hdferr<0) WRITE (error_unit, *) 'Error while creating file: ', fname
  End Function mlf_hdf5_createFile

  type(c_ptr) Function c_hdf5_createFile(pfname, trunk) result(cptr) bind(C, name="mlf_hdf5_createFile")
    type(mlf_hdf5_file), pointer :: this
    class (*), pointer :: obj
    type(c_ptr), value :: pfname
    character(len=:, kind=c_char), allocatable, target :: fname
    integer(c_int), value :: trunk
    integer :: info
    ALLOCATE(this)
    cptr = c_null_ptr
    CALL mlf_stringFromC(pfname, fname)
    If(trunk == 0) then
      info = this%createFile(fname, H5F_ACC_EXCL_F)
    Else
      info = this%createFile(fname)
    Endif
    If(info<0) Then
      DEALLOCATE(this); RETURN
    Endif
    obj => this
    cptr = c_allocate(obj)
  End Function c_hdf5_createFile

  Integer Function mlf_hdf5_openFile(this, fname, access_flag) result(hdferr)
    class(mlf_hdf5_file), intent(inout) :: this
    character(LEN=*), intent(in) :: fname
    integer, intent(in), optional :: access_flag
    ! Close file if the handler as been previously used
    CALL this%finalize()
    If(PRESENT(access_flag)) then
      CALL H5Fopen_f(fname, access_flag, this%file_id, hdferr)
    Else
      ! By default, open in Read/Write mode
      CALL H5Fopen_f(fname, H5F_ACC_RDWR_F, this%file_id, hdferr)
    Endif
    If(hdferr<0) WRITE (error_unit, *) 'Error while opening file: ', fname
  End Function mlf_hdf5_openFile

  type(c_ptr) Function c_hdf5_openFile(pfname, rw) result(cptr) bind(C, name="mlf_hdf5_openFile")
    type(mlf_hdf5_file), pointer :: this
    class (*), pointer :: obj
    type(c_ptr), value :: pfname
    character(len=:, kind=c_char), allocatable, target :: fname
    integer(c_int), value :: rw
    integer :: info
    ALLOCATE(this)
    cptr = c_null_ptr
    CALL mlf_stringFromC(pfname, fname)
    If(rw == 0) Then
      info = this%openFile(fname, H5F_ACC_RDONLY_F)
    Else
      info = this%openFile(fname)
    Endif
    If(info<0) Then
      DEALLOCATE(this); RETURN
    Endif
    obj => this
    cptr = c_allocate(obj)
  End Function c_hdf5_openFile

  Subroutine mlf_hdf5_group_finalize(this)
    class(mlf_hdf5_group), intent(inout) :: this
    integer :: error
    If(this%group_id>=0) Then
      CALL H5Gclose_f(this%group_id, error)
      If(error<0) WRITE (error_unit, *) 'Error closing group (id): ', this%group_id
    Endif
    CALL mlf_obj_finalize(this)
  End Subroutine mlf_hdf5_group_finalize


  Subroutine mlf_hdf5_file_finalize(this)
    class(mlf_hdf5_file), intent(inout) :: this
    logical :: valid
    integer :: error
    If(this%file_id>=0) Then
      CALL H5Iis_valid_f(this%file_id, valid, error)
      If(valid) CALL H5Fclose_f(this%file_id, error)
      If(error<0) WRITE (error_unit, *) 'Error closing file (id): ', this%file_id
    Endif
    CALL mlf_obj_finalize(this)
  End Subroutine mlf_hdf5_file_finalize

  Integer Function mlf_hdf5_createOrWrite(file_id, dims, rname, h5_type, f_ptr, created) result(info)
    integer(HID_T), intent(in) :: file_id, h5_type
    character(len=*,kind=c_char), intent(in) :: rname
    integer(HSIZE_T), intent(in) :: dims(:)
    type(c_ptr), intent(in) :: f_ptr
    logical, intent(in) :: created
    integer(HID_T) :: space_id, data_id
    space_id = -1; data_id = -1
    info = createSpaceResource(file_id, dims, rname, h5_type, space_id, data_id, created)
    If(info >= 0) Then
      CALL H5Dwrite_f(data_id, h5_type, f_ptr, info)
      If(info<0) WRITE (error_unit, *) 'Error writing resource: '//rname
    Endif
    info = MIN(info, mlf_hdf5_closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_createOrWrite 

  Integer Function mlf_hdf5_closeHDFIds(data_id, space_id, file_id) result(info)
    integer(HID_T), intent(in), optional :: file_id, space_id, data_id
    info = 0
    If(PRESENT(data_id)) Then
      If(data_id>=0) CALL H5Sclose_f(space_id, info)
      If(CheckF(info, "Error closing dataset")) RETURN
    Endif
    If(PRESENT(space_id)) Then
      info = 0
      If(space_id>=0) CALL H5Dclose_f(data_id, info)
      If(CheckF(info, "Error closing dataspace")) RETURN
    Endif
    If(PRESENT(file_id)) Then
      info = 0
      If(file_id>=0) CALL H5Fclose_f(file_id, info)
      If(info<0) WRITE (error_unit, *) 'Error closing file (id): ', file_id
    Endif
  End Function mlf_hdf5_closeHDFIds

  Integer Function createSpaceResource(file_id, dims, rname, h5_type, space_id, data_id, created) &
      result(info)
    integer(HID_T), intent(in) :: file_id, h5_type
    character(len=*,kind=c_char), intent(in) :: rname
    integer(HSIZE_T), intent(in) :: dims(:)
    integer(HID_T), intent(inout) :: space_id, data_id
    logical, intent(in) :: created
    CALL H5Screate_simple_f(SIZE(dims), dims, space_id, info)
    If(CheckF(info, "Error creating dataspace: ", dims)) RETURN
    If(created) Then
      CALL H5Dopen_f(file_id, rname, data_id, info)
      If(CheckF(info, "Error opening resource "//rname)) RETURN
    Else
      CALL H5Dcreate_f(file_id, rname, h5_type, space_id, data_id, info)
      If(CheckF(info, "Error creating dataset: "//rname)) RETURN
    Endif
  End Function createSpaceResource

  Integer Function mlf_hdf5_openSpaceResource(file_id, dims, rname, space_id, data_id, fixed_dims) result(info)
    integer(HID_T), intent(in) :: file_id ! Or groupe id
    character(len=*,kind=c_char), intent(in) :: rname
    integer(HSIZE_T), intent(inout) :: dims(:)
    logical, intent(in), optional :: fixed_dims(:)
    integer(HSIZE_T) :: max_dims(SIZE(dims)), orig_dims(SIZE(dims))
    integer :: frank0, frank
    integer(HID_T), intent(inout) :: space_id, data_id
    info = 0
    frank = SIZE(dims)
    CALL H5Dopen_f(file_id, rname, data_id, info)
    If(CheckF(info, "Error opening resource")) RETURN
    CALL H5Dget_space_f(data_id, space_id, info)
    If(CheckF(info, "Error getting space_id")) RETURN
    CALL H5Sget_simple_extent_ndims_f(space_id, frank0, info)
    If(CheckF(info, "Error getting rank")) RETURN
    If(frank0 /= frank) Then
      WRITE (error_unit, *) "Error rank doesn't match:", frank, " /= ", frank0
      info = -1
      RETURN
    Endif
    If(.NOT. PRESENT(fixed_dims)) orig_dims = dims
    CALL H5Sget_simple_extent_dims_f(space_id, dims, max_dims, info)
    If(CheckF(info, "Error getting dimensions")) RETURN
    If(.NOT. PRESENT(fixed_dims)) RETURN
    If(.NOT. ANY((orig_dims /= dims).AND.fixed_dims)) Then
      WRITE (error_unit, *) "Dimensions don't match:", PACK(orig_dims, fixed_dims), &
        " /= ", PACK(dims, fixed_dims)
      info = -1
      RETURN
    Endif
  End Function mlf_hdf5_openSpaceResource

  Integer Function mlf_hdf5_getDouble1d(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_double), intent(out), target, allocatable :: rsc_data(:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(1)
    character(len=*,kind=c_char), intent(in) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id, data_id, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    info = 0; space_id = -1; data_id = -1
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    If(PRESENT(dims)) dims0 = dims
    info = mlf_hdf5_openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    If(info >= 0) Then
      ALLOCATE(rsc_data(dims0(1)))
      f_ptr = C_LOC(rsc_data)
      CALL H5Dread_f(data_id, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data: "//rsc_name)
    Endif
    info = MIN(info, mlf_hdf5_closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getDouble1d

  Integer Function mlf_hdf5_getDouble2d(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_double), intent(out), target, allocatable :: rsc_data(:,:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(2)
    character(len=*,kind=c_char), intent(in) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id, data_id, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    info = 0; space_id = -1; data_id = -1
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    If(PRESENT(dims)) dims0 = dims
    info = mlf_hdf5_openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    If(info >= 0) Then
      ALLOCATE(rsc_data(dims0(1), dims0(2)))
      f_ptr = C_LOC(rsc_data)
      CALL H5Dread_f(data_id, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data: "//rsc_name)
    Endif
    info = MIN(info, mlf_hdf5_closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getDouble2d

  Integer Function mlf_hdf5_getDouble3d(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_double), intent(out), target, allocatable :: rsc_data(:,:,:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(3)
    character(len=*,kind=c_char), intent(in) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id, data_id, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    info = 0; space_id = -1; data_id = -1
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    If(PRESENT(dims)) dims0 = dims
    info = mlf_hdf5_openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    If(info >= 0) Then
      ALLOCATE(rsc_data(dims0(1), dims0(2), dims0(3)))
      f_ptr = C_LOC(rsc_data)
      CALL H5Dread_f(data_id, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data: "//rsc_name)
    Endif
    info = MIN(info, mlf_hdf5_closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getDouble3d

  Integer Function mlf_hdf5_getFloat1d(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_float), intent(out), target, allocatable :: rsc_data(:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(1)
    character(len=*,kind=c_char), intent(in) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id, data_id, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    info = 0; space_id = -1; data_id = -1
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    If(PRESENT(dims)) dims0 = dims
    info = mlf_hdf5_openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    If(info >= 0) Then
      ALLOCATE(rsc_data(dims0(1)))
      f_ptr = C_LOC(rsc_data)
      CALL H5Dread_f(data_id, h5kind_to_type(c_float, H5_REAL_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data: "//rsc_name)
    Endif
    info = MIN(info, mlf_hdf5_closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getFloat1d

  Integer Function mlf_hdf5_getFloat2d(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_float), intent(out), target, allocatable :: rsc_data(:,:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(2)
    character(len=*,kind=c_char), intent(in) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id, data_id, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    info = 0; space_id = -1; data_id = -1
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    If(PRESENT(dims)) dims0 = dims
    info = mlf_hdf5_openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    If(info >= 0) Then
      ALLOCATE(rsc_data(dims0(1), dims0(2)))
      f_ptr = C_LOC(rsc_data)
      CALL H5Dread_f(data_id, h5kind_to_type(c_float, H5_REAL_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data: "//rsc_name)
    Endif
    info = MIN(info, mlf_hdf5_closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getFloat2d

  Integer Function mlf_hdf5_getFloat3d(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_float), intent(out), target, allocatable :: rsc_data(:,:,:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(3)
    character(len=*,kind=c_char), intent(in) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id, data_id, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    info = 0; space_id = -1; data_id = -1
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    If(PRESENT(dims)) dims0 = dims
    info = mlf_hdf5_openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    If(info >= 0) Then
      ALLOCATE(rsc_data(dims0(1), dims0(2), dims0(3)))
      f_ptr = C_LOC(rsc_data)
      CALL H5Dread_f(data_id, h5kind_to_type(c_float, H5_REAL_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data: "//rsc_name)
    Endif
    info = MIN(info, mlf_hdf5_closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getFloat3d

  Integer Function mlf_hdf5_getInt64(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    integer(c_int64_t), intent(out), target, allocatable :: rsc_data(:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(1)
    character(len=*,kind=c_char), intent(in) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id, data_id, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    info = 0; space_id = -1; data_id = -1
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    If(PRESENT(dims)) dims0 = dims
    info = mlf_hdf5_openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    If(info >= 0) Then
      ALLOCATE(rsc_data(dims0(1)))
      f_ptr = C_LOC(rsc_data)
      CALL H5Dread_f(data_id, h5kind_to_type(c_int64_t, H5_INTEGER_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data: "//rsc_name)
    Endif
    info = MIN(info, mlf_hdf5_closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getInt64
    
  Integer Function mlf_hdf5_getInt32(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    integer(c_int32_t), intent(out), target, allocatable :: rsc_data(:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(1)
    character(len=*,kind=c_char), intent(in) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id, data_id, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    info = 0; space_id = -1; data_id = -1
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    If(present(dims)) dims0 = dims
    info = mlf_hdf5_openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    If(info >= 0) Then
      ALLOCATE(rsc_data(dims0(1)))
      f_ptr = C_LOC(rsc_data)
      CALL H5Dread_f(data_id, h5kind_to_type(c_int32_t, H5_INTEGER_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data: "//rsc_name)
    Endif
    info = MIN(info, mlf_hdf5_closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getInt32

  Integer(c_int) Function mlf_pushdata_double2d(this, M, rname) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_double), target :: M(:,:)
    integer(HID_T) :: gid
    character(len=*,kind=c_char), intent(in) :: rname
    type(c_ptr) :: f_ptr
    info = 0
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    f_ptr = C_LOC(M)
    info = mlf_hdf5_createOrWrite(gid, SHAPE(M, kind=HSIZE_T), &
      rname, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, .FALSE.)
  End Function mlf_pushdata_double2d

  Integer(c_int) Function mlf_pushdata_double3d(this, M, rname) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_double), target :: M(:,:,:)
    integer(HID_T) :: gid
    character(len=*,kind=c_char), intent(in) :: rname
    type(c_ptr) :: f_ptr
    info = 0
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    f_ptr = C_LOC(M)
    info = mlf_hdf5_createOrWrite(gid, SHAPE(M, kind=HSIZE_T), &
      rname, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, .FALSE.)
  End Function mlf_pushdata_double3d

  Integer(c_int) Function mlf_pushdata_double1d(this, M, rname) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_double), target :: M(:)
    integer(HID_T) :: gid
    character(len=*,kind=c_char), intent(in) :: rname
    type(c_ptr) :: f_ptr
    info = 0
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    f_ptr = C_LOC(M)
    info = mlf_hdf5_createOrWrite(gid, [SIZE(M, kind=HSIZE_T)], &
      rname, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, .FALSE.)
  End Function mlf_pushdata_double1d

  Integer(c_int) Function mlf_pushdata_float2d(this, M, rname) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_float), target :: M(:,:)
    integer(HID_T) :: gid
    character(len=*,kind=c_char), intent(in) :: rname
    type(c_ptr) :: f_ptr
    info = 0
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    f_ptr = C_LOC(M)
    info = mlf_hdf5_createOrWrite(gid, SHAPE(M, kind=HSIZE_T), &
      rname, h5kind_to_type(c_float, H5_REAL_KIND), f_ptr, .FALSE.)
  End Function mlf_pushdata_float2d

  Integer(c_int) Function mlf_pushdata_float3d(this, M, rname) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_float), target :: M(:,:,:)
    integer(HID_T) :: gid
    character(len=*,kind=c_char), intent(in) :: rname
    type(c_ptr) :: f_ptr
    info = 0
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    f_ptr = C_LOC(M)
    info = mlf_hdf5_createOrWrite(gid, SHAPE(M, kind=HSIZE_T), &
      rname, h5kind_to_type(c_float, H5_REAL_KIND), f_ptr, .FALSE.)
  End Function mlf_pushdata_float3d

  Integer(c_int) Function mlf_pushdata_float1d(this, M, rname) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_float), target :: M(:)
    integer(HID_T) :: gid
    character(len=*,kind=c_char), intent(in) :: rname
    type(c_ptr) :: f_ptr
    info = 0
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    f_ptr = C_LOC(M)
    info = mlf_hdf5_createOrWrite(gid, [SIZE(M, kind=HSIZE_T)], &
      rname, h5kind_to_type(c_float, H5_REAL_KIND), f_ptr, .FALSE.)
  End Function mlf_pushdata_float1d

  Integer(c_int) Function mlf_pushdata_int32_1d(this, M, rname) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    integer(c_int32_t), target :: M(:)
    integer(HID_T) :: gid
    character(len=*,kind=c_char), intent(in) :: rname
    type(c_ptr) :: f_ptr
    info = 0
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    f_ptr = C_LOC(M)
    info = mlf_hdf5_createOrWrite(gid, [SIZE(M, kind=HSIZE_T)], &
      rname, h5kind_to_type(c_int32_t, H5_INTEGER_KIND), f_ptr, .FALSE.)
  End Function mlf_pushdata_int32_1d

  Integer(c_int) Function mlf_pushdata_int64_1d(this, M, rname) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    integer(c_int64_t), target :: M(:)
    integer(HID_T) :: gid
    character(len=*,kind=c_char), intent(in) :: rname
    type(c_ptr) :: f_ptr
    info = 0
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    f_ptr = C_LOC(M)
    info = mlf_hdf5_createOrWrite(gid, [SIZE(M, kind=HSIZE_T)], &
      rname, h5kind_to_type(c_int64_t, H5_INTEGER_KIND), f_ptr, .FALSE.)
  End Function mlf_pushdata_int64_1d

  Integer(c_int) Function Hdf5PushSubObject(this, obj, nameObj, override) Result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    class(mlf_obj), intent(in), target :: obj
    character(*, kind=c_char), intent(in) :: nameObj
    logical, intent(in), optional :: override
    class(mlf_hdf5_handler), pointer :: handler
    logical :: create
    info = -1
    create = .NOT. PresentAndTrue(override)
    handler => this%open_group(nameObj, create)
    If(.NOT. ASSOCIATED(handler)) RETURN
    info = handler%pushState(obj, override, subObjects = .TRUE.)
  End Function Hdf5PushSubObject

  Integer(c_int) Function mlf_hdf5_pushState(this, obj, override, subObjects) Result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    class(mlf_obj), intent(in), target :: obj
    integer :: i
    type(c_ptr) :: f_ptr
    integer(HID_T) :: gid
    logical, intent(in), optional :: override, subObjects
    logical :: ov
    info = 0
    gid = this%getId()
    If(gid<0) info=-1
    If(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    info = 0
    If(PresentAndTrue(subObjects) .AND. ALLOCATED(obj%sub_objects)) Then
      Do i = 1,SIZE(obj%sub_objects)
        If(.NOT. ASSOCIATED(obj%sub_objects(i)%elt)) CYCLE
        If(.NOT. ALLOCATED(obj%sub_objects(i)%obj_name)) CYCLE
        ASSOCIATE(p => obj%sub_objects(i)%elt, pname => obj%sub_objects(i)%obj_name)
          Select Type(p)
          Class is (mlf_obj)
            info = Hdf5PushSubObject(this, p, pname, override)
          End Select
        END ASSOCIATE
      End Do
    Endif
    If(.NOT. ALLOCATED(obj%v)) RETURN
    ov = PresentAndTrue(override)
    Do i=1,size(obj%v)
      If(.NOT. ALLOCATED(obj%v(i)%r)) CYCLE
      ASSOCIATE(x => obj%v(i)%r)
        SELECT TYPE(x)
        Class is (mlf_rsc_int64_1d)
          f_ptr = C_LOC(x%V)
          info = mlf_hdf5_createOrWrite(gid, [SIZE(x%V, KIND=HSIZE_T)], &
            obj%v(i)%r_name, h5kind_to_type(c_int64_t, H5_INTEGER_KIND), f_ptr, ov)
        Class is (mlf_rsc_int32_1d)
          f_ptr = C_LOC(x%V)
          info = mlf_hdf5_createOrWrite(gid, [SIZE(x%V, KIND=HSIZE_T)], &
            obj%v(i)%r_name, h5kind_to_type(c_int32_t, H5_INTEGER_KIND), f_ptr, ov)
        Class is (mlf_rsc_float1d)
          f_ptr = C_LOC(x%V)
          info = mlf_hdf5_createOrWrite(gid, [SIZE(x%V, KIND=HSIZE_T)], &
            obj%v(i)%r_name, h5kind_to_type(c_float, H5_REAL_KIND), f_ptr, ov)
        Class is (mlf_rsc_float2d)
          f_ptr = C_LOC(x%V)
          info = mlf_hdf5_createOrWrite(gid, SHAPE(x%V, KIND=HSIZE_T), &
            obj%v(i)%r_name, h5kind_to_type(c_float, H5_REAL_KIND), f_ptr, ov)
        Class is (mlf_rsc_float3d)
          f_ptr = C_LOC(x%V)
          info = mlf_hdf5_createOrWrite(gid, SHAPE(x%V, KIND=HSIZE_T), &
            obj%v(i)%r_name, h5kind_to_type(c_float, H5_REAL_KIND), f_ptr, ov)
        Class is (mlf_rsc_double1d)
          f_ptr = C_LOC(x%V)
          info = mlf_hdf5_createOrWrite(gid, [SIZE(x%V, KIND=HSIZE_T)], &
            obj%v(i)%r_name, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, ov)
        Class is (mlf_rsc_double2d)
          f_ptr = C_LOC(x%V)
          info = mlf_hdf5_createOrWrite(gid, SHAPE(x%V, KIND=HSIZE_T), &
            obj%v(i)%r_name, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, ov)
        Class is (mlf_rsc_double3d)
          f_ptr = C_LOC(x%V)
          info = mlf_hdf5_createOrWrite(gid, SHAPE(x%V, KIND=HSIZE_T), &
            obj%v(i)%r_name, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, ov)
        END SELECT
      END ASSOCIATE
    End Do
  End Function mlf_hdf5_pushState
End Module mlf_hdf5

