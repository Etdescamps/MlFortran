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
  Use hdf5
  IMPLICIT NONE
  PRIVATE

  integer, parameter :: h5_maxrank = 32

  Type, Public, abstract, extends(mlf_data_handler) :: mlf_hdf5_handler
  Contains
    procedure(mlf_hdf5_getId), deferred :: getId ! Get file or group identifier
    procedure :: open_group => mlf_hdf5_openGroup
    procedure :: get_subgroup => mlf_hdf5_getSubGroup
    procedure :: pushState => mlf_hdf5_pushState
    procedure :: getdata_double1d => mlf_hdf5_getDouble1d
    procedure :: getdata_double2d => mlf_hdf5_getDouble2d
    procedure :: getdata_double3d => mlf_hdf5_getDouble3d
    procedure :: getdata_int64 => mlf_hdf5_getInt64
    procedure :: getdata_int32 => mlf_hdf5_getInt32
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

  Function mlf_hdf5_openGroup(this, groupName, create) result(handler)
    class(mlf_hdf5_handler), intent(inout), target :: this
    class(mlf_hdf5_group), pointer :: handler
    class(mlf_obj), pointer :: np
    character(LEN=*), intent(in) :: groupName
    logical, optional :: create
    logical :: do_create = .FALSE.
    integer(HID_T) :: id
    integer :: hdferr
    handler => NULL()
    id = this%getId()
    if(id<0) RETURN
    if(present(create)) do_create = create
    if(do_create) then
      call H5GCreate_f(id, groupName, handler%group_id, hdferr)
    else
      call H5GOpen_f(id, groupName, handler%group_id, hdferr)
    endif
    if(hdferr<0) GOTO 10
    handler%obj_name = groupName // C_NULL_CHAR
    np => handler
    call this%add_subobject(np)
    RETURN
10  ALLOCATE(handler)
    DEALLOCATE(handler)
    handler => NULL()
  EndFunction mlf_hdf5_openGroup
  
  Function mlf_hdf5_getSubGroup(this, groupName) result(handler)
    class(mlf_hdf5_handler), intent(inout), target :: this
    class(mlf_hdf5_group), pointer :: handler
    character(LEN=*), intent(in) :: groupName
    class(mlf_obj), pointer :: np
    handler => NULL()
    np => this%get_subobject(groupName)
    if(ASSOCIATED(np)) then
      Select type(np)
      class is (mlf_hdf5_group)
        handler => np
      class default
        write (error_unit, *) 'Hdf5: subobject not from mlf_hdf5_group class', groupName
      End Select
    endif
    handler => this%open_group(groupName) ! Open existing group
    if(ASSOCIATED(handler)) RETURN
    handler => this%open_group(groupName, create= .TRUE.) ! Create new group
    if(ASSOCIATED(handler)) RETURN
  End Function mlf_hdf5_getSubGroup

  integer Function mlf_hdf5_createFile(this, fname, access_flag) result(hdferr)
    class(mlf_hdf5_file), intent(inout) :: this
    character(LEN=*), intent(in) :: fname
    integer, intent(in), optional :: access_flag
    ! Close file if the handler as been previously used
    call this%finalize()
    if(present(access_flag)) then
      call H5Fcreate_f(fname, access_flag, this%file_id, hdferr)
    else
      ! By default, delete the file
      call H5Fcreate_f(fname, H5F_ACC_TRUNC_F, this%file_id, hdferr)
    endif
    if(hdferr<0) write (error_unit, *) 'Error while creating file: ', fname
    this%obj_name = fname // C_NULL_CHAR
  End Function mlf_hdf5_createFile

  type(c_ptr) Function c_hdf5_createFile(pfname, trunk) result(cptr) bind(C, name="mlf_hdf5_createFile")
    type(mlf_hdf5_file), pointer :: this
    class (mlf_obj), pointer :: obj
    type(c_ptr), value :: pfname
    character(len=:, kind=c_char), allocatable, target :: fname
    integer(c_int), value :: trunk
    integer :: info
    ALLOCATE(this)
    cptr = c_null_ptr
    call mlf_stringFromC(pfname, fname)
    if(trunk == 0) then
      info = this%createFile(fname, H5F_ACC_EXCL_F)
    else
      info = this%createFile(fname)
    endif
    if(info<0) then
      deallocate(this); RETURN
    endif
    obj => this
    cptr = c_allocate(obj)
  End Function c_hdf5_createFile

  Integer Function mlf_hdf5_openFile(this, fname, access_flag) result(hdferr)
    class(mlf_hdf5_file), intent(inout) :: this
    character(LEN=*), intent(in) :: fname
    integer, intent(in), optional :: access_flag
    this%obj_name = fname
    ! Close file if the handler as been previously used
    call this%finalize()
    if(present(access_flag)) then
      call H5Fopen_f(fname, access_flag, this%file_id, hdferr)
    else
      ! By default, open in Read/Write mode
      call H5Fopen_f(fname, H5F_ACC_RDWR_F, this%file_id, hdferr)
    endif
    if(hdferr<0) write (error_unit, *) 'Error while opening file: ', fname
    this%obj_name = fname // C_NULL_CHAR
  End Function mlf_hdf5_openFile

  type(c_ptr) Function c_hdf5_openFile(pfname, rw) result(cptr) bind(C, name="mlf_hdf5_openFile")
    type(mlf_hdf5_file), pointer :: this
    class (mlf_obj), pointer :: obj
    type(c_ptr), value :: pfname
    character(len=:, kind=c_char), allocatable, target :: fname
    integer(c_int), value :: rw
    integer :: info
    ALLOCATE(this)
    cptr = c_null_ptr
    call mlf_stringFromC(pfname, fname)
    if(rw == 0) then
      info = this%openFile(fname, H5F_ACC_RDONLY_F)
    else
      info = this%openFile(fname)
    endif
    if(info<0) then
      deallocate(this); RETURN
    endif
    obj => this
    cptr = c_allocate(obj)
  End Function c_hdf5_openFile

  Subroutine mlf_hdf5_group_finalize(this)
    class(mlf_hdf5_group), intent(inout) :: this
    integer :: error
    if(this%group_id>=0) then
      call H5Gclose_f(this%group_id, error)
      if(error<0) write (error_unit, *) 'Error closing group (id): ', this%group_id
    endif
    call mlf_obj_finalize(this)
  End Subroutine mlf_hdf5_group_finalize


  Subroutine mlf_hdf5_file_finalize(this)
    class(mlf_hdf5_file), intent(inout) :: this
    logical :: valid
    integer :: error
    if(this%file_id>=0) then
      call H5Iis_valid_f(this%file_id, valid, error)
      if(valid) call H5Fclose_f(this%file_id, error)
      if(error<0) write (error_unit, *) 'Error closing file (id): ', this%file_id
    endif
    call mlf_obj_finalize(this)
  End Subroutine mlf_hdf5_file_finalize

  Integer Function mlf_hdf5_createOrWrite(file_id, dims, rname, h5_type, f_ptr, created) result(info)
    integer(HID_T), intent(in) :: file_id, h5_type
    character(len=*,kind=c_char), intent(in) :: rname
    integer(HSIZE_T), intent(in) :: dims(:)
    type(c_ptr), intent(in) :: f_ptr
    logical, intent(inout) :: created
    integer(HID_T) :: space_id = -1, data_id = -1
    info = createSpaceResource(file_id, dims, rname, h5_type, space_id, data_id, created)
    if(info >= 0) then
      call H5Dwrite_f(data_id, h5_type, f_ptr, info)
      if(info<0) write (error_unit, *) 'Error writing resource: '//rname
    endif
    info = min(info, closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_createOrWrite 

  Integer Function closeHDFIds(data_id, space_id, file_id) result(info)
    integer(HID_T), intent(in), optional :: file_id, space_id, data_id
    info = 0
    if(present(data_id)) then
      if(data_id>=0) call H5Sclose_f(space_id, info)
      if(CheckF(info, "Error closing dataset")) RETURN
    endif
    if(present(space_id)) then
      info = 0
      if(space_id>=0) call H5Dclose_f(data_id, info)
      if(CheckF(info, "Error closing dataspace")) RETURN
    endif
    if(present(file_id)) then
      info = 0
      if(file_id>=0) call H5Fclose_f(file_id, info)
      if(info<0) write (error_unit, *) 'Error closing file (id): ', file_id
    endif
  End Function closeHDFIds

  Integer Function createSpaceResource(file_id, dims, rname, h5_type, space_id, data_id, created) &
      result(info)
    integer(HID_T), intent(in) :: file_id, h5_type
    character(len=*,kind=c_char), intent(in) :: rname
    integer(HSIZE_T), intent(in) :: dims(:)
    integer(HID_T), intent(inout) :: space_id, data_id
    logical, intent(inout) :: created
    if(created) then
      call H5Gunlink_f(file_id, rname, info)
      if(CheckF(info, "Error removing dataset"//rname)) RETURN
    endif
    call H5Screate_simple_f(size(dims), dims, space_id, info)
    if(CheckF(info, "Error creating dataspace", dims)) RETURN
    call H5Dcreate_f(file_id, rname, h5_type, space_id, data_id, info)
    if(CheckF(info, "Error creating dataset"//rname)) RETURN
    created = .TRUE.
  End Function createSpaceResource

  Integer Function openSpaceResource(file_id, dims, rname, space_id, data_id, fixed_dims) result(info)
    integer(HID_T), intent(in) :: file_id ! Or groupe id
    character(len=*,kind=c_char), intent(in) :: rname
    integer(HSIZE_T), intent(inout) :: dims(:)
    logical, intent(in), optional :: fixed_dims(:)
    integer(HSIZE_T) :: max_dims(size(dims)), orig_dims(size(dims))
    integer :: frank0, frank
    integer(HID_T), intent(inout) :: space_id, data_id
    frank = size(dims)
    call H5Dopen_f(file_id, rname, data_id, info)
    if(CheckF(info, "Error opening resource")) RETURN
    call H5Dget_space_f(data_id, space_id, info)
    if(CheckF(info, "Error getting space_id")) RETURN
    call H5Sget_simple_extent_ndims_f(space_id, frank0, info)
    if(CheckF(info, "Error getting rank")) RETURN
    if(frank0 /= frank) then
      write (error_unit, *) "Error rank doesn't match:", frank, " /= ", frank0
      info = -1
      RETURN
    endif
    if(.NOT. present(fixed_dims)) orig_dims = dims
    call H5Sget_simple_extent_dims_f(space_id, dims, max_dims, info)
    if(CheckF(info, "Error getting dimensions")) RETURN
    if(.NOT. present(fixed_dims)) RETURN
    if(.NOT. any((orig_dims /= dims).AND.fixed_dims)) then
      write (error_unit, *) "Dimensions don't match:", pack(orig_dims, fixed_dims), &
        " /= ", pack(dims, fixed_dims)
      info = -1
      RETURN
    endif
  End Function openSpaceResource

  Integer Function mlf_hdf5_getDouble1d(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_double), intent(out), target, allocatable :: rsc_data(:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(1)
    character(len=*,kind=c_char) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id = -1, data_id = -1, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    gid = this%getId()
    if(gid<0) info=-1
    if(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    if(present(dims)) dims0 = dims
    info = openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    if(info >= 0) then
      allocate(rsc_data(dims0(1)))
      f_ptr = C_LOC(rsc_data)
      call H5Dread_f(data_id, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data"//rsc_name)
    endif
    info = min(info, closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getDouble1d

  Integer Function mlf_hdf5_getDouble2d(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_double), intent(out), target, allocatable :: rsc_data(:,:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(2)
    character(len=*,kind=c_char) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id = -1, data_id = -1, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    gid = this%getId()
    if(gid<0) info=-1
    if(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    if(present(dims)) dims0 = dims
    info = openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    if(info >= 0) then
      allocate(rsc_data(dims0(1), dims0(2)))
      f_ptr = C_LOC(rsc_data)
      call H5Dread_f(data_id, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data"//rsc_name)
    endif
    info = min(info, closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getDouble2d

  Integer Function mlf_hdf5_getDouble3d(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    real(c_double), intent(out), target, allocatable :: rsc_data(:,:,:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(3)
    character(len=*,kind=c_char) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id = -1, data_id = -1, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    gid = this%getId()
    if(gid<0) info=-1
    if(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    if(present(dims)) dims0 = dims
    info = openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    if(info >= 0) then
      allocate(rsc_data(dims0(1), dims0(2), dims0(3)))
      f_ptr = C_LOC(rsc_data)
      call H5Dread_f(data_id, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data"//rsc_name)
    endif
    info = min(info, closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getDouble3d

  Integer Function mlf_hdf5_getInt64(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    integer(c_int64_t), intent(out), target, allocatable :: rsc_data(:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(1)
    character(len=*,kind=c_char) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id = -1, data_id = -1, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    gid = this%getId()
    if(gid<0) info=-1
    if(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    if(present(dims)) dims0 = dims
    info = openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    if(info >= 0) then
      allocate(rsc_data(dims0(1)))
      f_ptr = C_LOC(rsc_data)
      call H5Dread_f(data_id, h5kind_to_type(c_int64_t, H5_INTEGER_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data"//rsc_name)
    endif
    info = min(info, closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getInt64
    
  Integer Function mlf_hdf5_getInt32(this, rsc_data, rsc_name, dims, fixed_dims) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    integer(c_int32_t), intent(out), target, allocatable :: rsc_data(:)
    integer(c_int64_t), intent(in), optional :: dims(:)
    integer(HSIZE_T) :: dims0(1)
    character(len=*,kind=c_char) :: rsc_name
    logical, intent(in), optional :: fixed_dims(:)
    integer(HID_T) :: space_id = -1, data_id = -1, gid
    type(c_ptr) :: f_ptr
    logical :: dumb
    gid = this%getId()
    if(gid<0) info=-1
    if(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    if(present(dims)) dims0 = dims
    info = openSpaceResource(gid, dims0, rsc_name, space_id, data_id, fixed_dims)
    if(info >= 0) then
      allocate(rsc_data(dims0(1)))
      f_ptr = C_LOC(rsc_data)
      call H5Dread_f(data_id, h5kind_to_type(c_int32_t, H5_INTEGER_KIND), f_ptr, info)
      dumb = CheckF(info, "Error reading data"//rsc_name)
    endif
    info = min(info, closeHDFIds(data_id, space_id))
  End Function mlf_hdf5_getInt32

  Integer(c_int) Function mlf_hdf5_pushState(this, obj) result(info)
    class(mlf_hdf5_handler), intent(inout), target :: this
    class(mlf_obj), intent(inout), target :: obj
    integer :: i
    type(c_ptr) :: f_ptr
    integer(HID_T) :: gid
    logical :: created = .FALSE.
    gid = this%getId()
    if(gid<0) info=-1
    if(CheckF(info, "mlf_hdf5: error getting id")) RETURN
    info = 0
    if(.NOT. allocated(obj%v)) RETURN
    do i=1,size(obj%v)
      if(.NOT. allocated(obj%v(i)%r)) CYCLE
      created = obj%v(i)%present_in_file
      associate(x => obj%v(i)%r)
        select type(x)
        class is (mlf_rsc_int64_1d)
          f_ptr = c_loc(x%V)
          info = mlf_hdf5_createOrWrite(gid, [size(x%V, kind=HSIZE_T)], &
            obj%v(i)%r_name, h5kind_to_type(c_int64_t, H5_INTEGER_KIND), f_ptr, created)
        class is (mlf_rsc_int32_1d)
          f_ptr = c_loc(x%V)
          info = mlf_hdf5_createOrWrite(gid, [size(x%V, kind=HSIZE_T)], &
            obj%v(i)%r_name, h5kind_to_type(c_int32_t, H5_INTEGER_KIND), f_ptr, created)
        class is (mlf_rsc_double1d)
          f_ptr = c_loc(x%V)
          info = mlf_hdf5_createOrWrite(gid, [size(x%V, kind=HSIZE_T)], &
            obj%v(i)%r_name, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, created)
        class is (mlf_rsc_double2d)
          f_ptr = c_loc(x%V)
          info = mlf_hdf5_createOrWrite(gid, shape(x%V, kind=HSIZE_T), &
            obj%v(i)%r_name, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, created)
        class is (mlf_rsc_double3d)
          f_ptr = c_loc(x%V)
          info = mlf_hdf5_createOrWrite(gid, shape(x%V, kind=HSIZE_T), &
            obj%v(i)%r_name, h5kind_to_type(c_double, H5_REAL_KIND), f_ptr, created)
        end select
      end associate
      obj%v(i)%present_in_file = created
    end do
  End Function mlf_hdf5_pushState
End Module mlf_hdf5

