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


! Contains procedures required for interfacting with C languages

Module mlf_intf
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_cfuns
  Use hdf5
  IMPLICIT NONE
  PRIVATE
  
  Public :: c_dealloc, c_allocate, c_getrsc, c_getinfo, c_updatersc
  Public :: mlf_init, mlf_quit
  Public :: mlf_obj_finalize, mlf_getobjfromc, mlf_getrawfromc

  ! Datatype descriptor for C handling functions
  Enum, bind(C)
    enumerator :: mlf_BOOL = 1, mlf_INT = 2, mlf_INT64 = 3, mlf_SIZE = 4
    enumerator :: mlf_FLOAT = 5, mlf_DOUBLE = 6, mlf_SIZEPARAM = 7, mlf_RAW = 8
  End Enum
  Enum, bind(C)
    enumerator :: mlf_DIRECT = 1, mlf_INDIRECT = 2, mlf_READONLY = 3
    enumerator :: mlf_COPYONLY = 4, mlf_WRITEONLY = 5
  End Enum
  Enum, bind(C)
    enumerator :: mlf_NAME = 1, mlf_DESC = 2, mlf_FIELDS = 3
  End Enum
  Enum, bind(C)
    enumerator :: mlf_OK = 0, mlf_UNINIT = -1, mlf_FUNERROR = -2, mlf_WRONGTYPE = -3
    enumerator :: mlf_WRONGRANK = -4, mlf_FILENOTFOUND = -5, mlf_FILEERROR = -6
    enumerator :: mlf_OTHERERROR = -7
  End Enum

  Public :: mlf_BOOL, mlf_INT, mlf_INT64, mlf_SIZE, mlf_FLOAT, mlf_DOUBLE, mlf_SIZEPARAM, mlf_RAW
  Public :: mlf_DIRECT, mlf_INDIRECT, mlf_READONLY, mlf_COPYONLY, mlf_WRITEONLY
  Public :: mlf_OK, mlf_UNINIT, mlf_FUNERROR, mlf_WRONGTYPE, mlf_NAME, mlf_DESC, mlf_FIELDS

  Type, Public, bind(C) :: mlf_dt
    integer(c_short) :: dt
    integer(c_short) :: acc
  End Type mlf_dt

  Type, Public :: mlf_universal
    class(*), pointer :: elt => NULL()
  End Type mlf_universal

  ! Generic ressource handler with c pointer structure
  Type, Public, abstract :: mlf_rsc_intf
  Contains
    procedure (mlf_cgetdata), deferred :: getd
    procedure (mlf_cupdatedata), deferred :: updated
  End Type mlf_rsc_intf

  ! Ressource handler
  Type, Public :: mlf_rsc
    character(:,kind=c_char), allocatable :: r_name, r_desc, r_fields
    class(mlf_rsc_intf), allocatable :: r
  Contains
    procedure :: set_str => mlf_rsc_setstr
    procedure :: get_str => mlf_rsc_getstr
  End Type mlf_rsc

  ! Ressource handler abstract class
  Type, Public :: mlf_obj
    type(mlf_rsc), allocatable :: v(:)
    type(mlf_universal), allocatable :: sub_objects(:)
    integer :: num_subobjects = 0
    character(:,kind=c_char), allocatable :: obj_name
  Contains
    ! FINAL not yet correctly implemented in GNU Fortran
    procedure :: finalize => mlf_obj_finalize
    procedure :: add_subobject => obj_add_subobject
    procedure :: del_subobject => obj_del_subobject
    procedure :: get_subobject => obj_get_subobject
  End Type mlf_obj

  Type, Public, extends(mlf_obj), abstract :: mlf_rsc_handler
  Contains
    procedure (mlf_rsc_handler_push), deferred :: pushState
  End Type mlf_rsc_handler

  ! Universal container. Permits the use of polymorphism
  Type, Public :: mlf_cintf
    class (*), pointer :: obj
    logical :: do_deallocate = .TRUE.
  End Type mlf_cintf

  Abstract Interface
    Function mlf_cgetdata(this, nD, D, dt, ptr)
      use iso_c_binding
      import :: mlf_rsc_intf, mlf_dt
      class(mlf_rsc_intf), intent(in), target :: this
      integer(c_int), intent(out) :: nD, D(*)
      type(mlf_dt), intent(out), pointer :: dt
      type(c_ptr) :: mlf_cgetdata, ptr
    End Function mlf_cgetdata
    Function mlf_cupdatedata(this, ptr)
      use iso_c_binding
      import :: mlf_rsc_intf
      class(mlf_rsc_intf), intent(inout), target :: this
      type(c_ptr) :: ptr
      integer(c_int) :: mlf_cupdatedata
    End Function mlf_cupdatedata

    Function mlf_rsc_handler_push(this, obj, override)
      use iso_c_binding
      import :: mlf_rsc_handler, mlf_obj
      class(mlf_rsc_handler), intent(inout), target :: this
      class(mlf_obj), intent(in), target :: obj
      integer(c_int) :: mlf_rsc_handler_push
      logical, intent(in), optional :: override
    End Function mlf_rsc_handler_push
  End Interface

Contains
  Function obj_get_subobject(this, name_obj) result(ptr)
    class(mlf_obj), intent(inout) :: this
    class(mlf_obj), pointer :: ptr
    character(len=*,kind=c_char) :: name_obj
    character(len=:,kind=c_char), allocatable, target :: nobj
    type(c_ptr) :: p1, p2
    integer :: i
    integer(c_size_t) :: lng
    ptr => NULL()
    nobj = name_obj//C_NULL_CHAR
    lng = len(nobj)
    p1 = c_loc(nobj)
    if(.NOT. allocated(this%sub_objects)) RETURN
    Do i=1,this%num_subobjects
      if(.NOT. associated(this%sub_objects(i)%elt)) CYCLE
      Associate(U => this%sub_objects(i)%elt)
        select type(U)
        class is (mlf_obj)
          ! Object shall have a name for getting access
          if(.NOT. allocated(U%obj_name)) CYCLE
          p2 = c_loc(U%obj_name)
          if(c_strncmp(p1, p2, lng)/=0) CYCLE
          ptr => U
          RETURN
        end select
      End Associate
    End Do
  End Function obj_get_subobject

  Integer Function obj_del_subobject(this, name_obj) result(info)
    class(mlf_obj), intent(inout) :: this
    class(mlf_obj), pointer :: ptr
    character(len=*,kind=c_char) :: name_obj
    character(len=:,kind=c_char), allocatable, target :: nobj
    type(c_ptr) :: p1, p2
    integer :: i
    integer(c_size_t) :: lng
    info = -1
    ptr => NULL()
    nobj = name_obj//C_NULL_CHAR
    lng = len(nobj)
    p1 = c_loc(nobj)
    if(.NOT. allocated(this%sub_objects)) RETURN
    Do i=1,this%num_subobjects
      if(.NOT. associated(this%sub_objects(i)%elt)) CYCLE
      Associate(U => this%sub_objects(i)%elt)
        select type(U)
        class is (mlf_obj)
          ! Object shall have a name for getting access
          if(.NOT. allocated(U%obj_name)) CYCLE
          p2 = c_loc(U%obj_name)
          if(c_strncmp(p1, p2, lng)/=0) CYCLE
          call U%finalize()
          DEALLOCATE(U)
          this%sub_objects(i)%elt => this%sub_objects(this%num_subobjects)%elt
          this%sub_objects(this%num_subobjects)%elt => NULL()
          this%num_subobjects = this%num_subobjects-1
          info = 0
          RETURN
        end select
      End Associate
    End Do
  End Function obj_del_subobject

  Subroutine obj_add_subobject(this, obj)
    class(mlf_obj), intent(inout) :: this
    class(mlf_obj), pointer :: obj
    integer :: nobj
    type(mlf_universal), allocatable :: temp(:)
    if(.NOT. allocated(this%sub_objects)) ALLOCATE(this%sub_objects(64))
    nobj = size(this%sub_objects)
    if(this%num_subobjects == nobj) then
      ALLOCATE(temp, source = this%sub_objects)
      nobj = nobj*2
      DEALLOCATE(this%sub_objects)
      ALLOCATE(this%sub_objects(nobj), source = temp)
    endif
    this%num_subobjects = this%num_subobjects+1
    this%sub_objects(this%num_subobjects)%elt => obj
  End Subroutine obj_add_subobject

  ! Workaround for permitting finalization of objects
  Subroutine mlf_obj_finalize(this)
    class(mlf_obj), intent(inout) :: this
    integer :: i
    if(.NOT. allocated(this%sub_objects)) RETURN
    Do i=1,this%num_subobjects
      if(.NOT. associated(this%sub_objects(i)%elt)) CYCLE
      Associate(U => this%sub_objects(i)%elt)
        select type(U)
        class is (mlf_obj)
          call U%finalize()
        end select
      End Associate
      deallocate(this%sub_objects(i)%elt)
    End Do
  End Subroutine mlf_obj_finalize

  ! Universal deallocator for C interface
  integer(c_int) Function c_dealloc(cptr) bind(C, name="mlf_dealloc")
    type(c_ptr), value :: cptr
    type(mlf_cintf), pointer :: this
    class (*), pointer :: obj
    c_dealloc = -1
    ! print *, "mlf_dealloc"
    if(.NOT. C_ASSOCIATED(cptr)) RETURN
    call C_F_POINTER(cptr, this)
    obj => this%obj
    select type(obj)
      class is (mlf_obj)
        ! FINAL not yet correctly implemented in GNU Fortran
        if(this%do_deallocate .AND. ASSOCIATED(this%obj)) then
          CALL obj%finalize()
          DEALLOCATE(obj)
        endif
      class default
        if(this%do_deallocate .AND. ASSOCIATED(this%obj)) DEALLOCATE(obj)
    end select
    DEALLOCATE(this)
    c_dealloc = 0
  End Function c_dealloc

  ! Utility function for C wrappers
  Function c_allocate(obj, do_deallocate) result(cptr)
    type(c_ptr) :: cptr
    type(mlf_cintf), pointer :: this
    class (*), pointer :: obj
    logical, optional, intent(in) :: do_deallocate
    ALLOCATE(this)
    this%obj => obj
    if(present(do_deallocate)) then
      this%do_deallocate = do_deallocate
    else
      this%do_deallocate = .TRUE. ! Default behaviour
    endif
    cptr = C_LOC(this)
  End Function c_allocate

  Function mlf_getrawfromc(cptr) result(obj)
    type(c_ptr) :: cptr
    type(mlf_cintf), pointer :: this
    class (*), pointer :: obj
    obj => NULL()
    if(.NOT. C_ASSOCIATED(cptr)) RETURN
    call C_F_POINTER(cptr, this)
    obj => this%obj
  End Function mlf_getrawfromc

  Function mlf_getobjfromc(cptr) result(obj)
    type(c_ptr) :: cptr
    type(mlf_cintf), pointer :: this
    class (mlf_obj), pointer :: obj
    obj => NULL()
    if(.NOT. C_ASSOCIATED(cptr)) RETURN
    call C_F_POINTER(cptr, this)
    if(.NOT. ASSOCIATED(this%obj)) RETURN
    associate(o => this%obj)
      select type (o)
        class is (mlf_obj)
          obj => o
      end select
    end associate
  End Function mlf_getobjfromc

  ! Get ressource handler for most of the classes
  integer(c_int) Function c_getnumrsc(cptr) result(r) bind(C, name="mlf_getnumrsc")
    type(c_ptr), value :: cptr
    class (mlf_obj), pointer :: obj
    obj => mlf_getobjfromc(cptr)
    r = -1
    if(.NOT. associated(obj)) RETURN
    if(.NOT. allocated(obj%v)) RETURN
    r = size(obj%v)
  End Function c_getnumrsc


  ! Get ressource handler for most of the classes
  type(c_ptr) Function c_getrsc(cptr, id, dt0, nD, D, dataptr) bind(C, name="mlf_getrsc")
    type(c_ptr), value :: cptr, dataptr, dt0
    class (mlf_obj), pointer :: obj
    integer(c_int), value :: id
    type(mlf_dt), pointer :: dt
    integer(c_int), intent(out) :: nD, D(*)
    obj => mlf_getobjfromc(cptr)
    c_getrsc = C_NULL_PTR
    if(.NOT. associated(obj)) RETURN
    if(.NOT. allocated(obj%v)) RETURN
    call C_F_POINTER(dt0, dt)
    Associate(v => obj%v)
      if(id < LBOUND(v,1) .OR. id > UBOUND(v,1)) RETURN
      if(.NOT. allocated(v(id)%r)) RETURN
      c_getrsc = v(id)%r%getd(nD, D, dt, dataptr)
    End Associate
  End Function c_getrsc
  
  ! C wrapper for getting information contained into strings
  type(c_ptr) Function c_getinfo(cptr, id, itype) bind(C, name="mlf_getinfo")
    type(c_ptr), value :: cptr
    integer(c_int), value :: id, itype
    class (mlf_obj), pointer :: obj
    c_getinfo = C_NULL_PTR
    obj => mlf_getobjfromc(cptr)
    if(.NOT. associated(obj)) RETURN
    if(.NOT. allocated(obj%v)) RETURN
    Associate(v => obj%v)
      if(id < LBOUND(v,1) .OR. id > UBOUND(v,1)) RETURN
      c_getinfo = v(id)%get_str(itype)
    End Associate
  End Function c_getinfo

  ! C wrapper for updating ressource content
  integer(c_int) Function c_updatersc(cptr, id, dataptr) bind(C, name="mlf_updatersc")
    type(c_ptr), value :: cptr, dataptr
    integer(c_int), value :: id
    class (mlf_obj), pointer :: obj
    c_updatersc = -1
    obj => mlf_getobjfromc(cptr)
    if(.NOT. associated(obj)) RETURN
    if(.NOT. allocated(obj%v)) RETURN
    Associate(v => obj%v)
      if(id < LBOUND(v,1) .OR. id > UBOUND(v,1)) RETURN
      if(.NOT. allocated(v(id)%r)) RETURN
      c_updatersc = v(id)%r%updated(dataptr)
    End Associate
  End Function c_updatersc

  integer(c_int) Function c_pushstate(cptr, cobj, override) bind(C, name="mlf_pushState")
    type(c_ptr), value :: cptr, cobj
    class (mlf_obj), pointer :: this, obj
    integer(c_int), value :: override
    c_pushstate = -1
    obj => mlf_getobjfromc(cobj)
    this => mlf_getobjfromc(cptr)
    if((.NOT. associated(obj)) .OR. (.NOT. associated(this))) RETURN
    select type(this)
    class is (mlf_rsc_handler)
      c_pushstate = this%pushState(obj, override /= 0)
    end select
  End Function c_pushstate

  ! Set a peculiar description string of a ressource
  Subroutine mlf_rsc_setstr(this, itype, str)
    class(mlf_rsc), intent(inout) :: this
    character(len=*,kind=c_char) :: str
    integer(c_int), intent(in) :: itype
    select case(itype)
    case(mlf_NAME)
      this%r_name = str//C_NULL_CHAR
    case(mlf_DESC)
      this%r_desc = str//C_NULL_CHAR
    case(mlf_FIELDS)
      this%r_fields = str//C_NULL_CHAR
    end select
  End Subroutine mlf_rsc_setstr

  ! Get a description string of a ressource
  type(c_ptr) Function mlf_rsc_getstr(this, itype) result(ptr)
    class(mlf_rsc), intent(in), target :: this
    integer(c_int), intent(in) :: itype
    ptr = C_NULL_PTR
    select case(itype)
    case(mlf_NAME)
      if(allocated(this%r_name)) ptr = C_LOC(this%r_name(1:1))
    case(mlf_DESC)
      if(allocated(this%r_desc)) ptr = C_LOC(this%r_desc(1:1))
    case(mlf_FIELDS)
      if(allocated(this%r_fields)) ptr = C_LOC(this%r_fields(1:1))
    end select
  End Function mlf_rsc_getstr

  Integer(c_int) Function mlf_init() result(error) bind(C, name="mlf_init")
    integer :: err
    call h5open_f(err)
    error = err
  End Function mlf_init

  Integer(c_int) Function mlf_quit() result(error) bind(C, name="mlf_quit")
    integer :: err
    call h5close_f(err)
    error = err
  End Function mlf_quit
End Module mlf_intf

