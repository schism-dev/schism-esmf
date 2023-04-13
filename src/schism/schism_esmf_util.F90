! This code is part of the SCHISM-ESMF interface, it defines utility
! functions used both by the NUOPC and ESMF caps
!
! @copyright 2022 Virginia Institute of Marine Science
! @copyright 2021-2022 Helmholtz-Zentrum Hereon
! @copyright 2018-2021 Helmholtz-Zentrum Geesthacht
!
! @author Carsten Lemmen <carsten.lemmen@hereon.de>
! @author Joseph Zhang <yjzhang@vims.edu>
! @author Richard Hofmeister
!
! @license Apache License, 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
! 		http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!

#define ESMF_CONTEXT  line=__LINE__,file=ESMF_FILENAME,method=ESMF_METHOD
#define ESMF_ERR_PASSTHRU msg="SCHISM subroutine call returned error"
#undef ESMF_FILENAME
#define ESMF_FILENAME "schism_esmf_util.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism_esmf_util

  use esmf
  implicit none

! The internal state saves data across ESMF phases and is
! persistent throught the lifetime of an instance.  Here, we
! only provide a boilerplate implementation of an empty internal state

  type type_PtrMap
#ifndef ESMF_NO_SEQUENCE
    sequence 
#endif
    integer(ESMF_KIND_I4), pointer  :: iarrayPtr1(:) => null()
    real(ESMF_KIND_R8), pointer  :: farrayPtr1(:) => null()
    real(ESMF_KIND_R8), pointer  :: farrayPtr2(:,:) => null()
    real(ESMF_KIND_R8), pointer  :: farrayPtr3(:,:,:) => null()
    character(len=ESMF_MAXSTR)   :: name
  end type

  type type_InternalState
#ifndef ESMF_NO_SEQUENCE
    sequence 
#endif   
    integer(ESMF_KIND_I4)          :: numForeignNodes=0
    integer(ESMF_KIND_I4)          :: numOwnedNodes=0
    integer(ESMF_KIND_I4), pointer :: ownedNodeIds(:) => null()
    integer(ESMF_KIND_I4), pointer :: foreignNodeIds(:) => null()
    type(ESMF_RouteHandle)         :: haloHandle
    type(type_PtrMap), allocatable :: ptrMap(:)
  end type

  type type_InternalStateWrapper
#ifndef ESMF_NO_SEQUENCE
    sequence ! why is this needed here? taken from documentation
#endif
    type(type_InternalState), pointer :: wrap => null()
  end type

  type(ESMF_MeshLoc) :: meshloc

  public meshloc
  public clockCreateFrmParam, SCHISM_FieldRealize
  public type_InternalState, type_InternalStateWrapper
  public SCHISM_StateFieldCreateRealize,SCHISM_StateFieldCreate
  !public SCHISM_StateGetField, 
  public SCHISM_StateImportWaveTensor 
  public SCHISM_MeshCreateNode
  public SCHISM_MeshCreateElement
  public SCHISM_InitializePtrMap
  !public SCHISM_FieldGet, SCHISM_FieldPut, 
  public SCHISM_StateUpdate, SCHISM_FieldPtrUpdate
  private 

  interface SCHISM_StateUpdate
    module procedure SCHISM_StateUpdate1
    module procedure SCHISM_StateUpdate2
    module procedure SCHISM_StateUpdate3
    module procedure SCHISM_StateUpdate4
  end interface 

contains

! #undef  ESMF_METHOD
! #define ESMF_METHOD "SCHISM_StateGetField"
! subroutine SCHISM_StateGetField(state, itemName, kwe, farrayPtr, field,  rc)

!   type(ESMF_State), intent(in)                        :: state
!   character(len=ESMF_MAXSTR), intent(in)              :: itemName
!   type(ESMF_KeywordEnforcer), intent(in), optional    :: kwe
!   type(ESMF_Field), intent(out), optional             :: field
!   real(ESMF_KIND_R8), pointer, intent(out), optional  :: farrayPtr(:)
!   integer(ESMF_KIND_I4), intent(out), optional        :: rc

!   integer(ESMF_KIND_I4)          :: rc_, localrc
!   type(ESMF_Field)               :: field_
!   character(len=ESMF_MAXSTR)     :: message
!   type(ESMF_StateItem_Flag)      :: itemType

!   localrc = ESMF_SUCCESS 
!   if (present(rc)) rc = ESMF_SUCCESS
!   if (.not.(present(field).or.(present(farrayPtr)))) return

!   call ESMF_StateGet(state, itemname=trim(itemName), itemType=itemType, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!   if (itemType /= ESMF_STATEITEM_FIELD) then 
!     write(message,'(A)') '--- field '//trim(itemName)//' could not be found in state'
!     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
!     localrc = ESMF_RC_NOT_FOUND
!     if (present(rc)) rc = localrc
!     return
!     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
!   endif
  
!   call ESMF_StateGet(state, itemname=trim(itemName), field=field_, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!   if (present(field)) field=field_
!   if (.not.present(farrayPtr)) return

!   call ESMF_FieldGet(field_, farrayptr=farrayPtr, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

! end subroutine SCHISM_StateGetField

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_InitializePtrMap"
subroutine SCHISM_InitializePtrMap(comp, kwe, rc)

  use schism_glbl, only: dav, pr2, tr_nd, eta2, windx2, windy2, npa
  use schism_glbl, only: uu2, vv2, srad, shum2, airt2, idry_e

  type(ESMF_GridComp), intent(inout)                  :: comp
  type(ESMF_KeywordEnforcer), intent(in), optional    :: kwe
  integer(ESMF_KIND_I4), intent(out), optional        :: rc

  integer(ESMF_KIND_I4)           :: rc_, localrc, i, ip 
  character(len=ESMF_MAXSTR)      :: message, name
  type(type_InternalStateWrapper) :: isPtr

  localrc = ESMF_SUCCESS 
  if (present(rc)) rc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, name=name, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_GridCompGetInternalState(comp, isPtr, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(isPtr%wrap%ptrMap(14))

  isPtr%wrap%ptrMap(1)%name = 'depth-averaged_x-velocity'
  !isPtr%wrap%ptrMap(1)%farrayPtr1 => dav(1,:)

  !> the following is not valid
  !do ip=1, npa
  !  isPtr%wrap%ptrMap(1)%farrayPtr1(ip)=>dav(1,ip)
  !enddo

  isPtr%wrap%ptrMap(2)%name = 'depth-averaged_y-velocity'
  isPtr%wrap%ptrMap(2)%farrayPtr1 => dav(2,:)

  isPtr%wrap%ptrMap(3)%name = 'air_pressure_at_sea_level'
  isPtr%wrap%ptrMap(3)%farrayPtr1 => pr2

  isPtr%wrap%ptrMap(4)%name = 'elevation_at_sea_level'
  isPtr%wrap%ptrMap(4)%farrayPtr1 => eta2

  isPtr%wrap%ptrMap(5)%name = 'air_temperature_at_sea_level'
  isPtr%wrap%ptrMap(5)%farrayPtr1 => airt2

  isPtr%wrap%ptrMap(6)%name = 'inst_zonal_wind_height10m'
  isPtr%wrap%ptrMap(6)%farrayPtr1 => windx2

  isPtr%wrap%ptrMap(7)%name = 'inst_zonal_wind_height10m'
  isPtr%wrap%ptrMap(7)%farrayPtr1 => windx2

  isPtr%wrap%ptrMap(8)%name = 'air_specific_humidity_at_sea_level'
  isPtr%wrap%ptrMap(8)%farrayPtr1 => shum2

  isPtr%wrap%ptrMap(9)%name = 'downwelling_shortwave_radiation_at_sea_level'
  isPtr%wrap%ptrMap(9)%farrayPtr1 => srad

  isPtr%wrap%ptrMap(11)%name = 'y-velocity_in_water'
  isPtr%wrap%ptrMap(11)%farrayPtr2 => vv2

  isPtr%wrap%ptrMap(10)%name = 'x-velocity_in_water'
  isPtr%wrap%ptrMap(10)%farrayPtr2 => uu2

  isPtr%wrap%ptrMap(12)%name = 'elevation_at_water_surface'
  isPtr%wrap%ptrMap(12)%farrayPtr1 => eta2

  isPtr%wrap%ptrMap(13)%name = 'tracer_concentration_in_water'
  isPtr%wrap%ptrMap(13)%farrayPtr3 => tr_nd

  isPtr%wrap%ptrMap(14)%name = 'ocean_mask'
  isPtr%wrap%ptrMap(14)%iarrayPtr1 => idry_e
  
end subroutine SCHISM_InitializePtrMap

! #undef  ESMF_METHOD
! #define ESMF_METHOD "SCHISM_FieldGet"
! subroutine SCHISM_FieldGet(field, isPtr, kwe, rc)

!   use schism_glbl, only : np, nvrt

!   type(ESMF_Field), intent(inout)                   :: field
!   type(type_InternalState), pointer, intent(in)     :: isPtr
!   type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe  
!   integer(ESMF_KIND_I4), intent(out), optional      :: rc

!   integer(ESMF_KIND_I4)          :: rc_, localrc, ip, i, rank
!   character(len=ESMF_MAXSTR)     :: message, name
!   logical                        :: isPresent
!   real(ESMF_KIND_R8), pointer    :: schismPtr1(:) => null()
!   real(ESMF_KIND_R8), pointer    :: schismPtr2(:) => null()
!   real(ESMF_KIND_R8), pointer    :: farrayPtr1(:) => null()
!   real(ESMF_KIND_R8), pointer    :: farrayPtr2(:,:) => null()
!   type(ESMF_TypeKind_Flag)       :: typeKind

!   localrc = ESMF_SUCCESS 
!   if (present(rc)) rc = ESMF_SUCCESS

!   isPresent = ESMF_FieldIsCreated(field, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!   if (.not.isPresent) then 
!     write(message, '(A)') '--- tried to access a field that is not created.'
!     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
!     localrc = ESMF_RC_ARG_BAD
!     if (present(rc)) rc = localrc
!     return
!     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
!   endif

!   call ESMF_FieldGet(field, name=name, typeKind=typeKind, rank=rank, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!   do i=1, ubound(isPtr%ptrMap,1)

!     write(message, '(A)') 'SCHISM_FieldGet '//trim(name)//' ?= '//trim(isPtr%ptrMap(i)%name)
!     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
  
!     if (trim(isPtr%ptrMap(i)%name) /= trim(name)) cycle 

!     if (rank == 1) then 
!       schismPtr1 => isPtr%ptrMap(i)%farrayPtr1
!     elseif (rank == 2) then 
!       !schismPtr2 => isPtr%ptrMap(i)%farrayPtr2
!     endif

!     exit
!   enddo 

!   if (.not.associated(schismPtr1)) then 
!     write(message,'(A)') '--- could not find ptrMap for '//trim(name)
!     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
!     return 
!   endif

!   write(message,'(A)') '--- found ptrMap for '//trim(name)
!   call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

!   isPresent = ESMF_RouteHandleIsCreated(isPtr%haloHandle, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!   if (isPresent) then 
!     !> @todo do we do this after or before assigning the variable
!     call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
!     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!     write(message,'(A)') '--- obtained halo route for field '//trim(name)
!     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

!   endif    

!   call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  
!   if (associated(schismPtr1)) then 
!     do ip = 1, isPtr%numOwnedNodes
!       schismPtr1(isPtr%ownedNodeIds(ip)) = farrayPtr1(ip)
!     end do
!   else
!     do ip = 1, isPtr%numOwnedNodes
!       !schismPtr2(isPtr%ownedNodeIds(ip),1:nvrt) = farrayPtr2(ip,1:nvrt)
!     end do
!   endif

!   write(message,'(A,5(X,F7.2))') 'FieldGet '//trim(name), schismPtr1(1:5)
!   call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

! end subroutine SCHISM_FieldGet

! #undef  ESMF_METHOD
! #define ESMF_METHOD "SCHISM_TracerPut"
! subroutine SCHISM_TracerPut(field, isPtr,  kwe, rc)

!   use schism_glbl, only : np, nvrt, tr_nd

!   type(ESMF_Field), intent(inout)                   :: field
!   type(type_InternalState), pointer, intent(in)     :: isPtr
!   type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe  
!   integer(ESMF_KIND_I4), intent(out), optional      :: rc

!   integer(ESMF_KIND_I4)          :: rc_, localrc, ip, i, rank
!   character(len=ESMF_MAXSTR)     :: message, name
!   logical                        :: isPresent
!   real(ESMF_KIND_R8), pointer    :: farrayPtr2(:,:) => null()
!   type(ESMF_TypeKind_Flag)       :: typeKind

!   localrc = ESMF_SUCCESS 
!   if (present(rc)) rc = ESMF_SUCCESS

!   isPresent = ESMF_FieldIsCreated(field, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!   if (.not.isPresent) then 
!     write(message, '(A)') '--- tried to access a field that is not created.'
!     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
!     localrc = ESMF_RC_ARG_BAD
!     if (present(rc)) rc = localrc
!     return
!     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
!   endif

!   call ESMF_FieldGet(field, name=name, typeKind=typeKind, rank=rank, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!   if (trim(name) == 'temperature') then 
!     i=1
!   elseif (trim(name) == 'salinity') then
!     i=2
!   else 
! #ifdef USE_FABM
!     !> @todo add code for FABM tracers
! #endif
!   endif

!   ! if (i<1 .or. i>ntracer) then 
!   !   write(message, '(A)') 'SCHISM_TracerPut '
!   !   call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
!   !   return
!   ! endif 
  
!   call ESMF_FieldGet(field, farrayPtr=farrayPtr2, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  
!   do ip=1, isPtr%numOwnedNodes
!     farrayPtr2(ip,1:nvrt) = tr_nd(i,isPtr%ownedNodeIds(ip),1:nvrt)  
!   end do

!   isPresent = ESMF_RouteHandleIsCreated(isPtr%haloHandle, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!   if (isPresent) then 
!     !> @todo do we do this after or before assigning the variable
!     call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
!     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
!   endif    

! end subroutine SCHISM_TracerPut

! #undef  ESMF_METHOD
! #define ESMF_METHOD "SCHISM_FieldPut"
! subroutine SCHISM_FieldPut(field, isPtr, kwe, rc)

!   type(ESMF_Field), intent(inout)                   :: field
!   type(type_InternalState), pointer, intent(in)     :: isPtr
!   type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe  
!   integer(ESMF_KIND_I4), intent(out), optional      :: rc

!   integer(ESMF_KIND_I4)          :: rc_, localrc, ip, i, rank
!   character(len=ESMF_MAXSTR)     :: message, name
!   logical                        :: isPresent
!   real(ESMF_KIND_R8), pointer    :: schismPtr1(:) => null()
!   real(ESMF_KIND_R8), pointer    :: farrayPtr1(:) => null()

!   localrc = ESMF_SUCCESS 
!   if (present(rc)) rc = ESMF_SUCCESS

!   isPresent = ESMF_FieldIsCreated(field, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!   if (.not.isPresent) then 
!     write(message, '(A)') 'Tried to access a field that is not created.'
!     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
!     localrc = ESMF_RC_ARG_BAD
!     if (present(rc)) rc = localrc
!     return
!     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
!   endif

!   call ESMF_FieldGet(field, name=name, rank=rank, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!   schismPtr1 => null()
!   do i=1, ubound(isPtr%ptrMap,1)

!     write(message, '(A)') 'SCHISM_FieldPut '//trim(name)//' ?= '//trim(isPtr%ptrMap(i)%name)
!     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
  
!     if (trim(isPtr%ptrMap(i)%name) /= trim(name)) cycle 

!     if (rank == 1) then 
!       schismPtr1 => isPtr%ptrMap(i)%farrayPtr1
!     elseif (rank == 2) then 
!       !schismPtr2 => isPtr%ptrMap(i)%farrayPtr2
!     endif

!     exit
!   enddo 

!   write(message,'(A)') '--- found ptrMap for '//trim(name)
!   call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

!   if (.not.associated(schismPtr1)) then 
!     write(message,'(A)') '--- could not find ptrMap for '//trim(name)
!     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
!     return 
!   endif

!   call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  
!   do ip = 1, isPtr%numOwnedNodes
!     farrayPtr1(ip) = schismPtr1(isPtr%ownedNodeIds(ip))
!   end do

!   isPresent = ESMF_RouteHandleIsCreated(isPtr%haloHandle, rc=localrc)
!   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!   if (isPresent) then 
!     call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
!     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!     write(message,'(A)') '--- obtained halo route for field '//trim(name)
!     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

!   endif    

! end subroutine SCHISM_FieldPut

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_FieldPtrUpdate"
subroutine SCHISM_FieldPtrUpdate(field, farray, kwe, isPtr, rc)

  type(ESMF_Field), intent(inout)                   :: field
  real(ESMF_KIND_R8),  intent(inout), target        :: farray(:)
  type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe
  type(type_InternalState), pointer, intent(in)     :: isPtr
  
  integer(ESMF_KIND_I4), intent(out), optional      :: rc

  real(ESMF_KIND_R8), pointer    :: farrayPtr1(:) => null()
  integer(ESMF_KIND_I4)          :: rc_, localrc, ip
  character(len=ESMF_MAXSTR)     :: message, name
  logical                        :: isPresent

  localrc = ESMF_SUCCESS 
  if (present(rc)) rc = ESMF_SUCCESS

  isPresent = ESMF_FieldIsCreated(field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (.not.isPresent) then 
    write(message, '(A)') 'Tried to access a field that is not created.'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_ARG_BAD
    if (present(rc)) rc = localrc
    return
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  call ESMF_FieldGet(field, name=name, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isPresent = ESMF_RouteHandleIsCreated(isPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (isPresent) then 
    call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A)') '--- obtained halo route for field '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  endif    

  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  
  do ip = 1, isPtr%numOwnedNodes
    farray(isPtr%ownedNodeIds(ip)) = farrayPtr1(ip)
  end do

  !> The following seems to overwrite valid data, so we cannot do this
  !do ip = 1,isPtr%numForeignNodes
  !  farray(isPtr%foreignNodeIds(ip)) = farrayPtr1(ip+isPtr%numOwnedNodes)
  !end do

  if (isPresent) then 
    !> @todo do we do this after or before assigning the variable
    call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A)') '--- obtained halo route for field '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  endif    

end subroutine SCHISM_FieldPtrUpdate

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateUpdate1"
subroutine SCHISM_StateUpdate1(state, name, farray, kwe, isPtr, rc)

  type(ESMF_State), intent(inout)                   :: state
  character(len=*), intent(in)                      :: name
  real(ESMF_KIND_R8),  intent(inout), target        :: farray(:)
  type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe
  type(type_InternalState), pointer, intent(in)     :: isPtr
  
  integer(ESMF_KIND_I4), intent(out), optional      :: rc

  type(ESMF_Field)               :: field
  real(ESMF_KIND_R8), pointer    :: farrayPtr1(:) => null()
  integer(ESMF_KIND_I4)          :: rc_, localrc, ip
  character(len=ESMF_MAXSTR)     :: message
  logical                        :: isPresent
  type(ESMF_StateIntent_Flag)    :: intent
  type(ESMF_StateItem_Flag)      :: itemType

  localrc = ESMF_SUCCESS 
  if (present(rc)) rc = ESMF_SUCCESS

  call ESMF_StateGet(state, itemname=trim(name), itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (itemType /= ESMF_STATEITEM_FIELD) then 
    write(message,'(A)') '--- SCHISM_StateUpdate1 skipped non-field '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    return
  endif 

  call ESMF_StateGet(state, itemname=trim(name), field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  
  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isPresent = ESMF_RouteHandleIsCreated(isPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_StateGet(state, stateintent=intent, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (isPresent) then 
    call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif    

  if (intent == ESMF_STATEINTENT_IMPORT) then 

    do ip = 1, isPtr%numOwnedNodes
      farray(isPtr%ownedNodeIds(ip)) = farrayPtr1(ip)
    end do

    write(message,'(A)') '--- SCHISM_StateUpdate1 imported '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  elseif (intent == ESMF_STATEINTENT_EXPORT) then 

    do ip = 1, isPtr%numOwnedNodes
       farrayPtr1(ip) = farray(isPtr%ownedNodeIds(ip))
    end do

    write(message,'(A)') '--- SCHISM_StateUpdate1 exported '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  else 
    write(message,'(A)') '--- SCHISM_StateUpdate1 skipped unspecified intent'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif    

  if (isPresent) then 
    call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif 

end subroutine SCHISM_StateUpdate1

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateUpdate2"
subroutine SCHISM_StateUpdate2(state, name, farray, kwe, isPtr, rc)

  type(ESMF_State), intent(inout)                   :: state
  character(len=*), intent(in)                      :: name
  real(ESMF_KIND_R8),  intent(inout), target        :: farray(:,:)
  type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe
  type(type_InternalState), pointer, intent(in)     :: isPtr
  
  integer(ESMF_KIND_I4), intent(out), optional      :: rc

  type(ESMF_Field)               :: field
  real(ESMF_KIND_R8), pointer    :: farrayPtr1(:) => null()
  integer(ESMF_KIND_I4)          :: rc_, localrc, ip
  character(len=ESMF_MAXSTR)     :: message
  logical                        :: isPresent
  type(ESMF_StateIntent_Flag)    :: intent
  type(ESMF_StateItem_Flag)      :: itemType

  localrc = ESMF_SUCCESS 
  if (present(rc)) rc = ESMF_SUCCESS

  if (ubound(farray,1) - lbound(farray,1) > 0) then 
    write(message,'(A)') '--- SCHISM_StateUpdate2 skipped non-degenerate '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    return
  endif 

  call ESMF_StateGet(state, itemname=trim(name), itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (itemType /= ESMF_STATEITEM_FIELD) then 
    write(message,'(A)') '--- SCHISM_StateUpdate2 skipped non-field '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    return
  endif 

  call ESMF_StateGet(state, itemname=trim(name), field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  
  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isPresent = ESMF_RouteHandleIsCreated(isPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_StateGet(state, stateintent=intent, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (intent == ESMF_STATEINTENT_IMPORT) then 

    if (isPresent) then 
      call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
    endif    

    do ip = 1, isPtr%numOwnedNodes
      farray(1,isPtr%ownedNodeIds(ip)) = farrayPtr1(ip)
    end do

    write(message,'(A)') '--- SCHISM_StateUpdate2 imported '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  elseif (intent == ESMF_STATEINTENT_EXPORT) then 

    do ip = 1, isPtr%numOwnedNodes
       farrayPtr1(ip) = farray(1,isPtr%ownedNodeIds(ip))
    end do

    if (isPresent) then 
      call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
    endif 

    write(message,'(A)') '--- SCHISM_StateUpdate2 exported '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  else 
    write(message,'(A)') '--- SCHISM_StateUpdate2 skipped unspecified intent'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif    

end subroutine SCHISM_StateUpdate2

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateUpdate3"
subroutine SCHISM_StateUpdate3(state, name, farray, kwe, isPtr, rc)

  use schism_glbl, only: nvrt 
  implicit none 

  type(ESMF_State), intent(inout)                   :: state
  character(len=*), intent(in)                      :: name
  real(ESMF_KIND_R8),  intent(inout), target        :: farray(:,:,:)
  type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe
  type(type_InternalState), pointer, intent(in)     :: isPtr
  
  integer(ESMF_KIND_I4), intent(out), optional      :: rc

  type(ESMF_Field)               :: field
  real(ESMF_KIND_R8), pointer    :: farrayPtr2(:,:) => null()
  integer(ESMF_KIND_I4)          :: rc_, localrc, ip
  character(len=ESMF_MAXSTR)     :: message
  logical                        :: isPresent
  type(ESMF_StateIntent_Flag)    :: intent
  type(ESMF_StateItem_Flag)      :: itemType

  localrc = ESMF_SUCCESS 
  if (present(rc)) rc = ESMF_SUCCESS

  if (ubound(farray,1) - lbound(farray,1) > 0) then 
    write(message,'(A)') '--- SCHISM_StateUpdate3 skipped non-degenerate '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    return
  endif 

  call ESMF_StateGet(state, itemname=trim(name), itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (itemType /= ESMF_STATEITEM_FIELD) then 
    write(message,'(A)') '--- SCHISM_StateUpdate3 skipped non-field '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    return
  endif 

  call ESMF_StateGet(state, itemname=trim(name), field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  
  call ESMF_FieldGet(field, farrayPtr=farrayPtr2, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isPresent = ESMF_RouteHandleIsCreated(isPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_StateGet(state, stateintent=intent, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (intent == ESMF_STATEINTENT_IMPORT) then 

    if (isPresent) then 
      call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
    endif    

    do ip = 1, isPtr%numOwnedNodes
      farray(1,isPtr%ownedNodeIds(ip),1:nvrt) = farrayPtr2(ip,1:nvrt)
    end do

    write(message,'(A)') '--- SCHISM_StateUpdate3 imported '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  elseif (intent == ESMF_STATEINTENT_EXPORT) then 

    do ip = 1, isPtr%numOwnedNodes
       farrayPtr2(ip,1:nvrt) = farray(1,isPtr%ownedNodeIds(ip),1:nvrt)
    end do

    if (isPresent) then 
      call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
    endif 

    write(message,'(A)') '--- SCHISM_StateUpdate exported '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  else 
    write(message,'(A)') '--- SCHISM_StateUpdate skipped unspecified intent'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif    

end subroutine SCHISM_StateUpdate3

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateUpdate4"
subroutine SCHISM_StateUpdate4(state, name, farray, kwe, isPtr, rc)

  type(ESMF_State), intent(inout)                   :: state
  character(len=*), intent(in)                      :: name
  integer(ESMF_KIND_I4),  intent(inout), target     :: farray(:)
  type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe
  type(type_InternalState), pointer, intent(in)     :: isPtr

  integer(ESMF_KIND_I4), intent(out), optional      :: rc

  type(ESMF_Field)               :: field
  integer(ESMF_KIND_I4), pointer :: farrayPtr1(:) => null()
  integer(ESMF_KIND_I4)          :: rc_, localrc, ip
  character(len=ESMF_MAXSTR)     :: message
  logical                        :: isPresent
  type(ESMF_StateIntent_Flag)    :: intent
  type(ESMF_StateItem_Flag)      :: itemType

  localrc = ESMF_SUCCESS
  if (present(rc)) rc = ESMF_SUCCESS

  call ESMF_StateGet(state, itemname=trim(name), itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (itemType /= ESMF_STATEITEM_FIELD) then
    write(message,'(A)') '--- SCHISM_StateUpdate1 skipped non-field '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    return
  endif

  call ESMF_StateGet(state, itemname=trim(name), field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isPresent = ESMF_RouteHandleIsCreated(isPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_StateGet(state, stateintent=intent, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (isPresent) then
    call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  if (intent == ESMF_STATEINTENT_IMPORT) then

    do ip = 1, isPtr%numOwnedNodes
      farray(isPtr%ownedNodeIds(ip)) = farrayPtr1(ip)
    end do

    write(message,'(A)') '--- SCHISM_StateUpdate4 imported '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  elseif (intent == ESMF_STATEINTENT_EXPORT) then

    do ip = 1, isPtr%numOwnedNodes
       farrayPtr1(ip) = farray(isPtr%ownedNodeIds(ip))
    end do

    write(message,'(A)') '--- SCHISM_StateUpdate4 exported '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  else
    write(message,'(A)') '--- SCHISM_StateUpdate4 skipped unspecified intent'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif

  if (isPresent) then
    call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

end subroutine SCHISM_StateUpdate4

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_GetPtr"
subroutine SCHISM_GetPtr(name, isPtr, farrayPtr, kwe, rc)

  character(len=ESMF_MAXSTR), intent(in)            :: name
  type(type_InternalState), pointer, intent(in)     :: isPtr
  real(ESMF_KIND_R8), pointer, intent(out)          :: farrayPtr(:)
  type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe
  integer(ESMF_KIND_I4), intent(out), optional      :: rc

  integer(ESMF_KIND_I4)          :: rc_, localrc, i
  character(len=ESMF_MAXSTR)     :: message

  localrc = ESMF_SUCCESS 
  if (present(rc)) rc = ESMF_SUCCESS
  farrayPtr => null()

  if (.not.associated(isPtr)) return
  if (.not.allocated(isPtr%ptrMap)) return

  do i=1, ubound(isPtr%ptrMap,1)
    if (trim(isPtr%ptrMap(i)%name) /= trim(name)) cycle
    farrayPtr => isPtr%ptrMap(i)%farrayPtr1
  enddo

end subroutine SCHISM_GetPtr

!> @todo separate into ESMF and NUOPC parts that reside in different source files
#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_FieldUpdate"
subroutine SCHISM_FieldUpdate(field, isPtr, kwe, rc)

  type(ESMF_Field), intent(inout)                   :: field
  type(type_InternalState), pointer, intent(in)     :: isPtr
  type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe
  integer(ESMF_KIND_I4), intent(out), optional      :: rc

  real(ESMF_KIND_R8), pointer    :: farrayPtr1(:) => null()
  real(ESMF_KIND_R8), pointer    :: schismPtr1(:) => null()
  integer(ESMF_KIND_I4)          :: rc_, localrc, ip
  character(len=ESMF_MAXSTR)     :: message, name
  logical                        :: isPresent

  localrc = ESMF_SUCCESS 
  if (present(rc)) rc = ESMF_SUCCESS

  isPresent = ESMF_FieldIsCreated(field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (.not.isPresent) then 
    write(message, '(A)') 'Tried to access a field that is not created.'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_ARG_BAD
    if (present(rc)) rc = localrc
    return
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  call ESMF_FieldGet(field, name=name, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isPresent = ESMF_RouteHandleIsCreated(isPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (isPresent) then 
    !> @todo do we do this after or before assigning the variable
    call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif    

  call SCHISM_GetPtr(name, isPtr, farrayPtr=schismPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  
  if (associated(farrayPtr1) .and. associated(schismPtr1)) then 
    do ip = 1, isPtr%numOwnedNodes
      schismPtr1(isPtr%ownedNodeIds(ip)) = farrayPtr1(ip)
    end do
  endif 

  if (isPresent) then 
    !> @todo do we do this after or before assigning the variable
    call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A)') '--- obtained halo route for field '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  endif    

end subroutine SCHISM_FieldUpdate

!> @todo separate into ESMF and NUOPC parts that reside in different source files
#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateFieldCreateRealize"
subroutine SCHISM_StateFieldCreateRealize(comp, state, name, field, kwe, rc)

  use schism_glbl, only: nws,np
  use NUOPC, only: NUOPC_Realize, NUOPC_IsConnected

  implicit none 

  type(ESMF_GridComp), intent(inout)               :: comp
  type(ESMF_State), intent(inout)                  :: state
  character(len=*), intent(in)                     :: name
  type(ESMF_Field), intent(out)                    :: field
  type(ESMF_KeywordEnforcer), intent(in), optional :: kwe
  integer(ESMF_KIND_I4), intent(out), optional     :: rc

  type(ESMF_Mesh)                    :: mesh
  type(ESMF_DistGrid)                :: distgrid
  integer(ESMF_KIND_I4)              :: localrc, rc_
  type(ESMF_Array)                   :: array 
  character(len=ESMF_MAXSTR)         :: compName, message
  type(type_InternalStateWrapper)    :: internalState
  type(type_InternalState), pointer  :: isDataPtr => null()
  type(ESMF_StateItem_Flag)          :: itemType
  logical                            :: isPresent
  
  rc_ = ESMF_SUCCESS
  localrc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, mesh=mesh, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (meshloc == ESMF_MESHLOC_NODE) then
    call ESMF_MeshGet(mesh, nodalDistgrid=distgrid, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  else
    call ESMF_MeshGet(mesh, elementDistgrid=distgrid, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  end if

  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isDataPtr => internalState%wrap
 
  write(message,'(A)') trim(compName)//' created array '//trim(name)//' on'
  if (isDataPtr%numForeignNodes > 0 .and. associated(isDataPtr%foreignNodeIds) .and. &
      meshloc == ESMF_MESHLOC_NODE) then 
    array = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, &
    haloSeqIndexList=isDataPtr%foreignNodeIds, &
    name=trim(name), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A,I5,A,I5,A)') trim(message)//' ',isDataPtr%numOwnedNodes,'/', np, &
    ' owned / resident nodes'
  else
    array = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, &
    name=trim(name), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A,I5,A)') trim(message)//' all resident nodes'
  endif
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  field = ESMF_FieldCreate(name=trim(name), mesh=mesh, array=array, &
    meshloc=meshloc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isPresent = ESMF_RouteHandleIsCreated(isDataPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (.not. isPresent) then 
    call ESMF_FieldHaloStore(field, routehandle=isDataPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A)') trim(compName)//' created halo route'
  endif    

  write(message,'(A)') trim(compName)//' created field '//trim(name)
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  call ESMF_StateGet(state, trim(name), itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (itemType /= ESMF_STATEITEM_FIELD) then 
    write(message,'(A)') trim(compName)//' skipped non-advertised field '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    return
  endif

  call ESMF_FieldHalo(field, routehandle=isDataPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message,'(A)') trim(compName)//' added halo route to field '//trim(name)
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  call NUOPC_Realize(state, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message,'(A)') trim(compName)//' realized field '//trim(name)
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  if (NUOPC_IsConnected(field, rc=localrc) .and. nws /= 3) then
    write(message, '(A,I1,A)') trim(compName)//' connected field '//trim(name)// &
      '  not used with nws = ', nws ,' (needs nws = 3)'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif

  if (present(rc)) rc = rc_

end subroutine SCHISM_StateFieldCreateRealize

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateFieldCreate"
subroutine SCHISM_StateFieldCreate(comp, state, name, field, kwe, rc)

  use schism_glbl, only: nws,np

  implicit none 

  type(ESMF_GridComp), intent(inout)               :: comp
  type(ESMF_State), intent(inout)                  :: state
  character(len=*), intent(in)                     :: name
  type(ESMF_Field), intent(out)                    :: field
  type(ESMF_KeywordEnforcer), intent(in), optional :: kwe
  integer(ESMF_KIND_I4), intent(out), optional     :: rc

  type(ESMF_Mesh)                    :: mesh
  type(ESMF_DistGrid)                :: distgrid
  integer(ESMF_KIND_I4)              :: localrc, rc_
  type(ESMF_Array)                   :: array 
  character(len=ESMF_MAXSTR)         :: compName, message
  type(type_InternalStateWrapper)    :: internalState
  type(type_InternalState), pointer  :: isDataPtr => null()
  type(ESMF_StateItem_Flag)          :: itemType
  logical                            :: isPresent
  
  rc_ = ESMF_SUCCESS
  localrc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, mesh=mesh, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (meshloc == ESMF_MESHLOC_NODE) then
    call ESMF_MeshGet(mesh, nodalDistgrid=distgrid, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  else
    call ESMF_MeshGet(mesh, elementDistgrid=distgrid, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  end if

  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isDataPtr => internalState%wrap
 
  write(message,'(A)') trim(compName)//' created array '//trim(name)//' on'
  if (isDataPtr%numForeignNodes > 0 .and. associated(isDataPtr%foreignNodeIds)) then 
    array = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, &
    haloSeqIndexList=isDataPtr%foreignNodeIds, &
    name=trim(name), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A,I5,A,I5,A)') trim(message)//' ',isDataPtr%numOwnedNodes,'/', np, &
    ' owned / resident nodes'
  else
    array = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, &
    name=trim(name), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A,I5,A)') trim(message)//' all resident nodes'
  endif
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  field = ESMF_FieldCreate(name=trim(name), mesh=mesh, array=array, &
    meshloc=meshloc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isPresent = ESMF_RouteHandleIsCreated(isDataPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (.not. isPresent) then 
    call ESMF_FieldHaloStore(field, routehandle=isDataPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A)') trim(compName)//' created halo route'
  endif    

  write(message,'(A)') trim(compName)//' created field '//trim(name)
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  call ESMF_StateAddReplace(state, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_StateGet(state, trim(name), itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (itemType /= ESMF_STATEITEM_FIELD) then 
    write(message,'(A)') trim(compName)//' skipped non-advertised field '//trim(name)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    return
  endif

  call ESMF_FieldHalo(field, routehandle=isDataPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message,'(A)') trim(compName)//' added halo route to field '//trim(name)
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  write(message,'(A)') trim(compName)//' created and added field '//trim(name)
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  if (present(rc)) rc = rc_

end subroutine SCHISM_StateFieldCreate

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_MeshCreateElement"
subroutine SCHISM_MeshCreateElement(comp, kwe, rc)
  use schism_glbl, only: pi
  use schism_glbl, only: np, npg, npa
  use schism_glbl, only: ne, neg, nea
  use schism_glbl, only: ylat, xlon
  use schism_glbl, only: ynd, xnd
  use schism_glbl, only: iplg, ipgl, idry, ics
  use schism_glbl, only: ielg, i34, idry_e, area, elnode

  implicit none

  type(ESMF_GridComp), intent(inout)               :: comp
  type(ESMF_KeywordEnforcer), intent(in), optional :: kwe
  integer(ESMF_KIND_I4), intent(out), optional     :: rc

  type(ESMF_Mesh)          :: mesh2d
  type(ESMF_DistGrid)      :: nodeDistgrid, elementDistgrid
  type(ESMF_CoordSys_Flag) :: coordsys
  type(ESMF_Field)         :: field
  type(ESMF_State)         :: exportState

  integer :: localPet
  integer :: indx, ip, ie, ii, nvcount
  character(len=ESMF_MAXSTR) :: message
  character(len=ESMF_MAXSTR) :: compName
  integer, dimension(:), allocatable            :: nodeids, elementids, nv
  real(ESMF_KIND_R8), dimension(:), allocatable :: nodecoords2d, elementcoords2d
  integer, dimension(:), allocatable            :: nodeowners, elementtypes
  integer, dimension(:), allocatable            :: nodemask, elementmask
  real(ESMF_KIND_R8), dimension(:), allocatable :: elementarea
  integer, dimension(1:4)                       :: elLocalNode
  integer :: rank2, localrc
  integer :: ownedCount, foreignCount
  real(ESMF_KIND_R8), parameter :: rad2deg=180.0d0/pi
  integer(ESMF_KIND_I4) :: rc_

  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()

  rc_ = ESMF_SUCCESS
  if (present(rc)) rc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, name=compName, localPet=localPet, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  ! allocate arrays, do not include ghosts
  allocate(  &
    nodeids(np), &
    nodecoords2d(2*np), &
    nodeowners(np), &
    nodemask(np), & 
    stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  allocate( &
    elementids(ne), &
    elementtypes(ne), & 
    elementmask(ne), &
    elementarea(ne), &
    elementcoords2d(2*ne), &
    stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)  

  ! set coordinate system type
  if (ics == 2) then
    coordsys=ESMF_COORDSYS_SPH_DEG
  else
    write(message, '(A)') trim(compName)//' uses a cartesian coordinate system'// &
       ' which may not be suitable for coupling'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    coordsys=ESMF_COORDSYS_CART
  endif

  ! fill arrays, nodes
  indx = 0
  do ip = 1, npa
    ! non-ghost nodes
    if (ip <= np) then
      ! ids
      indx = indx+1
      nodeids(indx) = iplg(ip)

      ! coordinates
      if (ics==2) then
        ! if geographical coordinates present
        nodecoords2d(2*indx-1) = rad2deg*xlon(ip) 
        nodecoords2d(2*indx)   = rad2deg*ylat(ip)
      else
        ! use cartesian coordinates
        nodecoords2d(2*indx-1) = xnd(ip)
        nodecoords2d(2*indx)   = ynd(ip)
      end if

      ! owner
      rank2 = ipgl(iplg(ip))%rank
      nodeowners(indx) = rank2
      if (associated(ipgl(iplg(ip))%next)) then
        if (ipgl(iplg(ip))%next%rank < rank2) then
          nodeowners(indx) = ipgl(iplg(ip))%next%rank
        end if
      end if

      ! mask
      nodemask(indx) = idry(ip)

      ! print out node related variables, just for debugging
      if (.false.) then
        write(message,'(A,2I8,2F10.3,I3)') trim(compName)//': nodeids, owner, x, y, mask = ', &
          nodeids(indx), nodeowners(indx), nodecoords2d(2*indx-1), nodecoords2d(2*indx), nodemask(indx)
        call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
      end if
    end if
  end do

  ! allocate array to store element connections
  allocate(nv(sum(i34(1:ne))), stat=localrc)

  ! fill arrays, elements
  indx = 0
  nvcount = 0
  do ie = 1, nea
    ! non-ghost elements
    if (ie <= ne) then
      ! ids
      indx = indx+1
      elementids(indx) = ielg(ie)

      ! element type
      if (i34(ie) == 3) then
        elementtypes(indx) = ESMF_MESHELEMTYPE_TRI
      else
        elementtypes(indx) = ESMF_MESHELEMTYPE_QUAD
      endif

      ! coordinates
      do ii = 1, i34(ie)
        elLocalNode(ii)=elnode(ii,ie)
        nvcount = nvcount+1
        nv(nvcount) = elnode(ii,ie) 
      end do
      elementcoords2d(2*indx-1) = sum(nodecoords2d(2*elLocalNode(1:i34(ie))-1))/i34(ie)
      elementcoords2d(2*indx)   = sum(nodecoords2d(2*elLocalNode(1:i34(ie))))/i34(ie)

      ! mask
      elementmask(indx) = idry_e(ie)

      ! area
      elementarea(indx) = area(ie)

      ! print out element related variables, just for debugging
      if (.false.) then
        write(message,'(A,I8,2F10.3,I3)') trim(compName)//': elementids, x, y, mask = ', &
          elementids(indx), elementcoords2d(2*indx-1), elementcoords2d(2*indx), elementmask(indx)
        call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
      end if
    end if
  end do

  ! create node distgrid
  nodeDistgrid = ESMF_DistgridCreate(nodeids,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  ! create element distgrid (distribute)
  elementDistgrid = ESMF_DistgridCreate(elementids,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  ! create mesh
  mesh2d = ESMF_MeshCreate(parametricDim=2, spatialDim=2, &
    nodeIds=nodeids, &
    nodeCoords=nodecoords2d, &
    nodeOwners=nodeOwners, &
    nodeMask=nodemask, &
    nodalDistgrid=nodeDistgrid, &
    elementIds=elementids, &
    elementTypes=elementtypes, &
    elementConn=nv, &
    elementMask=elementmask, &
    elementArea=elementarea, &
    elementCoords=elementcoords2d, &
    elementDistgrid=elementDistgrid, &
    coordSys=coordsys, &
    rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  ! write mesh in VTK format, just for debugging
  if (.false.) then 
    call ESMF_MeshWrite(mesh2d, filename="schism_mesh", rc=rc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  end if

  ! set component mesh
  call ESMF_GridCompSet(comp, mesh=mesh2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  write(message, '(A)') trim(compName)//' added mesh (element) to component'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! As the list of owned and non-owned nodes is not preserved in the ESMF_Mesh
  ! structure, we need to save this information to an internal state, for later 
  ! use in Array/Field creation.
  !call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !isDataPtr => internalState%wrap

  !isDataPtr%numOwnedNodes = 0
  !isDataPtr%numForeignNodes = 0
  !do ip = 1, npa
  !  ! non-ghost nodes
  !  if (ip <= np) then
  !    if (nodeowners(ip) == localPet) then 
  !      isDataPtr%numOwnedNodes = isDataPtr%numOwnedNodes + 1
  !    else 
  !      isDataPtr%numForeignNodes = isDataPtr%numForeignNodes + 1
  !    end if
  !  end if
  !end do

  !allocate(isDataPtr%ownedNodeIds(isDataPtr%numOwnedNodes), stat=localrc)
  !allocate(isDataPtr%foreignNodeIds(isDataPtr%numForeignNodes), stat=localrc)

  !ownedCount = 0
  !foreignCount = 0
  !do ip = 1, npa
  !  ! non-ghost nodes
  !  if (ip <= np) then
  !    if (nodeowners(ip) == localPet) then
  !      ownedCount=ownedCount + 1
  !      isDataPtr%ownedNodeIds(ownedCount) = ip
  !    else
  !      foreignCount=foreignCount + 1
  !      isDataPtr%foreignNodeIds(foreignCount) = ip
  !    endif
  !  end if
  !enddo

  ! add metadata
  !field = ESMF_FieldEmptyCreate(name='mesh_topology', rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_) 

  !call ESMF_AttributeSet(field, 'cf_role', 'mesh_topology', rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !call ESMF_AttributeSet(field, 'topology_dimension', 2, rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !call ESMF_AttributeSet(field, 'node_coordinates', 'mesh_node_lon mesh_node_lat', rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !call ESMF_AttributeSet(field, 'face_node_connectivity', 'mesh_element_node_connectivity', rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !call ESMF_GridCompGet(comp, exportState=exportState, rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !field = ESMF_FieldCreate(mesh2d, name='mesh_global_node_id', meshloc=ESMF_MESHLOC_NODE, typeKind=ESMF_TYPEKIND_I4, rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !call ESMF_FieldGet(field, farrayPtr=farrayPtrI41, rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !farrayPtrI41 = isDataPtr%ownedNodeIds

  !call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  ! clean up
  deallocate(nodeids, stat=localrc)
  deallocate(nodecoords2d, stat=localrc)
  deallocate(nodeowners, stat=localrc)
  deallocate(nodemask, stat=localrc)
  deallocate(elementids, stat=localrc)
  deallocate(elementtypes, stat=localrc)
  deallocate(elementmask, stat=localrc)
  deallocate(elementarea, stat=localrc)
  deallocate(elementCoords2d, stat=localrc)
  deallocate(nv, stat=localrc)

end subroutine SCHISM_MeshCreateElement

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_MeshCreateNode"
subroutine SCHISM_MeshCreateNode(comp, kwe, rc)
  ! Define ESMF domain partition
  use schism_glbl, only: pi, llist_type, elnode, i34, ipgl
  use schism_glbl, only: iplg, ielg, idry_e, idry, ynd, xnd
  !use schism_glbl, only: ylat, xlon, np, ne, ics
  use schism_glbl, only: ylat, xlon, np, npa, nea,  ics
  use schism_glbl, only:  nvrt, ne

  implicit none

  type(ESMF_GridComp), intent(inout)               :: comp
  type(ESMF_KeywordEnforcer), intent(in), optional :: kwe
  integer(ESMF_KIND_I4), intent(out), optional     :: rc

  type(ESMF_Mesh)          :: mesh2d
  type(ESMF_DistGrid)      :: nodeDistgrid, elementDistgrid, distgrid
  type(ESMF_CoordSys_Flag) :: coordsys
  integer, dimension(:), allocatable            :: nodeids, elementids, nv
  real(ESMF_KIND_R8), dimension(:), allocatable :: nodecoords2d
  real(ESMF_KIND_R8), dimension(:), allocatable :: elementcoords2d
  integer, dimension(:), allocatable            :: nodeowners, elementtypes
  integer, dimension(:), allocatable            :: nodemask, elementmask
  integer, dimension(:), allocatable            :: tmpIdx, tmpIdx2, localNodes, nodeHaloIdx
  integer, dimension(:), allocatable            :: schismTolocalNodes,testids
  integer, dimension(1:4)                       :: elLocalNode
  integer               :: numNodeHaloIdx
  integer               :: i,n,nvcount
  integer               :: ii,ip,ie, localrc, rc_
  integer               :: npo,neo,nef,npf,rank2
  type(llist_type),pointer :: nextp=>null()

  integer(ESMF_KIND_I4), pointer, dimension(:)  :: farrayPtrI41 => null()
  integer(ESMF_KIND_I4), pointer, dimension(:,:):: farrayPtrI42 => null()

  real(ESMF_KIND_R8), parameter :: rad2deg=180.0d0/pi

  integer(ESMF_KIND_I4),dimension(1:1) :: maskValues=(/1/)
  character(len=ESMF_MAXSTR)           :: compName, message, fieldName
  type(ESMF_State)                     :: exportState
  type(ESMF_Field)                     :: field

  integer(ESMF_KIND_I4) :: numOwnedNodes, numForeignNodes, localPet
  integer(ESMF_KIND_I4) :: ownedCount, foreignCount

  integer(ESMF_KIND_I4), allocatable, target :: ownedNodeIdx(:)
  integer(ESMF_KIND_I4), allocatable, target :: foreignNodeIdx(:)

  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()

!  write(0,*)'__LINE__ inside addSchismMesh'
!  call ESMF_Finalize() 

  rc_ = ESMF_SUCCESS
  if (present(rc)) rc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, name=compName, localPet=localPet, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !> @todo all non-ESMF stuff should be outsourced to schism_bmi.F90
  ! prepare mesh
  ! a) take local elements and get number of nodes for definition
  ! b) get number of augmented nodes, which are not connected to local elements
  ! c) allocate node arrays such that nodeids gives all nods belonging to local
  !    elements, all other nodes as part of the augmented domain are outside the
  !    exclusive region
  ! d) allocate element arrays such that local elements (ne) are in array and
  !    augmented elements are defined in the computational domain outside the
  !    exclusive domain

  !> A node is owned by same rank across PETs; interface nodes are owned by min rank
  !> @todo something is still off with npa = np
  allocate(  &
    nodeids(npa), &
    nodecoords2d(2*npa), &
    nodeowners(npa), &
    nodemask(npa), & 
    stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  allocate( &
    elementids(nea), &
    elementtypes(nea), & 
    elementmask(nea), &
    elementcoords2d(2*nea), &
    stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  ! nv (elemConn): 1D array for connectivity (packed from 2D array elnode).
  ! Outputs local node # 
  allocate(nv(sum(i34(1:nea))), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  ! set ESMF coordSys type
  if (ics==2) then
    coordsys=ESMF_COORDSYS_SPH_DEG
  else
    write(message, '(A)') trim(compName)//' uses a cartesian coordinate system'// &
       ' which may not be suitable for coupling'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    coordsys=ESMF_COORDSYS_CART
  endif

  do ip=1, npa
    nodeids(ip)=iplg(ip) !global node #
    if (ics==2) then
      ! if geographical coordinates present
      nodecoords2d(2*ip-1) = rad2deg*xlon(ip)
      nodecoords2d(2*ip)   = rad2deg*ylat(ip)
    else
      ! use cartesian coordinates
      nodecoords2d(2*ip-1) = xnd(ip)
      nodecoords2d(2*ip)   = ynd(ip)
    end if

    !nodeowners must be unique cross all PETs
    ! In Shinecook we have 3070 nodes, on two PET that should be 1558 and 1512 owned
    ! We actually get that correct with np, BUT with npa we get 1600 owned on PET0
 
    ! Resident nodes (ip <= np) can be owned or foreign, i.e. owned by the PET that has 
    ! the lowest rank in the linked list ipgl(iplg(ip))%next
    if (ip <= np) then 

      rank2=ipgl(iplg(ip))%rank
 
      ! If it is a resident node, we only need the first next, since the list is ascending
      nodeowners(ip) = rank2 !init 
      if(associated(ipgl(iplg(ip))%next)) then !interface or ghost node
        if(ipgl(iplg(ip))%next%rank<rank2) then
          nodeowners(ip) = ipgl(iplg(ip))%next%rank
        endif
      endif

    else ! ip > np ! ghost node, i.e. definitely not owned

      rank2 = huge(1) ! or at least number of processes + 1
      nextp => ipgl(iplg(ip))
      do while (associated(nextp%next))
        nextp => nextp%next
        if (nextp%rank < rank2) rank2 = nextp%rank 
      end do
    
      if (rank2 > huge(1) - 1 ) then ! greater than number of PET
        write(message,'(A,I5)') trim(compName)//' could not find owner of node ', ip
        localrc = ESMF_RC_ARG_SIZE
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
      endif
    
      if (rank2 == localPet) then
        write(message,'(A,I5)') trim(compName)//' found myself wrongly as owner of node', ip
        localrc = ESMF_RC_ARG_SIZE
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
      endif
    
      nodeowners(ip) = rank2
    endif 

    nodemask(ip) = idry(ip)
  end do ! ip...npa

  ! As the list of owned and non-owned nodes is not preserved in the ESMF_Mesh
  ! structure, we need to save this information to an internal state, for later 
  ! use in Array/Field creation.
  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isDataPtr => internalState%wrap

  isDataPtr%numOwnedNodes = 0
  isDataPtr%numForeignNodes = 0
  do ip=1, npa
      if (nodeowners(ip) == localPet) then 
      isDataPtr%numOwnedNodes = isDataPtr%numOwnedNodes + 1
    else 
      isDataPtr%numForeignNodes = isDataPtr%numForeignNodes + 1
    endif
  enddo

  if (isDataPtr%numForeignNodes + isDataPtr%numOwnedNodes /= npa) then 
    localrc = ESMF_RC_ARG_SIZE
    write(message, '(A,I4.4,A,I4.4,A,I4.4,A)') trim(compName)//' mesh with '// &
      'mismatching number of augmented npa=',npa,', owned=',isDataPtr%numOwnedNodes, &
      ' and foreign=', isDataPtr%numForeignNodes,' nodes'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  write(message,'(A,I6,A,I6,A,I6,A)') trim(compName)//' mesh with '// &
    'matching number of augmented npa=',npa,', owned=',isDataPtr%numOwnedNodes, &
    ' and foreign=', isDataPtr%numForeignNodes,' nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  
  allocate(isDataPtr%ownedNodeIds(isDataPtr%numOwnedNodes), stat=localrc)
  allocate(isDataPtr%foreignNodeIds(isDataPtr%numForeignNodes), stat=localrc)

  ownedCount = 0
  foreignCount = 0
  do ip=1,npa
    if (nodeowners(ip) == localPet) then
      ownedCount=ownedCount + 1
      isDataPtr%ownedNodeIds(ownedCount) = ip
    else
      foreignCount=foreignCount + 1
      isDataPtr%foreignNodeIds(foreignCount) = ip
    endif
  enddo

  if (isDataPtr%numForeignNodes /= foreignCount) then 
    localrc = ESMF_RC_ARG_SIZE
    write(message, '(A,I6,A,I6,A)') trim(compName)//' mesh with '// &
      'mismatching count of foreign nodes ',isDataPtr%numForeignNodes,' /= ', &
      foreignCount
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  if (isDataPtr%numOwnedNodes /= ownedCount) then 
    localrc = ESMF_RC_ARG_SIZE
    write(message, '(A,I6,A,I6,A)') trim(compName)//' mesh with '// &
      'mismatching count of owned nodes ',isDataPtr%numOwnedNodes,' /= ', &
      ownedCount
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  nvcount=0
  do i=1, nea
    elementids(i)=ielg(i)
    if(i34(i)==3) then
      elementtypes(i)=ESMF_MESHELEMTYPE_TRI
    else
      elementtypes(i)=ESMF_MESHELEMTYPE_QUAD
    endif
    elementmask(i)=idry_e(i)
    do ii=1,i34(i)
      elLocalNode(ii)=elnode(ii,i)
      nvcount = nvcount+1
      nv(nvcount) =elnode(ii,i) 
    end do

    elementcoords2d(2*i-1)=sum(nodecoords2d(2*elLocalNode(1:i34(i))-1))/i34(i)
    elementcoords2d(2*i)=sum(nodecoords2d(2*elLocalNode(1:i34(i))))/i34(i)
  end do ! i=1, nea

  if(ubound(nv,1) /= nvcount) then
    localrc=ESMF_RC_ARG_SIZE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc) 
  endif

  !> This is a three-part mesh generation, with later addition of node 
  !> and element information
  mesh2d = ESMF_MeshCreate(parametricDim=2, spatialdim=2, coordSys=coordsys, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_) 

  !> We can pass a nodalDistgrid (optional) 
  call ESMF_MeshAddNodes(mesh2d, nodeIds=nodeids, nodeCoords=nodecoords2d, &
    nodeOwners=nodeowners, nodeMask=nodemask, rc=localrc) 
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

!> @todo We need a generic handling of intformat
  write(message, '(A,I5,A,I5,A,I5,A,I5,A)') trim(compName)//' created mesh from ', &
   npa, '/', np, ' nodes with ', ownedCount, ' own and ', foreigncount, ' foreign'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  !> optional arguments are elementArea and elementDistgrid
  call ESMF_MeshAddElements(mesh2d, elementIds=elementids, elementTypes=elementtypes, &
    elementConn=nv, elementMask=elementmask, elementCoords=elementcoords2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_MeshGet(mesh2d, numOwnedNodes=npo, numOwnedElements=neo, &
    rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  
  if (isDataPtr%numOwnedNodes /= npo) then 
    localrc = ESMF_RC_ARG_SIZE
    write(message, '(A,I6,A,I6,A)') trim(compName)//' mesh with '// &
      'mismatching count of owned nodes ',isDataPtr%numOwnedNodes,' /= ', &
      npo
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  write(message, '(A,I5,A,I5,A,I5,A)') trim(compName)//' added to mesh from ', nea, &
    ' augmented and ', ne, ' resident elements ', neo, ' owned elements'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  call ESMF_GridCompSet(comp, mesh=mesh2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  write(message, '(A)') trim(compName)//' added mesh (node) to component'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  !> @todo the following steps don't work in the NUOPC cap yet

  !> Create fields for export to describe mesh (this information is not yet
  !> accessible with ESMF_MeshGet calls)
  !> @todo remove this part of the code once there is a suitable ESMF implementation

  !> Create a dummy field to satisfy ugrid conventions
  field = ESMF_FieldEmptyCreate(name='mesh_topology', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_AttributeSet(field, 'cf_role', 'mesh_topology', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_AttributeSet(field, 'topology_dimension', 2, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_AttributeSet(field, 'node_coordinates', &
   'mesh_node_lon mesh_node_lat', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_AttributeSet(field, 'face_node_connectivity', 'mesh_element_node_connectivity', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_GridCompGet(comp, exportState=exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !WIP here
  !call SCHISM_StateFieldCreate(comp, exportState, 'mesh_global_node_id', field, rc=localrc)
  !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  fieldName = 'mesh_global_node_id'
  field = ESMF_FieldCreate(mesh2d, name=fieldName,  &
    meshloc=ESMF_MESHLOC_NODE, typeKind=ESMF_TYPEKIND_I4, rc=localrc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtrI41, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  farrayPtrI41 = isDataPtr%ownedNodeIds

  call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_GridCompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  fieldName = 'mesh_global_element_id'
  field = ESMF_FieldCreate(mesh2d, name=fieldName,  &
    meshloc=ESMF_MESHLOC_ELEMENT, typeKind=ESMF_TYPEKIND_I4, rc=localrc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtrI41, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  farrayPtrI41 = elementIds(1:nea)

  call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on elements'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  nullify(farrayPtrI41)

  fieldName = 'mesh_element_node_connectivity'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
    meshloc=ESMF_MESHLOC_ELEMENT, ungriddedLBound=(/1/), ungriddedUBound=(/4/), &
    typeKind=ESMF_TYPEKIND_I4, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_FieldGet(field, farrayPtr=farrayPtrI42, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  do i=1,nea
    do n=1,i34(i)
      farrayPtrI42(i,n) = iplg(elnode(n,i))
    end do
  end do

  call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on elements'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  nullify(farrayPtrI42)

  ! clean up
  deallocate(nodeids, stat=localrc)
  deallocate(nodecoords2d, stat=localrc)
  deallocate(nodeowners, stat=localrc)
  deallocate(nodemask, stat=localrc)
  deallocate(elementids, stat=localrc)
  deallocate(elementmask, stat=localrc)
  deallocate(elementtypes, stat=localrc)
  if (allocated(elementCoords2d)) deallocate(elementCoords2d, stat=localrc)
  deallocate(nv, stat=localrc)

end subroutine SCHISM_MeshCreateNode

subroutine schism_esmf_add_bottom_tracer(name, mesh2d, tr_id, exportState, &
    importState,add_ws,rc)

  use schism_glbl, only: tr_el, kbe, wsett, rkind, nea, ne, nvrt
  implicit none

  character(len=ESMF_MAXSTR), intent(in)  :: name
  type(ESMF_Mesh), intent(in)      :: mesh2d
  integer, intent(in)              :: tr_id
  type(ESMF_State), intent(inout)  :: exportState
  type(ESMF_State), intent(inout), optional :: importState
  integer, intent(inout), optional :: rc
  logical, intent(in), optional    :: add_ws
  logical                          :: add_ws_
  real(kind=rkind), pointer        :: schism_ptr2d(:)
  integer                          :: localrc, rc_, i
  integer(ESMF_KIND_I4),dimension(1:1) :: maskValues=(/1/)
  type(ESMF_Field)                 :: field

  rc_ = ESMF_SUCCESS
  add_ws_ = .false.
  if (present(add_ws)) add_ws_ = add_ws

  field = ESMF_FieldCreate(mesh2d, name=trim(name)//'_at_soil_surface', &
                           typekind=ESMF_TYPEKIND_R8, &
                           meshloc=meshloc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  ! add maskValues to be used in regridding
  !call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !   initialize
  call ESMF_FieldGet(field,farrayPtr=schism_ptr2d,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  ! Consider all local elements from own and adjacent foreign nodes
  ! (nea) as these are the ones that are considered "owned" by ESMF
  do i=1,ne
    schism_ptr2d(i) = tr_el(tr_id,max(1,kbe(i)+1),i)
  end do

  call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (add_ws_) then
    field = ESMF_FieldCreate(mesh2d, name=trim(name)//'_z_velocity_at_soil_surface', &
                           typekind=ESMF_TYPEKIND_R8, &
                           meshloc=meshloc, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
    ! add maskValues to be used in regridding
    !call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    !   initialize
    call ESMF_FieldGet(field,farrayPtr=schism_ptr2d,rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    do i=1,ne
      ! to-do: requires to interpolate wsett to nodes
      schism_ptr2d(i) = wsett(tr_id,max(kbe(i)+1,1),i)
    end do

    call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  end if

  if (present(importState)) then
  ! add upward flux fields into importState
    field = ESMF_FieldCreate(mesh2d, name = trim(name)//'_upward_flux_at_soil_surface', &
                           typekind=ESMF_TYPEKIND_R8, &
                           meshloc=meshloc, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
    ! add maskValues to be used in regridding
    !call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    !   initialize
    call ESMF_FieldGet(field,farrayPtr=schism_ptr2d,rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
    schism_ptr2d = 0.0d0

    call ESMF_StateAddReplace(importState, (/field/), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  end if


  if (present(rc)) rc = rc_

end subroutine schism_esmf_add_bottom_tracer

subroutine schism_esmf_update_bottom_tracer(exportstate,name,tr_id,importState,dt,rc)
  use schism_glbl, only: nea, ne, kbe, tr_el, rkind, ze, idry_e
  implicit none

  type(ESMF_State), intent(inout)         :: exportState
  character(len=ESMF_MAXSTR), intent(in)  :: name
  integer, intent(in)                     :: tr_id
  type(ESMF_State), intent(inout),optional :: importState
  real(rkind), intent(in),optional        :: dt
  integer, intent(inout), optional        :: rc
  real(kind=ESMF_KIND_R8), pointer        :: ptr2d(:)
  integer                                 :: localrc, rc_, i
  type(ESMF_Field)                        :: field

  rc_ = ESMF_SUCCESS

  ! integrate bottom flux
  if (present(importState)) then
    call ESMF_StateGet(importState, trim(name)//'_upward_flux_at_soil_surface', field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    call ESMF_FieldGet(field, farrayPtr=ptr2d,rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    do i=1,ne
      if (idry_e(i)==0) then
        tr_el(tr_id,kbe(i)+1,i) = ptr2d(i)*dt/(ze(kbe(i)+1,i)-ze(kbe(i),i)) + tr_el(tr_id,max(kbe(i)+1,1),i)
      end if
    end do
  end if

  !   bottom temp, copy uncontiguous data into field in ESMF
  call ESMF_StateGet(exportState, trim(name)//'_at_soil_surface', field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_FieldGet(field, farrayPtr=ptr2d,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  do i=1,ne
    ptr2d(i) = tr_el(tr_id,max(kbe(i)+1,1),i)
  end do

  if (present(rc)) rc = rc_

end subroutine schism_esmf_update_bottom_tracer

!! Add a CIM component attribute package to a component
#undef ESMF_METHOD
#define ESMF_METHOD "addCIM"
subroutine addCIM(comp, rc)

  implicit none

  type(ESMF_GridComp)                          :: comp
  integer(ESMF_KIND_I4), intent(out), optional :: rc

  integer(ESMF_KIND_I4)      :: localrc, petCount
  type(ESMF_VM)              :: vm
  character(len=ESMF_MAXSTR) :: message, convention, purpose
  type(ESMF_State)           :: importState

  if (present(rc)) rc=ESMF_SUCCESS

  convention='CIM 1.5'
  purpose='ModelComp'

  !call ESMF_AttributeAdd(comp, convention=convention, &
  !  purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !call ESMF_AttributeSet(comp, 'ShortName', 'schism', &
  !  convention=convention, purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !call ESMF_AttributeSet(comp, 'LongName', 'schism', convention=convention, &
  !  purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !call ESMF_AttributeSet(comp, 'ModelType', 'ocean', convention=convention, &
  !  purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGet(comp, importState=importState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  call ESMF_AttributeGet(importState, name='simulation_start', value=message, defaultvalue='Untitled', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !call ESMF_AttributeSet(comp, 'SimulationStartDate', message, convention=convention, &
  !  purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  call ESMF_AttributeGet(importState, name='simulation_stop', value=message, defaultvalue='Untitled', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !call ESMF_AttributeSet(comp, 'SimulationDuration', message, convention=convention, &
  !  purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message,'(I4)') petCount
  !call ESMF_AttributeSet(comp, 'SimulationNumberOfProcessingElements', message, convention=convention, &
  !  purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  purpose='Platform'
  !call ESMF_AttributeGetAttPack(comp, convention, purpose, attpack=attpack, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !call ESMF_AttributeSet(comp, 'MachineName', 'unknown', convention=convention, &
  !  purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine addCIM

function clockCreateFrmParam(filename, kwe, relaxedFlag, rc) result(clock)

  use esmf
  implicit none

  character(len=ESMF_MAXSTR), intent(in)           :: filename
  type(ESMF_KeywordEnforcer), intent(in), optional :: kwe
  logical, intent(in), optional                    :: relaxedFlag
  integer(ESMF_KIND_I4), intent(out), optional     :: rc
  type(ESMF_Clock)                                 :: clock

  logical                 :: isPresent
  integer(ESMF_KIND_I4)   :: unit, localrc, rc_
  type(ESMF_Time)         :: stopTime, startTime
  type(ESMF_TimeInterval) :: timeStep
  character(len=ESMF_MAXSTR) :: message

  integer(ESMF_KIND_I4)  :: start_year=2000, start_month=1, start_day=1
  integer(ESMF_KIND_I4)  :: start_hour=0, rnday=2
  namelist /global/ start_year, start_month, start_day, start_hour, rnday

  if (present(rc)) rc = ESMF_SUCCESS

  inquire(file=filename, exist=isPresent)
  if (present(relaxedFlag)) then
    if (.not.relaxedFlag .and. .not.isPresent) then
      write(message, '(A)') '-- cannot find required '//trim(filename)
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)

      if (.not.present(rc)) then
        localrc = ESMF_RC_FILE_OPEN
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
      else
        rc = ESMF_RC_FILE_OPEN
        return
      endif
    endif
  endif

  if (isPresent) then

    call ESMF_UtilIOUnitGet(unit, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    open(unit, file=filename, iostat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    read(unit, nml=global, iostat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    close(unit)
  else
    write(message, '(A)') '-- assumes default clock 2000-01-01T00:00:00 + 2 days'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif

  write(message,'(I4,A,I2,A,I2,A,I2,A)') start_year,'-',start_month,'-', &
    start_day, 'T', start_hour,':00:00'

  ! Set day as timestep temporarily to count later to stop time
  !call ESMF_TimeSet(startTime, yy=start_year, mm=start_month, dd=start_day, &
  !  h=start_hour, rc=localrc)
  call ESMF_TimeSet(startTime, trim(message), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_TimeIntervalSet(timeStep, h=rnday, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  stopTime = startTime + timeStep

  ! Only now define the coupling timestep as fraction of full timeStep
  timeStep = timeStep / 24

  clock = ESMF_ClockCreate(timeStep, startTime, stopTime=stopTime, &
    name='main clock', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

end function clockCreateFrmParam

#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_FieldRealize"
subroutine SCHISM_FieldRealize(state, itemName, kwe, grid, mesh, typekind, rc)

  use NUOPC, only: NUOPC_Realize

  type(ESMF_State), intent(inout)    :: state
  character(len=*), intent(in)       :: itemName
  type(ESMF_KeywordEnforcer), intent(in), optional :: kwe
  type(ESMF_Grid), intent(in), optional :: grid
  type(ESMF_Mesh), intent(in), optional :: mesh
  type(ESMF_TypeKind_Flag), intent(in), optional :: typekind
  integer(ESMF_KIND_I4), intent(out), optional :: rc

  integer(ESMF_KIND_I4)       :: rc_, localrc
  character(len=ESMF_MAXSTR)  :: message
  type(ESMF_Field)            :: field
  type(ESMF_FieldStatus_Flag) :: fieldStatus
  type(ESMF_StateItem_Flag)   :: itemType

  rc_ = ESMF_SUCCESS

  if (present(grid).and.present(mesh)) then
    write(message, '(A)') '-- does not accept both mesh and grid'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    rc_ = ESMF_RC_ARG_WRONG
  elseif (.not.present(mesh).and..not.present(grid)) then
    write(message, '(A)') '-- needs either mesh or grid as argument'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    rc_ = ESMF_RC_ARG_WRONG
  endif

  if (rc_ /= ESMF_SUCCESS) then
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    if (present(rc)) then
      rc = rc_
      return
    endif
    localrc = rc_
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  call ESMF_StateGet(state, itemName=trim(itemName), itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (itemType /= ESMF_STATEITEM_FIELD) then 
    localrc = ESMF_RC_ARG_WRONG
    write(message, '(A)') '-- only accepts fields, obtained wrong type in item '// &
      trim(itemName)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    if (present(rc)) then
      rc = localrc
      return
    endif
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif
   
  call ESMF_StateGet(state, field=field, itemName=trim(itemName), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_FieldGet(field, status=fieldStatus, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (fieldStatus == ESMF_FIELDSTATUS_EMPTY .and. present(grid)) then
    call ESMF_FieldEmptySet(field, grid=grid, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  elseif (fieldStatus == ESMF_FIELDSTATUS_EMPTY .and. present(mesh)) then
    call ESMF_FieldEmptySet(field, mesh=mesh, meshloc=ESMF_MESHLOC_ELEMENT ,rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  call ESMF_FieldGet(field, status=fieldStatus, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (fieldStatus == ESMF_FIELDSTATUS_GRIDSET .and. present(typekind)) then
    call ESMF_FieldEmptyComplete(field, typekind=typekind, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  elseif (fieldStatus == ESMF_FIELDSTATUS_GRIDSET) then 
    call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  ! There is not need to formally call Realize() when completing the
  ! adverised field directly. However, calling Realize() also works.
  call NUOPC_Realize(state, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

end subroutine SCHISM_FieldRealize

subroutine schism_esmf_topbottom_tracer(name, mesh2d, tr_id, exportState, importState, &
    add_ws, rc)

  use schism_glbl, only: tr_el, tr_nd, kbe, wsett, rkind, npa, np, nea, ne, nvrt
  implicit none

  character(len=ESMF_MAXSTR), intent(in) :: name
  type(ESMF_Mesh), intent(in)            :: mesh2d
  integer, intent(in)                    :: tr_id
  type(ESMF_State), intent(inout)        :: exportState
  type(ESMF_State), intent(inout), optional :: importState
  integer, intent(inout), optional :: rc
  logical, intent(in), optional    :: add_ws
  logical                          :: add_ws_
  real(kind=rkind), pointer        :: schism_ptr2d(:)
  integer                          :: localrc, rc_, i
  integer(ESMF_KIND_I4),dimension(1:1) :: maskValues=(/1/)
  type(ESMF_Field)                 :: field

  rc_ = ESMF_SUCCESS
  add_ws_ = .false.
  if (present(add_ws)) add_ws_ = add_ws

  field = ESMF_FieldCreate(mesh2d, name=trim(name)//'_at_soil_surface', &
                           typekind=ESMF_TYPEKIND_R8, &
                           meshloc=meshloc, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  ! add maskValues to be used in regridding
  !> @todo re-enable when converted to eSMF_Info
  !call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  !   initialize
  call ESMF_FieldGet(field,farrayPtr=schism_ptr2d,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  ! Consider all local elements from own and adjacent foreign nodes
  ! (nea) as these are the ones that are considered "owned" by ESMF
  do i=1,ne
    schism_ptr2d(i) = tr_el(tr_id,max(1,kbe(i)+1),i)
  end do

  call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (add_ws_) then
    field = ESMF_FieldCreate(mesh2d, name=trim(name)//'_z_velocity_at_soil_surface', &
                           typekind=ESMF_TYPEKIND_R8, &
                           meshloc=meshloc, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
    ! add maskValues to be used in regridding
    !call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    !   initialize
    call ESMF_FieldGet(field,farrayPtr=schism_ptr2d,rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    do i=1,ne
      ! to-do: requires to interpolate wsett to nodes
      schism_ptr2d(i) = wsett(tr_id,max(kbe(i)+1,1),i)
    end do

    call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  end if

  if (present(importState)) then
  ! add upward flux fields into importState
    field = ESMF_FieldCreate(mesh2d, name = trim(name)//'_upward_flux_at_soil_surface', &
                           typekind=ESMF_TYPEKIND_R8, &
                           meshloc=meshloc, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    ! add maskValues to be used in regridding
    !call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    !   initialize
    call ESMF_FieldGet(field,farrayPtr=schism_ptr2d,rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
    schism_ptr2d = 0.0d0

    call ESMF_StateAddReplace(importState, (/field/), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  end if


  if (present(rc)) rc = rc_

end subroutine schism_esmf_topbottom_tracer


#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateImportWaveTensor"
subroutine SCHISM_StateImportWaveTensor(state, isPtr, rc)

  use schism_glbl, only: np
  implicit none

  type(ESMF_State), intent(in)                  :: state
  type(type_InternalState), pointer, intent(in) :: isPtr
  integer(ESMF_KIND_I4), intent(out), optional  :: rc

  logical                    :: isPresent
  integer(ESMF_KIND_I4)      :: localrc, rc_, i
  character(len=ESMF_MAXSTR) :: message
  type(ESMF_Field)           :: field
  type(ESMF_StateItem_Flag)  :: itemType

  real(ESMF_KIND_R8), pointer :: farrayPtr1(:) => null()
  real(ESMF_KIND_R8), pointer :: eastward_wave_radiation_stress(:) => null()
  real(ESMF_KIND_R8), pointer :: eastward_northward_wave_radiation_stress(:) => null()
  real(ESMF_KIND_R8), pointer :: northward_wave_radiation_stress(:) => null()

  if (present(rc)) rc=ESMF_SUCCESS
  localrc = ESMF_SUCCESS

  allocate(farrayPtr1(np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateGet(state, itemname='eastward_wave_radiation_stress', itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (itemType /= ESMF_STATEITEM_FIELD) then
    write(message,'(A)') '--- skipped computing of wave stress'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    return
  endif

  if (itemType == ESMF_STATEITEM_FIELD) then 
    call ESMF_StateGet(state, itemname='eastward_wave_radiation_stress', field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call SCHISM_FieldPtrUpdate(field, farrayPtr1, isPtr=isPtr, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    allocate(eastward_wave_radiation_stress(isPtr%numOwnedNodes), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    eastward_wave_radiation_stress(:) = farrayPtr1(1:isPtr%numOwnedNodes)
  
  endif 

  call ESMF_StateGet(state, itemname='eastward_northward_wave_radiation_stress', itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (itemType == ESMF_STATEITEM_FIELD) then 
    call ESMF_StateGet(state, itemname='eastward_northward_wave_radiation_stress', field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call SCHISM_FieldPtrUpdate(field, farrayPtr1, isPtr=isPtr, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    allocate(eastward_northward_wave_radiation_stress(isPtr%numOwnedNodes), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    eastward_northward_wave_radiation_stress(:) = farrayPtr1(1:isPtr%numOwnedNodes)
  endif 

  call ESMF_StateGet(state, itemname='northward_wave_radiation_stress', itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (itemType == ESMF_STATEITEM_FIELD) then 
    call ESMF_StateGet(state, itemname='northward_wave_radiation_stress', field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call SCHISM_FieldPtrUpdate(field, farrayPtr1, isPtr=isPtr, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    allocate(northward_wave_radiation_stress(isPtr%numOwnedNodes), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    northward_wave_radiation_stress(:) = farrayPtr1(1:isPtr%numOwnedNodes)
  endif 

  call compute_wave_force_lon(eastward_wave_radiation_stress, & 
    eastward_northward_wave_radiation_stress,northward_wave_radiation_stress)

  deallocate(eastward_wave_radiation_stress)
  deallocate(eastward_northward_wave_radiation_stress)
  deallocate(northward_wave_radiation_stress)
  deallocate(farrayPtr1)
  
end subroutine SCHISM_StateImportWaveTensor

end module schism_esmf_util
