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
! @license Apache License, Version 2.0 (the "License");
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

  type type_InternalState
    sequence ! why is this needed here? taken from documentation
    ! Store the number of and indices in the 1:np resident nodes
    integer(ESMF_KIND_I4)          :: numForeignNodes=0
    integer(ESMF_KIND_I4)          :: numOwnedNodes=0
    integer(ESMF_KIND_I4), pointer :: ownedNodeIds(:) => null()
    integer(ESMF_KIND_I4), pointer :: foreignNodeIds(:) => null()
    type(ESMF_RouteHandle)         :: haloHandle
  end type

  type type_InternalStateWrapper
    sequence ! why is this needed here? taken from documentation
    type(type_InternalState), pointer :: wrap => null()
  end type

!  public addSchismMesh, clockCreateFrmParam, SCHISM_FieldRealize
  public clockCreateFrmParam, SCHISM_FieldRealize
  public type_InternalState, type_InternalStateWrapper
  public SCHISM_StateFieldCreateRealize, SCHISM_StateGetField, SCHISM_FieldPtrUpdate
  private

contains

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateGetField"
subroutine SCHISM_StateGetField(state, itemName, kwe, farrayPtr, field,  rc)

  type(ESMF_State), intent(in)                        :: state
  character(len=*), intent(in)                        :: itemName
  type(ESMF_KeywordEnforcer), intent(in), optional    :: kwe
  type(ESMF_Field), intent(out), optional             :: field
  real(ESMF_KIND_R8), pointer, intent(out), optional  :: farrayPtr(:)
  integer(ESMF_KIND_I4), intent(out), optional        :: rc

  integer(ESMF_KIND_I4)          :: rc_, localrc
  type(ESMF_Field)               :: field_
  character(len=ESMF_MAXSTR)     :: message
  type(ESMF_StateItem_Flag)      :: itemType

  localrc = ESMF_SUCCESS 
  if (present(rc)) rc = ESMF_SUCCESS
  if (.not.(present(field).or.(present(farrayPtr)))) return

  call ESMF_StateGet(state, itemname=trim(itemName), itemType=itemType, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (itemType /= ESMF_STATEITEM_FIELD) then 
    write(message,'(A)') '--- field '//trim(itemName)//' could not be found in state'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    localrc = ESMF_RC_NOT_FOUND
    if (present(rc)) rc = localrc
    return
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif
  
  call ESMF_StateGet(state, itemname=trim(itemName), field=field_, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (present(field)) field=field_
  if (.not.present(farrayPtr)) return

  call ESMF_FieldGet(field_, farrayptr=farrayPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine SCHISM_StateGetField

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_FieldPtrUpdate"
subroutine SCHISM_FieldPtrUpdate(field, farray, kwe, isPtr, rc)

  type(ESMF_Field), intent(inout)                   :: field
  real(ESMF_KIND_R8),  intent(inout), target        :: farray(:)
  type(ESMF_KeywordEnforcer), intent(in), optional  :: kwe
  type(type_InternalState), pointer, intent(in)     :: isPtr
  
  integer(ESMF_KIND_I4), intent(out), optional      :: rc

  real(ESMF_KIND_R8), pointer    :: farrayPtr1(:)
  integer(ESMF_KIND_I4)          :: rc_, localrc, ip
  character(len=ESMF_MAXSTR)     :: message
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

  isPresent = ESMF_RouteHandleIsCreated(isPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (isPresent) then 
    call ESMF_FieldHalo(field, routehandle=isPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif    

  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  
  do ip = 1, isPtr%numOwnedNodes
    farray(isPtr%ownedNodeIds(ip)) = farrayPtr1(ip)
  end do
  do ip = 1,isPtr%numForeignNodes
    farray(isPtr%foreignNodeIds(ip)) = farrayPtr1(ip+isPtr%numOwnedNodes)
  end do

end subroutine SCHISM_FieldPtrUpdate

!> @todo separate into ESMF and NUOPC parts that reside in different source files
#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateFieldCreateRealize"
subroutine SCHISM_StateFieldCreateRealize(comp, state, name, field, kwe, rc)

  use schism_glbl, only: nws 
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

  call ESMF_MeshGet(mesh, nodalDistgrid=distgrid, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isDataPtr => internalState%wrap
 
  if (isDataPtr%numForeignNodes > 0 .and. associated(isDataPtr%foreignNodeIds)) then 
    array = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, &
    haloSeqIndexList=isDataPtr%foreignNodeIds, &
    name=trim(name), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  else
    array = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, &
    name=trim(name), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  field = ESMF_FieldCreate(name=trim(name), mesh=mesh, array=array, &
    meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  isPresent = ESMF_RouteHandleIsCreated(isDataPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  if (.not. isPresent) then 
    call ESMF_FieldHaloStore(field, routehandle=isDataPtr%haloHandle, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
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

  call NUOPC_Realize(state, field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldHalo(field, routehandle=isDataPtr%haloHandle, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (NUOPC_IsConnected(field, rc=localrc) .and. nws /= 3) then
    write(message, '(A,I1,A)') trim(compName)//' connected field '//trim(name)// &
      '  not used with nws=', nws ,' (needs nws = 3)'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif

  if (present(rc)) rc = rc_

end subroutine SCHISM_StateFieldCreateRealize

#undef  ESMF_METHOD
#define ESMF_METHOD "addSchismMesh"
subroutine addSchismMesh(comp, rc)
! Define ESMF domain partition
  !> @todo apply only filter to 'use schism_glbl'
  use schism_glbl, only: pi, llist_type, elnode, i34, ipgl
  use schism_glbl, only: iplg, ielg, idry_e, idry, ynd, xnd
  use schism_glbl, only: ylat, xlon, npa, np, nea, ne, ics
  use schism_glbl, only:  nvrt

  implicit none

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Mesh)          :: mesh2d
  type(ESMF_DistGrid)      :: elementDistgrid, distgrid
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
  integer               :: ii,ip,ie, localrc
  integer               :: mynp,myne,rank2
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

  write(0,*)'__LINE__ inside addSchismMesh'
  call ESMF_Finalize() 

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, name=compName, localPet=localPet, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGetInternalState(comp, internalState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  isDataPtr => internalState%wrap

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
!  numNodeHaloIdx=0
!  npa=npa

!  allocate(localNodes(npa), stat=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!  do i=1,npa
!    localNodes(i)=i
!  end do

  ! get inverse matching
!  allocate(schismToLocalNodes(npa))
!  schismToLocalNodes(:) = -1 ! initialize to -1
!  do i=1,npa
!    schismToLocalNodes(localNodes(i))=i
!  end do

  ! define mesh
  !Points to global node #
  allocate(nodeids(np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !Global coordinates
  allocate(nodecoords2d(2*np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !A node is owned by same rank across PETs; interface nodes are owned by min rank
  allocate(nodeowners(np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(nodemask(np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !Points to global elem #
  allocate(elementids(ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(elementtypes(ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(elementmask(ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(elementcoords2d(2*ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! nv (elemConn): 1D array for connectivity (packed from 2D array elnode).
  ! Outputs local node # 
  allocate(nv(sum(i34(1:ne))), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! set ESMF coordSys type
  if (ics==2) then
    coordsys=ESMF_COORDSYS_SPH_DEG
  else
    write(message, '(A)') trim(compName)//' uses a cartesian coordinate system'// &
       ' which may not be suitable for coupling'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    coordsys=ESMF_COORDSYS_CART
  endif

  do ip=1, np
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
    rank2=ipgl(iplg(ip))%rank
    nodeowners(ip) =rank2 !init 
    if(associated(ipgl(iplg(ip))%next)) then !interface or ghost node
      if(ipgl(iplg(ip))%next%rank<rank2) then
        nodeowners(ip) =ipgl(iplg(ip))%next%rank
      endif
    endif

!    if (ip<=np) then
!      ! if iplg() is resident, ipgl(iplg(i))%rank=myrank
!      nodeowners(ip) = ipgl(iplg(ip))%rank
!    else !not executed at the moment
!      ! get owner of foreign node, the SCHISM manual tells us to not use
!      ! ipgl for foreign nodes ???
!      ! ipgl%next%next%next.... is the linked list, with ranks in ascending order.  Unless ipgb is an interface node (i.e., resident in more than 1 process), the list has only 1 entry
!      nextp => ipgl(iplg(ip))
!      ! We here advance to the end of the list, Joseph suggested to just take the next element
!      do while (associated(nextp%next))
!        nextp => nextp%next
!      end do
!      nodeowners(ip) = nextp%rank
!    end if
    nodemask(ip)         = idry(ip)
  end do !ip

  isDataPtr%numOwnedNodes = 0
  isDataPtr%numForeignNodes = 0
  do ip=1,np
    if (nodeowners(ip) == localPet) then 
      isDataPtr%numOwnedNodes = isDataPtr%numOwnedNodes + 1
    else 
      isDataPtr%numForeignNodes = isDataPtr%numForeignNodes + 1
    endif
  enddo

  if (isDataPtr%numForeignNodes + isDataPtr%numOwnedNodes /= np) then 
    localrc = ESMF_RC_ARG_SIZE
    write(message, '(A,I4.4,A,I4.4,A,I4.4,A)') trim(compName)//' mesh with '// &
      'mismatching number of resident np=',np,', owned=',isDataPtr%numOwnedNodes, &
      ' and foreign=', isDataPtr%numForeignNodes,' nodes'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  allocate(isDataPtr%ownedNodeIds(isDataPtr%numOwnedNodes), stat=localrc)
  allocate(isDataPtr%foreignNodeIds(isDataPtr%numForeignNodes), stat=localrc)

  ownedCount = 0
  foreignCount = 0
  do i=1,np
    if (nodeowners(ip) == localPet) then
      ownedCount=ownedCount + 1
      isDataPtr%ownedNodeIds(ownedCount) = ip
    else
      foreignCount=foreignCount + 1
      isDataPtr%foreignNodeIds(foreignCount) = ip
    endif
  enddo

  write(message, '(A,I4.4,A,I4.4,A,I4.4,A)') trim(compName)//' mesh with '// &
    'number of resident np=',np,', owned=',isDataPtr%numOwnedNodes, &
    ' and foreign=', isDataPtr%numForeignNodes,' nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  write(*,*) 'Owned nodes on PET ',localPet,isDataPtr%ownedNodeIds
  write(*,*) 'Foreign nodes on PET ',localPet,isDataPtr%foreignNodeIds

  nvcount=0
  do i=1,ne
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

    ! Element coord is the sum of nodes divided by number of nodes
    elementcoords2d(2*i-1)=sum(nodecoords2d(2*elLocalNode(1:i34(i))-1))/i34(i)
    elementcoords2d(2*i)=sum(nodecoords2d(2*elLocalNode(1:i34(i))))/i34(i)
  end do !i

  if(ubound(nv,1)/=nvcount) then
    localrc=ESMF_RC_ARG_SIZE
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc) 
  endif

  ! create element distgrid (distribute)
  elementDistgrid = ESMF_DistgridCreate(elementids,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  mesh2d = ESMF_MeshCreate(parametricDim=2,spatialdim=2,nodeIds=nodeids, &
             nodeCoords=nodecoords2d,nodeOwners=nodeowners, &
             coordSys=coordsys, &
             nodeMask=nodemask,elementMask=elementmask, &
             elementIds=elementids, elementTypes=elementtypes, &
             elementCoords=elementcoords2d, &
             elementDistgrid=elementDistgrid, &
             elementConn=nv, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_MeshGet(mesh2d, numOwnedNodes=mynp, numOwnedElements=myne, elementDistgrid=distgrid, &
    rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(testids(myne))
  call ESMF_DistGridGet(distgrid,localDE=0,seqIndexList=testids,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,I3.3,A,I3.3,A)') trim(compName)//' created mesh from "', np, &
    ' resident nodes and ', myne, ' resident elements in SCHISM'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)

  write(message, '(A,I3.3,A,I3.3,A)') trim(compName)//' created mesh with "', mynp, &
    'owned nodes and ', myne, ' owned elements'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)

  !> @todo the following might overflow the message buffer easily ...
  write(message,*) 'elementIds:',elementIds
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)

  write(message,*) 'distgridElementIds:',testids
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)

  deallocate(testids)

  call ESMF_GridCompSet(comp, mesh=mesh2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo the following steps don't work in the NUOPC cap yet

  !> Create fields for export to describe mesh (this information is not yet
  !> accessible with ESMF_MeshGet calls)
  !> @todo remove this part of the code once there is a suitable ESMF implementation

  !> Create a dummy field to satisfy ugrid conventions
  field = ESMF_FieldEmptyCreate(name='mesh_topology', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(field, 'cf_role', 'mesh_topology', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(field, 'topology_dimension', 2, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(field, 'node_coordinates', &
   'mesh_node_lon mesh_node_lat', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(field, 'face_node_connectivity', 'mesh_element_node_connectivity', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGet(comp, exportState=exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)


  call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  fieldName = 'mesh_global_node_id'
  field = ESMF_FieldCreate(mesh2d, name=fieldName,  &
    meshloc=ESMF_MESHLOC_NODE, typeKind=ESMF_TYPEKIND_I4, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtrI41, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  farrayPtrI41 = nodeIds(1:np)

  call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  fieldName = 'mesh_global_element_id'
  field = ESMF_FieldCreate(mesh2d, name=fieldName,  &
    meshloc=ESMF_MESHLOC_ELEMENT, typeKind=ESMF_TYPEKIND_I4, rc=localrc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtrI41, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  farrayPtrI41 = elementIds(1:ne)

  call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on elements'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  nullify(farrayPtrI41)

  fieldName = 'mesh_element_node_connectivity'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
    meshloc=ESMF_MESHLOC_ELEMENT, ungriddedLBound=(/1/), ungriddedUBound=(/4/), &
    typeKind=ESMF_TYPEKIND_I4, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtrI42, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i=1,ne
    do n=1,i34(i)
      farrayPtrI42(i,n) = iplg(elnode(n,i))
    end do
  end do

  call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

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

  write(message, '(A)') trim(compName)//' created 2D mesh"'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

end subroutine addSchismMesh

subroutine schism_esmf_add_bottom_tracer(name,mesh2d,tr_id,exportState, &
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
                           meshloc=ESMF_MESHLOC_ELEMENT, rc=localrc)
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
                           meshloc=ESMF_MESHLOC_ELEMENT, rc=localrc)
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
                           meshloc=ESMF_MESHLOC_ELEMENT, rc=localrc)
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
    call ESMF_FieldEmptySet(field, mesh=mesh, rc=localrc)
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
                           meshloc=ESMF_MESHLOC_ELEMENT, rc=localrc)
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
                           meshloc=ESMF_MESHLOC_ELEMENT, rc=localrc)
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
                           meshloc=ESMF_MESHLOC_ELEMENT, rc=localrc)
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




end module schism_esmf_util
