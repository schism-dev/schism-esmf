! This code is part of the SCHISM-ESMF interface, it defines utility
! functions used both by the NUOPC and ESMF caps
!
! @copyright (C) 2018, 2019, 2020 Helmholtz-Zentrum Geesthacht
! @author Carsten Lemmen carsten.lemmen@hzg.de
! @author Richard Hofmeister richard.hofmeister@hzg.de
!
! @license under the Apache License, Version 2.0 (the "License");
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
! This version accounts for halo(ghost) zone, because ESMF by default
! partitions among nodes instead of elements

#define ESMF_CONTEXT  line=__LINE__,file=ESMF_FILENAME,method=ESMF_METHOD
#define ESMF_ERR_PASSTHRU msg="SCHISM subroutine call returned error"
#undef ESMF_FILENAME
#define ESMF_FILENAME "schism_esmf_util.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism_esmf_util

  use esmf
  implicit none

contains

#undef  ESMF_METHOD
#define ESMF_METHOD "addSchismMesh"
subroutine addSchismMesh(comp, rc)

  !> @todo apply only filter to 'use schism_glbl'
  use schism_glbl, only: pi, llist_type, elnode, i34, ipgl
  use schism_glbl, only: iplg, ielg, idry_e, idry, ynd, xnd
  use schism_glbl, only: ylat, xlon, npa, np, nea, ne, lreadll
  use schism_glbl, only:  nvrt

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Mesh)       :: mesh2d, mesh3d
  type(ESMF_DistGrid)   :: elementDistgrid, distgrid
  type(ESMF_CoordSys_Flag) :: coordsys
  integer, dimension(:), allocatable            :: nodeids, elementids, nv
  real(ESMF_KIND_R8), dimension(:), allocatable :: nodecoords2d, nodecoords3d
  real(ESMF_KIND_R8), dimension(:), allocatable :: elementcoords2d, elementcoords3d
  integer, dimension(:), allocatable            :: nodeowners, elementtypes
  integer, dimension(:), allocatable            :: nodemask, elementmask
  integer, dimension(:), allocatable            :: tmpIdx, tmpIdx2, localNodes, nodeHaloIdx
  integer, dimension(:), allocatable            :: schismTolocalNodes
  integer, dimension(1:4)                       :: elLocalNode
  integer               :: numLocalNodes, numNodeHaloIdx
  integer               :: i,n,nvcount
  integer               :: ii,ip,ie, localrc
  integer               :: mynp,myne
  type(llist_type),pointer :: nextp=>null()

  integer(ESMF_KIND_I4), pointer, dimension(:)  :: farrayPtrI41 => null()
  integer(ESMF_KIND_I4), pointer, dimension(:,:):: farrayPtrI42 => null()

  real(ESMF_KIND_R8), parameter :: rad2deg=180.0d0/pi

  integer(ESMF_KIND_I4),dimension(1:1) :: maskValues=(/1/)
  character(len=ESMF_MAXSTR)           :: compName, message, fieldName
  type(ESMF_State)                     :: exportState
  type(ESMF_Field)                     :: field

  rc = ESMF_SUCCESS

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
  numNodeHaloIdx=0
  numLocalNodes=np

  allocate(tmpIdx(npa-np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(tmpIdx2(npa), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i=1,np
    tmpIdx2(i)=i
  end do
  do i=np+1,npa-np
    ! check existence in local elements
    do ie=1,ne
      if (any(elnode(1:i34(ie),ie)==i)) then
        numLocalNodes = numLocalNodes+1
        tmpIdx2(numLocalNodes)=i
      else
        numNodeHaloIdx = numNodeHaloIdx+1
        tmpIdx(numNodeHaloIdx)=i
      end if
    end do
  end do
  allocate(nodeHaloIdx(numNodeHaloIdx))
  nodeHaloIdx(:)=tmpIdx(1:numNodeHaloIdx)
  deallocate(tmpIdx)
  allocate(localNodes(numLocalNodes))
  localNodes(:)=tmpIdx2(1:numLocalNodes)
  deallocate(tmpIdx2)
  ! get inverse matching
  allocate(schismToLocalNodes(npa))
  schismToLocalNodes(:) = -1 ! initialize to -1
  do i=1,numLocalNodes
    schismToLocalNodes(localNodes(i))=i
  end do

  ! define mesh
  allocate(nodeids(numLocalNodes), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(nodecoords2d(2*numLocalNodes), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(nodecoords3d(3*numLocalNodes), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(nodeowners(numLocalNodes), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(nodemask(numLocalNodes), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(elementids(ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(elementtypes(ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(elementmask(ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(elementcoords2d(2*ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !allocate(elementcoords2d(3*nea), stat=localrc)
  allocate(nv(4*ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! set ESMF coordSys type
  if (lreadll) then
    coordsys=ESMF_COORDSYS_SPH_DEG
  else
    coordsys=ESMF_COORDSYS_CART
  endif

  do ip=1,numLocalNodes
    i = localNodes(ip)
    nodeids(ip)=iplg(i)
    if (lreadll) then
      ! if geographical coordinates present
      nodecoords2d(2*ip-1) = rad2deg*xlon(i)
      nodecoords2d(2*ip)   = rad2deg*ylat(i)
    else
      ! use cartesian coordinates
      nodecoords2d(2*ip-1) = xnd(i)
      nodecoords2d(2*ip)   = ynd(i)
    end if
    if (i<=np) then
      nodeowners(ip)       = ipgl(iplg(i))%rank
    else
      ! get owner of foreign node
      nextp => ipgl(iplg(i))
      do while (associated(nextp%next))
        nextp => nextp%next
      end do
      nodeowners(ip) = nextp%rank
    end if
    nodemask(ip)         = idry(i)
  end do

  nvcount=1

  do i=1,ne
    elementids(i)=ielg(i)
    elementtypes(i)=i34(i)
    elementmask(i)=idry_e(i)
    elLocalNode(:)=-1 ! get local elnode
    do ii=1,i34(i)
      elLocalNode(ii)=schismToLocalNodes(elnode(ii,i))
      nv(nvcount+ii-1) = schismTolocalNodes(elnode(ii,i))
    end do
    nvcount = nvcount+i34(i)

    ! Element coord is the sum of nodes divided by number of nodes
    elementcoords2d(2*i-1)=sum(nodecoords2d(2*elLocalNode(1:i34(i))-1))/i34(i)
    elementcoords2d(2*i)=sum(nodecoords2d(2*elLocalNode(1:i34(i))))/i34(i)
  end do

#if 0
  write(0,*) 'Local Nodes, np=',np
  do i=1,numLocalNodes
   write(0,*) 'localid, globalid, nodeowner',i,nodeids(i),nodeowners(i)
  end do
  nvcount=1
  do i=1,ne
    write(0,*) 'ne, conn',i,nv(nvcount:nvcount+i34(i)-1)
    nvcount = nvcount+i34(i)
  end do
#endif

  ! create element distgrid
  elementDistgrid = ESMF_DistgridCreate(elementids,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  mesh2d = ESMF_MeshCreate(parametricDim=2,spatialdim=2,nodeIds=nodeids, &
    nodeCoords=nodecoords2d,nodeOwners=nodeowners, &
    coordSys=coordsys, &
    nodeMask=nodemask,elementMask=elementmask, &
    elementIds=elementids, elementTypes=elementtypes, &
    elementCoords=elementcoords2d, &
    elementDistgrid=elementDistgrid, &
    elementConn=nv(1:nvcount-1), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

#if 0
  ! output mesh information from schism and esmf
  call ESMF_MeshGet(mesh2d, numOwnedNodes=mynp,numOwnedElements=myne, elementDistgrid=distgrid, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(testids(myne))
  call ESMF_DistGridGet(distgrid,localDE=0,seqIndexList=testids,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(0,*) 'esmf   owned nodes,elements',mynp,myne
  write(0,*) 'schism owned nodes,elements',np,ne
  write(0,*) 'schism   all nodes,elements',npa,nea
  write(0,*) 'elementIds',elementIds
  write(0,*) 'distgridElementIds',testids
  deallocate(testids)
#endif

  call ESMF_GridCompSet(comp, mesh=mesh2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  return
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

  !call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
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

  !call ESMF_StateAddReplace(exportstate, (/field/), rc=localrc)
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
  deallocate(LocalNodes, stat=localrc)
  deallocate(schismToLocalNodes, stat=localrc)

end subroutine addSchismMesh

subroutine schism_esmf_add_bottom_tracer(name,mesh2d,tr_id,exportState,importState,add_ws,rc)
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
  call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
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
    call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
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
    call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
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

  call ESMF_AttributeAdd(comp, convention=convention, &
    purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(comp, 'ShortName', 'schism', &
    convention=convention, purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(comp, 'LongName', 'schism', convention=convention, &
    purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(comp, 'ModelType', 'ocean', convention=convention, &
    purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGet(comp, importState=importState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeGet(importState, name='simulation_start', value=message, defaultvalue='Untitled', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(comp, 'SimulationStartDate', message, convention=convention, &
    purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeGet(importState, name='simulation_stop', value=message, defaultvalue='Untitled', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(comp, 'SimulationDuration', message, convention=convention, &
    purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message,'(I4)') petCount
  call ESMF_AttributeSet(comp, 'SimulationNumberOfProcessingElements', message, convention=convention, &
    purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  purpose='Platform'
  !call ESMF_AttributeGetAttPack(comp, convention, purpose, attpack=attpack, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeSet(comp, 'MachineName', 'unknown', convention=convention, &
    purpose=purpose, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine addCIM

end module schism_esmf_util
