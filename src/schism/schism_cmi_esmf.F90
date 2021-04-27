! This code is part of the SCHISM-ESMF interface
!
! @copyright (C) 2021 Helmholtz-Zentrum Hereon
! @copyright (C) 2018--2021 Helmholtz-Zentrum Geesthacht
!
! @author Carsten Lemmen <carsten.lemmen@hereon.de>
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
! This version accounts for halo(ghost) zone, because ESMF by default
! partitions among nodes instead of elements

#define ESMF_CONTEXT  line=__LINE__,file=ESMF_FILENAME,method=ESMF_METHOD
#define ESMF_ERR_PASSTHRU msg="SCHISM subroutine call returned error"
#undef ESMF_FILENAME
#define ESMF_FILENAME "schism_cmi_esmf.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism_cmi_esmf

  use schism_bmi
  use schism_esmf_util
  use esmf

  implicit none

  public SetServices

  !Saved cohort arrays for schism_get and _save, and all that require resetting
  !among different cohorts
  integer, save, allocatable :: idry_e_c(:,:),idry_s_c(:,:),idry_c(:,:), istack(:)
  integer, save, allocatable :: ncid_c(:) !CWB2021-1
  real(8), save, allocatable ::we_c(:,:,:),tr_el_c(:,:,:,:), su2_c(:,:,:),sv2_c(:,:,:), &
 &eta2_c(:,:),tr_nd_c(:,:,:,:),tr_nd0_c(:,:,:,:),q2_c(:,:,:),xl_c(:,:,:),dfv_c(:,:,:), &
 &dfh_c(:,:,:),dfq1_c(:,:,:),dfq2_c(:,:,:), uu2_c(:,:,:),vv2_c(:,:,:)


contains

#undef  ESMF_METHOD
#define ESMF_METHOD "SetServices"
subroutine SetServices(comp, rc)

  type(ESMF_GridComp)                :: comp
  integer(ESMF_KIND_I4), intent(out) :: rc

  integer(ESMF_KIND_I4)              :: localrc
! integer                            :: pdaf_stat

  rc = ESMF_SUCCESS

  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phase=0, userRoutine=InitializeP0, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phase=1, userRoutine=InitializeP1, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
    userRoutine=Run, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
    userRoutine=Finalize, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine SetServices

#undef  ESMF_METHOD
#define ESMF_METHOD "InitializeP0"
  subroutine InitializeP0(gridComp, importState, exportState, &
    parentClock, rc)

    implicit none

    type(ESMF_GridComp)         :: gridComp
    type(ESMF_State)            :: importState
    type(ESMF_State)            :: exportState
    type(ESMF_Clock)            :: parentClock
    integer, intent(out)        :: rc

    character(len=10)           :: InitializePhaseMap(1)
    character(len=ESMF_MAXSTR)  :: myName
    type(ESMF_Time)             :: currTime
    integer                     :: localrc

    rc=ESMF_SUCCESS

    InitializePhaseMap(1) = "IPDv00p1=1"

    call ESMF_AttributeAdd(gridComp, convention="NUOPC", &
      purpose="General", &
      attrList=(/"InitializePhaseMap"/), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_AttributeSet(gridComp, name="InitializePhaseMap", valueList=InitializePhaseMap, &
      convention="NUOPC", purpose="General", rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end subroutine InitializeP0

#undef  ESMF_METHOD
#define ESMF_METHOD "InitializeP1"
subroutine InitializeP1(comp, importState, exportState, clock, rc)

  !> @todo apply only filter to 'use schism_glbl', in the mid-term
  !> much of this code should go to schism_esmf_util and schism_bmi
  use schism_glbl, only: pi, llist_type, elnode, i34, ipgl
  use schism_glbl, only: iplg, ielg, idry_e, idry, ynd, xnd
  use schism_glbl, only: ylat, xlon, npa, np, nea, ne, ics
  use schism_glbl, only: windx2, windy2, pr2, airt2, shum2
  use schism_glbl, only: srad, fluxevp, fluxprc, tr_nd, uu2
  use schism_glbl, only: dt, rnday, vv2, nvrt,ifile,nc_out
!  use schism_msgp, only: schism_mpi_comm=>comm
  use schism_msgp, only: parallel_init
  use schism_io, only: fill_nc_header!ncid
! use mod_assimilation, only: outf !PDAF module
#ifdef USE_FABM
  use fabm_schism, only: fabm_istart=>istart, fs
#endif

  type(ESMF_GridComp)   :: comp
  type(ESMF_State)      :: importState
  type(ESMF_State)      :: exportState
  type(ESMF_Clock)      :: clock
  integer, intent(out)  :: rc

  type(ESMF_Field)      :: field
  type(ESMF_Clock)      :: schismClock
  type(ESMF_VM)         :: vm
  type(ESMF_Mesh)       :: mesh2d,mesh3d
  type(ESMF_DistGrid)   :: elementDistgrid,distgrid
  type(ESMF_CoordSys_Flag) :: coordsys
  integer, dimension(:), allocatable            :: nodeids,elementids,nv
  real(ESMF_KIND_R8), dimension(:), allocatable :: nodecoords2d, nodecoords3d
  real(ESMF_KIND_R8), dimension(:), allocatable :: elementcoords2d, elementcoords3d
  real(ESMF_KIND_R8), dimension(:), pointer     :: schism_windx,schism_windy
  real(ESMF_KIND_R8), dimension(:), pointer     :: schism_ptr2d
  integer, dimension(:), allocatable            :: nodeowners, elementtypes
  integer, dimension(:), allocatable            :: nodemask, elementmask
  integer, dimension(:), allocatable            :: tmpIdx, tmpIdx2, localNodes, nodeHaloIdx
  integer, dimension(:), allocatable            :: schismTolocalNodes
  integer, dimension(1:4)                       :: elLocalNode
  integer               :: numLocalNodes, numNodeHaloIdx
  integer               :: mpi_comm
  integer               :: i, n, nvcount
  integer               :: ii,ip,ie, iths=0, ntime=0
  integer               :: mynp,myne
  integer, dimension(:), allocatable :: testids
  real(ESMF_KIND_R8), parameter :: rad2deg=180.0d0/pi
  type(ESMF_TimeInterval) :: schism_dt
  type(llist_type),pointer :: nextp=>null()
  integer(ESMF_KIND_I4),dimension(1:1) :: maskValues=(/1/)
  integer(ESMF_KIND_I4), pointer, dimension(:)  :: farrayPtrI41 => null()
  integer(ESMF_KIND_I4), pointer, dimension(:,:):: farrayPtrI42 => null()

  character(len=ESMF_MAXSTR)  :: message, name, compName, fieldName
  integer(ESMF_KIND_I4)       :: localrc, petCount,localPet,cohortIndex,ncohort,runhours,schism_dt2
  logical                     :: isPresent
  character(len=ESMF_MAXSTR)  :: configFileName, simulationDirectory
  character(len=ESMF_MAXSTR)  :: currentDirectory
  type(ESMF_Config)           :: config

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initializing component ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  if (.not.ESMF_StateIsCreated(importState)) then
    importState=ESMF_StateCreate(name=trim(compName)//'Import', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  if (.not.ESMF_StateIsCreated(exportState)) then
    importState=ESMF_StateCreate(name=trim(compName)//'Export', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  ! Read the configuration for this component from file if not
  ! already present in the component
  call ESMF_GridCompGet(comp, configIsPresent=isPresent, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (isPresent) then
    call ESMF_GridCompGet(comp, config=config, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message, '(A)') trim(compName)//' uses internal configuration'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  else
    configfilename=trim(compName)//'.cfg'
    inquire(file=trim(configfilename), exist=isPresent)

    config = ESMF_ConfigCreate(rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    if (isPresent) then
      call ESMF_ConfigLoadFile(config, trim(configfilename), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      write(message,'(A)')  trim(compName)//' read configuration from '//trim(configFileName)
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    else
      write(message,'(A)')  trim(compName)//' has no configuration; use global config'
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    endif

    call ESMF_GridCompSet(comp, config=config, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  ! Make a local copy of the clock if there isn't one already
  call ESMF_GridCompGet(comp, clockIsPresent=isPresent, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (isPresent) then
    call ESMF_GridCompGet(comp, clock=schismClock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  if ((isPresent .and. .not.ESMF_ClockIsCreated(schismClock)) .or. &
  (.not. isPresent)) then
    schismClock = ESMF_ClockCreate(clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_GridCompSet(comp, clock=schismClock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message, '(A)') trim(compName)//' created clock as copy of parent clock'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  endif

  ! Get VM for this component
  call ESMF_GridCompGet(comp, vm=vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_VMGet(vm, mpiCommunicator=mpi_comm, petCount=petCount, &
    localPet=localPet, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> If we run cohorts concurrently on different PET, we don't want to initialize schism
  !> for each member of a sequence but only for the first one.  To detect this,
  !> we ask for an attribute of the component
  call ESMF_AttributeGet(comp, name='cohort_index', value=cohortIndex, defaultValue=0, rc=localrc)
  call ESMF_AttributeGet(comp, name='ncohort', value=ncohort, defaultValue=1, rc=localrc)
  call ESMF_AttributeGet(comp, name='runhours', value=runhours, defaultValue=1, rc=localrc)
  call ESMF_AttributeGet(comp, name='schism_dt2', value=schism_dt2, defaultValue=1, rc=localrc)

!Debug
!    write(message,*) trim(compName)//' cohort_index=',cohortIndex
!    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

#ifndef ESMF_MPIUNI
  if (cohortIndex == 0) then
    ! Not serial mode; initialize schism's MPI
!    call MPI_Comm_dup(mpi_comm, schism_mpi_comm, rc)
!    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call parallel_init(communicator=mpi_comm)
  endif !cohortIndex
#endif

!Debug
!    write(message,*) trim(compName)//' reached here'
!    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  !> The input directory is by default '.'.  If present as an attribute,
  !> this one is used, and if also present in the config file, the latter
  !> one is used.
  call ESMF_AttributeGet(comp, name='input_directory', &
      value=simulationDirectory, defaultValue='.', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!Debug
!  write(message,'(A)')'Input dir='//trim(simulationDirectory)
!  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  call ESMF_ConfigGetAttribute(config, simulationDirectory, &
    label='simulationDirectory:', default=trim(simulationDirectory), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

! Construct the full path to the simulation directory if it is not an
! absolute path starting with slash or backslash
  call ESMF_UtilIOGetCWD(currentDirectory, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (simulationDirectory(1:1) /= '/' .and. simulationDirectory(1:1) /= '\\') then
    simulationDirectory = trim(currentDirectory)//'/'//trim(simulationDirectory)
  endif

  ! call initialize model with parameters iths=0, ntime=0
  write(message, '(A,I4,A,A)') trim(compName)//' initializing science model on ', &
    petCount, ' PET in ', trim(simulationDirectory)
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  if (cohortIndex == 0) then !First task in the sequence that uses same PETs; init
    call schism_init(0,trim(simulationDirectory), iths, ntime)
  else !not first task
    !Bypass alloc and decomp
    call schism_init(1,trim(simulationDirectory), iths, ntime)
  endif
#ifdef USE_PDAF
! CWB2021-1 start
! ADDing PDAF C-pre, same for other changes
! debug
  write(message,*) 'init_P1: ifile=',ifile
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
! Check
  if (nc_out>0) then
      write(message,*) 'Change nc_out to 0 in param.nml'
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end if
! write(*,*) 'outf=',outf
! if ((outf==2).or.(outf==3)) then ! control output
  ! outf is initialized (in init_pdaf) after this routine, so this if statement
  ! has to be removed!
     call fill_nc_header(0) !ncfile output init
! end if
! CWB2021-1 end
#endif

  !Check consistency in inputs
  if(abs(runhours-rnday*24)>1.e-5.or.abs(schism_dt2-dt)>1.e-5) then
    !write(message,*) 'init_P1: Check rnday, dt;',runhours,rnday*24,schism_dt2,dt
    write(message,*) 'init_P1: Check rnday, dt;',runhours,rnday*24,schism_dt2,dt
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_VAL_OUTOFRANGE
    !_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  !Save init cohort states to prep for stepping
  call schism_alloc_arrays(ncohort)
  call schism_save_state(cohortIndex)

  !Save output stack ID to avoid jumbled outputs
  if(.not.allocated(istack)) allocate(istack(0:ncohort-1))
  istack(cohortIndex)=ifile

  write(message, '(A)') trim(compName)//' initialized science model'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! Use schism time interval
  call ESMF_TimeIntervalSet(schism_dt, s_r8=dt, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockSet(schismClock, timeStep=schism_dt, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

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

  !allocate(tmpIdx(npa-np), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  allocate(tmpIdx2(npa), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !write(0,*) __LINE__, localPet, np, npa, ne, nea

  ! do i=1,np
  !   tmpIdx2(i)=i
  ! end do
  ! do i=np+1,npa-np
  !   ! check existence in local elements
  !   do ie=1,ne
  !     if (any(elnode(1:i34(ie),ie)==i)) then
  !       numLocalNodes = numLocalNodes+1
  !       tmpIdx2(numLocalNodes)=i
  !       write(0,*) __LINE__,i,ie
  !       write(0,*) 'We should never have arrived here'
  !       stop
  !     else
  !       write(0,*) __LINE__,i,ie, i-np, ubound(tmpIdx)
  !       numNodeHaloIdx = numNodeHaloIdx+1
  !       tmpIdx(numNodeHaloIdx)=i-np
  !     end if
  !   end do
  ! end do

  !allocate(nodeHaloIdx(numNodeHaloIdx))
  !nodeHaloIdx(:)=tmpIdx(1:numNodeHaloIdx)
  !deallocate(tmpIdx)
  allocate(localNodes(numLocalNodes))
  do i=1, numLocalNodes
    localNodes(i) = i
  end do

!  localNodes(:)=tmpIdx2(1:numLocalNodes)
  !deallocate(tmpIdx2)

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

  !allocate(elementcoords3d(3*nea), stat=localrc)
  allocate(nv(4*ne), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! set ESMF coordSys type
  if (ics==2) then
    coordsys=ESMF_COORDSYS_SPH_DEG
  else
    write(message, '(A)') trim(compName)//' uses a cartesian coordinate system'// &
       'which may not be suitable for coupling'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    coordsys=ESMF_COORDSYS_CART
  endif

  do ip=1, numLocalNodes
    i = localNodes(ip)
    ! iplg(i) is global node index of local node i in the augmented domain
    nodeids(ip)=iplg(i)
    if (ics==2) then
      ! if geographical coordinates present
      nodecoords2d(2*ip-1) = rad2deg*xlon(i)
      nodecoords2d(2*ip)   = rad2deg*ylat(i)
    else
      ! use cartesian coordinates
      nodecoords2d(2*ip-1) = xnd(i)
      nodecoords2d(2*ip)   = ynd(i)
    end if
    if (i<=np) then
      ! if iplg(i) is resident, ipgl(iplg(i))%rank=myrank
      nodeowners(ip) = ipgl(iplg(i))%rank
    else
      ! get owner of foreign node, the SCHISM manual tells us to not use
      ! ipgl for foreign nodes ???
      ! ipgl%next%next%next.... is the linked list, with ranks in ascending order.  Unless ipgb is an interface node (i.e., resident in more than 1 process), the list has only 1 entry
      nextp => ipgl(iplg(i))
      ! We here advance to the end of the list, Joseph suggested to just take the next element
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
  call ESMF_MeshGet(mesh2d,numOwnedNodes=mynp,numOwnedElements=myne,elementDistgrid=distgrid,rc=localrc)
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

  !> Create states if they were not created


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

  ! define fields for import and export
  schism_ptr2d => windx2(1:np)
  fieldName = 'wind_x-velocity_in_10m_height'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           farrayPtr=schism_ptr2d, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateAddReplace(importstate, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created import field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  schism_ptr2d => windy2(1:np)
  fieldName = 'wind_y-velocity_in_10m_height'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           farrayPtr=schism_ptr2d, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateAddReplace(importState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created import field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  schism_ptr2d => pr2(1:np)
  fieldName = 'air_pressure_at_water_surface'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           farrayPtr=schism_ptr2d, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateAddReplace(importState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created import field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  schism_ptr2d => airt2(1:np)
  fieldName = 'temperature_in_air'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           farrayPtr=schism_ptr2d, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateAddReplace(importState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created import field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  schism_ptr2d => shum2(1:np)
  fieldName = 'specific_humidity_in_air'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           farrayPtr=schism_ptr2d, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateAddReplace(importState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created import field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  schism_ptr2d => srad(1:np)
  fieldName = 'downwelling_short_wave_radiation_at_water_surface'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           farrayPtr=schism_ptr2d, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateAddReplace(importState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created import field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  schism_ptr2d => fluxevp(1:np)
  !> @todo ist this the correct name for fluxevp?
  fieldName = 'evaporation_flux_at_water_surface' !'temperature_at_water_surface'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           farrayPtr=schism_ptr2d, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateAddReplace(importState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created import field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  schism_ptr2d => fluxprc(1:np)
  fieldName = 'precipitation_flux_at_water_surface'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           farrayPtr=schism_ptr2d, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateAddReplace(importState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created import field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! fill export state

  fieldName = 'temperature_at_water_surface'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           typekind=ESMF_TYPEKIND_R8, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  ! add maskValues to be used in regridding
  call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !   initialize
  call ESMF_FieldGet(field,farrayPtr=schism_ptr2d,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  schism_ptr2d(1:np)=tr_nd(1,nvrt,1:np)
  call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  !Bottom T (for benthic model)
  name = 'temperature'
  call schism_esmf_add_bottom_tracer(name,mesh2d,1,exportState)
  write(message, '(A,A)') trim(compName)//' added as bottom tracer "', &
    trim(name)//'"'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  fieldName = 'water_x-velocity_at_water_surface'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           typekind=ESMF_TYPEKIND_R8, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  ! add maskValues to be used in regridding
  call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !   initialize
  call ESMF_FieldGet(field,farrayPtr=schism_ptr2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  schism_ptr2d(1:np)=uu2(nvrt,1:np)
  call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  fieldName = 'water_y-velocity_at_water_surface'
  field = ESMF_FieldCreate(mesh2d, name=fieldName, &
                           typekind=ESMF_TYPEKIND_R8, &
                           meshloc=ESMF_MESHLOC_NODE, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  ! add maskValues to be used in regridding
  call ESMF_AttributeSet(field, name="maskValues", valueList=maskValues, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !   initialize
  call ESMF_FieldGet(field,farrayPtr=schism_ptr2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  schism_ptr2d(1:np)=vv2(nvrt,1:np)
  call ESMF_StateAddReplace(exportState, (/field/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A,A)') trim(compName)//' created export field "', &
    trim(fieldName)//'" on nodes'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

#ifdef USE_FABM
  do i=1,fs%nvar
    name = trim(fs%model%state_variables(i)%name)
    call schism_esmf_add_bottom_tracer(name,mesh2d,fabm_istart-1+i, exportState, add_ws=.true., importState=importState)
    write(message, '(A,A)') trim(compName)//' added as bottom tracer "', &
      trim(name)//'"'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  end do
#endif

  call addCIM(comp, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

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

  write(message, '(A)') trim(compName)//' initialized.'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

end subroutine InitializeP1

#undef ESMF_METHOD
#define ESMF_METHOD "Run"
subroutine Run(comp, importState, exportState, parentClock, rc)

  use schism_glbl, only: rnday,dt, tr_nd, nvrt, npa, np, kbp,idry,uu2,vv2,eta2,&
     &in_dir,out_dir,len_in_dir,len_out_dir,iplg,ynd,xnd,windx,windy,ifile,ihfskip,nspool,&
     &nea,nsa,id_out_var,iof_hydro,idry_e,idry_s,znl,pr,ww2,prho,airt1,shum1,srad,fluxsu,fluxlu,&
     &hradu,hradd,sflux,fluxevp,fluxprc,tau_bot_node,tau,dav,dfh,dfv,q2,xl,su2,sv2,we,tr_el,it_main
  use schism_msgp, only: myrank,nproc
  use schism_io, only: writeout_nc,fill_nc_header!ncid
#ifdef USE_PDAF
  use mod_assimilation, only: outf !PDAF module
#endif
!  USE PDAF_interfaces_module

#ifdef USE_FABM
  use fabm_schism, only: fabm_istart=>istart, fs
#endif
  implicit none

  type(ESMF_GridComp)     :: comp
  type(ESMF_State)        :: importState
  type(ESMF_State)        :: exportState
  type(ESMF_Clock)        :: parentClock
  integer, intent(out)    :: rc

  type(ESMF_Clock)        :: schismClock
  type(ESMF_Time)         :: nextTime, currTime, parentCurrTime
  type(ESMF_TimeInterval) :: timeStep
  integer                 :: i,num_schism_steps
  integer, save           :: it=1
  type(ESMF_Field)        :: field
  real(ESMF_KIND_R8), pointer :: ptr2d(:)
  character(len=ESMF_MAXSTR)  :: message,name,compName,simulationDirectory,currentDirectory
  integer(ESMF_KIND_I4)   :: localrc,cohortIndex,analysis_step
  integer(ESMF_KIND_I8)   :: advanceCount
  real(ESMF_KIND_R8)      :: dt_coupling
  type(ESMF_Config)           :: config

  INTEGER :: doexit, steps,status_pdaf
  REAL    :: timenow
  character(len=72) :: fdb
  integer :: lfdb

#ifdef USE_PDAF
! External subroutines
! comment out --> into generic interface (assimilate_pdaf)
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time,
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_ens !, &            ! User supplied pre/poststep routine
!      collect_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf, &
!      init_obs_pdaf,prodRinvA_pdaf,init_obsvar_pdaf
#endif

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, name=compName, clock=schismClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockGet(schismClock, advanceCount=advanceCount, &
    currTime=currTime, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockGet(parentClock, currTime=parentCurrTime, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (currTime /= parentCurrTime) then
    write(message, '(A)') trim(compName)//' times not synchronized: '
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    call ESMF_TimePrint(currtime, unit=message, rc=localrc)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    call ESMF_TimePrint(currtime, unit=message, rc=localrc)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif

  write(message, '(A,I6.6)') trim(compName)//' running step ',advanceCount
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! get data from fields in import state
  !   here: nothing to be done for the wind velocities,
  !         since fields are created with ESMF_DATACOPY_REFERENCE

  ! run model to nextTime of parentClock by manipulating its
  ! stopTime and running its internal timestep until this stopTime.
  call ESMF_GridCompGet(comp, clock=schismClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockGetNextTime(parentClock, nextTime,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockSet(schismClock, stopTime=nextTime, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_AttributeGet(comp, name='cohort_index', &
    value=cohortIndex, defaultValue=0, rc=localrc)

  call ESMF_AttributeGet(comp, name='analysis_step', &
    value=analysis_step, defaultValue=0, rc=localrc)

  !Reinit input dir in case same PETs are doing a different run
  call ESMF_AttributeGet(comp, name='input_directory', &
      value=simulationDirectory, defaultValue='.', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

!  call ESMF_GridCompGet(comp, config=config, rc=localrc)
!    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!
!  call ESMF_ConfigGetAttribute(config, simulationDirectory, &
!    label='simulationDirectory:', default=trim(simulationDirectory), rc=localrc)
!  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

! Construct the full path to the simulation directory if it is not an
! absolute path starting with slash or backslash
  call ESMF_UtilIOGetCWD(currentDirectory, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (simulationDirectory(1:1) /= '/' .and. simulationDirectory(1:1) /= '\\') then
    simulationDirectory = trim(currentDirectory)//'/'//trim(simulationDirectory)
  endif

  !Update in/outdir for correct dir under same PETs
  in_dir=adjustl(trim(adjustl(simulationDirectory))//'/')
  len_in_dir=len_trim(in_dir)
  out_dir=adjustl(in_dir(1:len_in_dir)//'outputs/')
  len_out_dir=len_trim(out_dir)

  num_schism_steps=rnday*86400.d0/dt+0.5d0
  it=advanceCount+1 !SCHISM time step index

!new28: not sure needed
!!  call PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
!!       distribute_state_pdaf, prepoststep_ens, status_pdaf)

  call schism_get_state(cohortIndex)

#ifdef USE_PDAF
! Put PDAF_get_state here
  call PDAF_get_state(steps,timenow, doexit, next_observation_pdaf, distribute_state_pdaf, prepoststep_ens, status_pdaf)
#endif

  !Rewind clock for forcing
  call other_hot_init(dble(it-1)*dt)

  !Reset stack #
  ifile=istack(cohortIndex)

  do while (.not. ESMF_ClockIsStopTime(schismClock, rc=localrc))
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message,*) 'it = ',it,', elapsed  ',it*dt,'s'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

    call schism_step(it)

!Debug
!    if(it==advanceCount+1) then
!      write(message,*) trim(compName)//' wind before:',cohortIndex,it,minval(windx),maxval(windx)
!      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
!    endif

    call ESMF_ClockAdvance(schismClock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    ! Output
    if(it==num_schism_steps) then
      write(message,*) 'starting output...'
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

      fdb='surface_0000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
      open(10,file=out_dir(1:len_out_dir)//trim(adjustl(fdb)),status='replace')
      write(10,*)np,nproc
      do i=1,np
        write(10,'(i11,100(1x,e20.12))')iplg(i),xnd(i),ynd(i),tr_nd(1:2,nvrt,i),uu2(nvrt,i),vv2(nvrt,i),eta2(i)
      enddo !i
      close(10)

      fdb='bottom_0000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
      open(10,file=out_dir(1:len_out_dir)//trim(adjustl(fdb)),status='replace')
      write(10,*)np,nproc
      do i=1,np
        write(10,'(i11,100(1x,e20.12))')iplg(i),xnd(i),ynd(i),tr_nd(1:2,1,i),uu2(2,i),vv2(2,i),eta2(i)
      enddo !i
      close(10)
    endif !it==

#ifdef USE_PDAF
!CWB2021-1 start
!   Writeout nc (Hydro only)
    if ((outf==2).or.(outf==3)) then ! control output
     if(mod(it,nspool)==0) then
        write(message,*) 'writeout nc at it = ',it,', elapsed ',it*dt,'s' !,' it_main=',it_main
        call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
        call writeout_nc(id_out_var(1),'wetdry_node',1,1,npa,dble(idry))
        call writeout_nc(id_out_var(2),'wetdry_elem',4,1,nea,dble(idry_e))
        call writeout_nc(id_out_var(3),'wetdry_side',7,1,nsa,dble(idry_s))
       !zcor MUST be 1st 3D var output for combine scripts to work!
        call writeout_nc(id_out_var(4),'zcor',2,nvrt,npa,znl(:,:))
        if(iof_hydro(1)==1) call writeout_nc(id_out_var(5),'elev',1,1,npa,eta2)
        if(iof_hydro(2)==1) call writeout_nc(id_out_var(6),'air_pressure',1,1,npa,pr)
        if(iof_hydro(3)==1) call writeout_nc(id_out_var(7),'air_temperature',1,1,npa,airt1)
        if(iof_hydro(4)==1) call writeout_nc(id_out_var(8),'specific_humidity',1,1,npa,shum1)
        if(iof_hydro(5)==1) call writeout_nc(id_out_var(9),'solar_radiation',1,1,npa,srad)
        if(iof_hydro(6)==1) call writeout_nc(id_out_var(10),'sensible_flux',1,1,npa,fluxsu)
        if(iof_hydro(7)==1) call writeout_nc(id_out_var(11),'latent_heat',1,1,npa,fluxlu)
        if(iof_hydro(8)==1) call writeout_nc(id_out_var(12),'upward_longwave',1,1,npa,hradu)
        if(iof_hydro(9)==1) call writeout_nc(id_out_var(13),'downward_longwave',1,1,npa,hradd)
        if(iof_hydro(10)==1) call writeout_nc(id_out_var(14),'total_heat_flux',1,1,npa,sflux)
        if(iof_hydro(11)==1) call writeout_nc(id_out_var(15),'evaporation',1,1,npa,fluxevp)
        if(iof_hydro(12)==1) call writeout_nc(id_out_var(16),'precipitation',1,1,npa,fluxprc)
        if(iof_hydro(13)==1) call writeout_nc(id_out_var(17),'bottom_stress',1,1,npa,tau_bot_node(1,:),tau_bot_node(2,:)) !Cdp)
        if(iof_hydro(14)==1) call writeout_nc(id_out_var(18),'wind_speed',1,1,npa,windx,windy)
        if(iof_hydro(15)==1) call writeout_nc(id_out_var(19),'wind_stress',1,1,npa,tau(1,:),tau(2,:))
        if(iof_hydro(16)==1) call writeout_nc(id_out_var(20),'dahv',1,1,npa,dav(1,:),dav(2,:))
        if(iof_hydro(17)==1) call writeout_nc(id_out_var(21),'vertical_velocity',2,nvrt,npa,ww2)
        if(iof_hydro(18)==1) call writeout_nc(id_out_var(22),'temp',2,nvrt,npa,tr_nd(1,:,:))
        if(iof_hydro(19)==1) call writeout_nc(id_out_var(23),'salt',2,nvrt,npa,tr_nd(2,:,:))
        if(iof_hydro(20)==1) call writeout_nc(id_out_var(24),'water_density',2,nvrt,npa,prho)
        if(iof_hydro(21)==1) call writeout_nc(id_out_var(25),'diffusivity',2,nvrt,npa,dfh)
        if(iof_hydro(22)==1) call writeout_nc(id_out_var(26),'viscosity',2,nvrt,npa,dfv)
        if(iof_hydro(23)==1) call writeout_nc(id_out_var(27),'TKE',2,nvrt,npa,q2)
        if(iof_hydro(24)==1) call writeout_nc(id_out_var(28),'mixing_length',2,nvrt,npa,xl)
        if(iof_hydro(25)==1) call writeout_nc(id_out_var(29),'hvel',2,nvrt,npa,uu2,vv2)
        if(iof_hydro(26)==1) call writeout_nc(id_out_var(30),'hvel_side',8,nvrt,nsa,su2,sv2)
        if(iof_hydro(27)==1) call writeout_nc(id_out_var(31),'wvel_elem',5,nvrt,nea,we)
        if(iof_hydro(28)==1) call writeout_nc(id_out_var(32),'temp_elem',6,nvrt,nea,tr_el(1,:,:))
        if(iof_hydro(29)==1) call writeout_nc(id_out_var(33),'salt_elem',6,nvrt,nea,tr_el(2,:,:))
!       if(iof_hydro(30)==1) call writeout_nc(id_out_var(34),'pressure_gradient',7,1,nsa,bpgr(:,1),bpgr(:,2))
!       bpgr is not in schism_glbl, skip it 
     end if !mod

!   Close nc files
     if(mod(it,ihfskip)==0) then
       ifile=ifile+1  !output file #
       call fill_nc_header(1)
     endif !it==ifile*ihfskip
!    debug
!    write(message,*) 'Run: ifile=',ifile
!    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    end if ! outf
    
!CWB2021-1 end
#endif

    it=it+1
  end do !while

!Debug
!  write(message,*) trim(compName)//' wind after:',minval(windx),maxval(windx),minval(windy),maxval(windy)
!  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  !Check if it's time for analysis
!new28
#ifdef USE_PDAF
! if(analysis_step/=0) then
     write(message,*)trim(compName)//' entering PDAF assimilate, ',it
     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
!    Using assimilate PDAF interface
     call assimilate_pdaf() !
!    call assimilate_pdaf(pdaf_stat)
!    localrc=pdaf_stat !pdaf_stat is using i2, localrc is i4
     write(message,*)'Done PDAF assimilate'
     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

! endif !analysis_step/=0
#endif

! Move save_state after DA
  call schism_save_state(cohortIndex)

  !Save stack #
  istack(cohortIndex)=ifile


  !> Do a clock correction for non-integer internal timesteps and issue
  !> a warning.  We could improve this somewhat by choosing to optionally
  !> run one less timestep if that gives better accuracy.
  call ESMF_ClockGet(schismClock, currTime=currTime, rc=localrc)

  if (currTime /= nextTime) then

    timeStep = nextTime-currTime

    call ESMF_TimeIntervalGet(timeStep, s_r8=dt_coupling)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message, '(A,I4,A)') trim(compName)//' corrected clock by ', &
      int(dt_coupling), ' seconds'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
    call ESMF_ClockSet(schismClock, currTime=nextTime, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  timeStep = nextTime-parentCurrTime

  call ESMF_TimeIntervalGet(timeStep, s_r8=dt_coupling)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)


  ! update 3d grid
  ! (so far not necessary for 2d variables to be exchanged)

  ! put data to fields in export state
  !   SST, copy uncontiguous data into field in ESMF
  call ESMF_StateGet(exportState, 'temperature_at_water_surface', field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, farrayPtr=ptr2d,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ptr2d(1:np) = tr_nd(1,nvrt,1:np)

  name = 'temperature'
  call schism_esmf_update_bottom_tracer(exportState, name, 1)

#ifdef USE_FABM
  do i=1,fs%nvar
    name = trim(fs%model%state_variables(i)%name)
    call schism_esmf_update_bottom_tracer(exportState, name, fabm_istart-1+i, importState=importState, dt=dt_coupling)
  end do
#endif

  call ESMF_StateGet(exportState, 'water_x-velocity_at_water_surface', field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call ESMF_FieldGet(field, farrayPtr=ptr2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  ptr2d(1:np) = uu2(nvrt,1:np)

  call ESMF_StateGet(exportState, 'water_y-velocity_at_water_surface', field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call ESMF_FieldGet(field, farrayPtr=ptr2d, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  ptr2d(1:np) = vv2(nvrt,1:np)

!  write(0,*) '  Ended Run from SCHISM'
  write(message, '(A)') trim(compName)//' ran.'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

end subroutine Run

#undef ESMF_METHOD
#define ESMF_METHOD "Finalize"
  subroutine Finalize(comp, importState, exportState, clock, rc)

  implicit none

  type(ESMF_GridComp)   :: comp
  type(ESMF_State)      :: importState
  type(ESMF_State)      :: exportState
  type(ESMF_Clock)      :: clock
  integer, intent(out)  :: rc

  type(ESMF_Clock)      :: schismClock
  type(ESMF_Mesh)       :: mesh
  type(ESMF_Field)      :: field
  type(ESMF_Distgrid)   :: distgrid
  integer(ESMF_KIND_I4) :: localrc
  character(ESMF_MAXSTR):: compName, message

  rc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  write(message, '(A)') trim(compName)//' finalizing ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  call ESMF_StateGet(importState, 'wind_x-velocity_in_10m_height', field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldGet(field, mesh=mesh, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_FieldDestroy(field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_MeshGet(mesh,elementDistgrid=distgrid,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_MeshDestroy(mesh,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call ESMF_DistgridDestroy(distgrid,rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGet(comp, clock=schismClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockDestroy(schismClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' finalized.'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  end subroutine Finalize

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

subroutine schism_alloc_arrays(ncohort)
  use schism_glbl, only: nea,nsa,npa,nvrt,ntracers,nea2
  implicit none

  integer(ESMF_KIND_I4), intent(in) :: ncohort

  !ncohort>=1 checked
  if(.not.allocated(idry_e_c)) allocate(idry_e_c(nea,0:ncohort-1))
  if(.not.allocated(idry_s_c)) allocate(idry_s_c(nsa,0:ncohort-1))
  if(.not.allocated(idry_c)) allocate(idry_c(npa,0:ncohort-1))

  if(.not.allocated(we_c)) allocate(we_c(nvrt,nea,0:ncohort-1))
  if(.not.allocated(tr_el_c)) allocate(tr_el_c(ntracers,nvrt,nea2,0:ncohort-1))
  if(.not.allocated(su2_c)) allocate(su2_c(nvrt,nsa,0:ncohort-1))
  if(.not.allocated(sv2_c)) allocate(sv2_c(nvrt,nsa,0:ncohort-1))
  if(.not.allocated(eta2_c)) allocate(eta2_c(npa,0:ncohort-1))
  if(.not.allocated(tr_nd_c)) allocate(tr_nd_c(ntracers,nvrt,npa,0:ncohort-1))
  if(.not.allocated(tr_nd0_c)) allocate(tr_nd0_c(ntracers,nvrt,npa,0:ncohort-1))
  if(.not.allocated(q2_c)) allocate(q2_c(nvrt,npa,0:ncohort-1))
  if(.not.allocated(xl_c)) allocate(xl_c(nvrt,npa,0:ncohort-1))
  if(.not.allocated(dfv_c)) allocate(dfv_c(nvrt,npa,0:ncohort-1))
  if(.not.allocated(dfh_c)) allocate(dfh_c(nvrt,npa,0:ncohort-1))
  if(.not.allocated(dfq1_c)) allocate(dfq1_c(nvrt,npa,0:ncohort-1))
  if(.not.allocated(dfq2_c)) allocate(dfq2_c(nvrt,npa,0:ncohort-1))

  !Extra: somehow this helps
  if(.not.allocated(uu2_c)) allocate(uu2_c(nvrt,npa,0:ncohort-1))
  if(.not.allocated(vv2_c)) allocate(vv2_c(nvrt,npa,0:ncohort-1))
  !CWB2021-1 Add ncid
  if(.not.allocated(ncid_c)) allocate(ncid_c(0:ncohort-1))
end subroutine schism_alloc_arrays

subroutine schism_save_state(ci)
  use schism_glbl, only: nea,nsa,npa,nvrt,ntracers,idry_e,we,tr_el, &
 &idry_s,su2,sv2, idry,eta2,tr_nd,tr_nd0,q2,xl,dfv,dfh,dfq1,dfq2, &
 &uu2,vv2
  use schism_io, only:ncid_schism_io !CWB2021-1
  implicit none
  integer(ESMF_KIND_I4), intent(in) :: ci !cohort index (0-based)

  character(len=ESMF_MAXSTR)  :: message
  integer(ESMF_KIND_I4)              :: localrc,rc

!Debug
!  if(ci<0.or.ci>1) then
!    write(message,*) 'cohort index exceeded:',ci
!    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
!    localrc = ESMF_RC_VAL_OUTOFRANGE
!    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!  endif

  idry_e_c(:,ci)=idry_e
  idry_s_c(:,ci)=idry_s
  idry_c(:,ci)=idry

  we_c(:,:,ci)=we
  tr_el_c(:,:,:,ci)=tr_el
  su2_c(:,:,ci)=su2
  sv2_c(:,:,ci)=sv2
  eta2_c(:,ci)=eta2
  tr_nd_c(:,:,:,ci)=tr_nd
  tr_nd0_c(:,:,:,ci)=tr_nd0
  q2_c(:,:,ci)=q2
  xl_c(:,:,ci)=xl
  dfv_c(:,:,ci)=dfv
  dfh_c(:,:,ci)=dfh
  dfq1_c(:,:,ci)=dfq1
  dfq2_c(:,:,ci)=dfq2

  uu2_c(:,:,ci)=uu2
  vv2_c(:,:,ci)=vv2
  ncid_c(ci)=ncid_schism_io !CWB2021-1
end subroutine schism_save_state

subroutine schism_get_state(ci)
  use schism_glbl, only: nea,nsa,npa,nvrt,ntracers,idry_e,we,tr_el, &
 &idry_s,su2,sv2, idry,eta2,tr_nd,tr_nd0,q2,xl,dfv,dfh,dfq1,dfq2, &
 &uu2,vv2
  use schism_io, only:ncid_schism_io !CWB2021-1
  implicit none
  integer(ESMF_KIND_I4), intent(in) :: ci !cohort index (0-based)

  character(len=ESMF_MAXSTR)  :: message
  integer(ESMF_KIND_I4)              :: localrc,rc

!Debug
!  if(ci<0.or.ci>1) then
!    write(message,*) 'cohort index exceeded:',ci
!    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
!    localrc = ESMF_RC_VAL_OUTOFRANGE
!    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
!  endif

  idry_e=idry_e_c(:,ci)
  idry_s=idry_s_c(:,ci)
  idry=idry_c(:,ci)

  we=we_c(:,:,ci)
  tr_el=tr_el_c(:,:,:,ci)
  su2=su2_c(:,:,ci)
  sv2=sv2_c(:,:,ci)
  eta2=eta2_c(:,ci)
  tr_nd=tr_nd_c(:,:,:,ci)
  tr_nd0=tr_nd0_c(:,:,:,ci)
  q2=q2_c(:,:,ci)
  xl=xl_c(:,:,ci)
  dfv=dfv_c(:,:,ci)
  dfh=dfh_c(:,:,ci)
  dfq1=dfq1_c(:,:,ci)
  dfq2=dfq2_c(:,:,ci)

  uu2=uu2_c(:,:,ci)
  vv2=vv2_c(:,:,ci)
  ncid_schism_io=ncid_c(ci)
end subroutine schism_get_state

end module schism_cmi_esmf
