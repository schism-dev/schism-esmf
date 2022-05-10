! This code is part of the SCHISM-ESMF interface. It is a main
! program for running three schism components concurrently to
! a dummy atmosphere
!
! @copyright (C) 2018, 2019, 2020-2021 Helmholtz-Zentrum Geesthacht
! @author Richard Hofmeister richard.hofmeister@hereon.de
! @author Carsten Lemmen carsten.lemmen@hereon.de
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
#define ESMF_FILENAME "triple_schism.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#define ESMF_METHOD "main"
program main

  use esmf
  use schism_esmf_util, only: clockCreateFrmParam
  use schism_cmi_esmf, only: schismSetServices => SetServices
  use atmosphere_cmi_esmf,  only: atmosSetServices => SetServices

  implicit none

  type(ESMF_GridComp), allocatable :: schism_components(:)
  type(ESMF_GridComp)     :: atmos_component

  type(ESMF_State), allocatable   :: schism_imports(:), schism_exports(:)
  type(ESMF_State)        :: atmos_import, atmos_export

  type(ESMF_TimeInterval) :: timestep
  type(ESMF_Time)         :: start_time, stop_time
  type(ESMF_Clock)        :: clock

  type(ESMF_Field)               :: field, field_in
  type(ESMF_Field), allocatable  :: fields_out(:)

  type(ESMF_RouteHandle)  :: routehandle_sea2air
  type(ESMF_RouteHandle), allocatable  :: routehandles_air2sea(:)
  type(ESMF_Vm)           :: vm

  integer(ESMF_KIND_I4)       :: petCountLocal, schismCount=3
  integer(ESMF_KIND_I4)       :: rc, petCount, i, j, inum, localrc
  integer(ESMF_KIND_I4), allocatable    :: petlist(:)
  real(ESMF_KIND_R8), pointer :: ptr1d(:)
  logical                     :: isPresent
  character(len=ESMF_MAXSTR)  :: filename, message

  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Inquire the parallel environment about available
  ! resources, and partition the environment to use
  ! all but one PET for SCHISM, and last core for the
  ! dummy atmosphere component

  call ESMF_VMGetGlobal(vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_VMGet(vm, petCount=petCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Create all components on their respective parallel
  ! environment provided by each petList
  allocate(schism_components(schismCount))
  do i = 1, schismCount

    petCountLocal = petCount/max(1, schismCount)
    allocate(petlist(petCountLocal), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    do j=1, petCountLocal
      petList(j)=(i-1)*petCountLocal + j - 1 ! 0 based
    end do


    write(message, '(A,I1)') 'schism_', i
    !write(0, *) trim(message), 'list=', petList, 'petCount=', petCount, petCountLocal

    schism_components(i) = ESMF_GridCompCreate(name=trim(message), &
      petList=petlist, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    deallocate(petList)
  end do

  allocate(schism_exports(schismCount))
  allocate(schism_imports(schismCount))

  do i = 1, schismCount

    call ESMF_GridCompSetServices(schism_components(i), &
      schismSetServices, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    schism_exports(i) = ESMF_StateCreate(name='schism export state', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    schism_imports(i) = ESMF_StateCreate(name='schism import state', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do

  atmos_component = ESMF_GridCompCreate(name='atmosphere_1', &
    petList=(/petCount - 1/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  atmos_export = ESMF_StateCreate(name='atmos export state', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  atmos_import = ESMF_StateCreate(name='atmos import state', rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetServices(atmos_component, &
    atmosSetServices, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)


  filename = './global.nml'
  clock = clockCreateFrmParam(filename, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i = 1, schismCount

    call ESMF_GridCompInitialize(schism_components(i), importState=schism_imports(i), &
      exportState=schism_exports(i), clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do

  call ESMF_GridCompInitialize(atmos_component, importState=atmos_import, &
    exportState=atmos_export, clock=clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !Make sure that all states are reconciled across the
  ! entire VM
  do i = 1, schismCount

    call ESMF_StateReconcile(schism_imports(i), vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateReconcile(schism_exports(i), vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  end do

  call ESMF_StateReconcile(atmos_export, vm=vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Within the schism component, the following fields are
  ! defined for import and export (y-component not implemented yet)
  call ESMF_StateGet(atmos_export,'wind_x-velocity', &
    field=field_in, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(fields_out(schismCount))
  allocate(routehandles_air2sea(schismCount))

  do i = 1, schismCount
    call ESMF_StateGet(schism_imports(i), 'wind_x-velocity_in_10m_height', &
      field=fields_out(i), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    ! Precompute the weights
    call ESMF_FieldRegridStore(field_in, fields_out(i), &
      routehandle=routehandles_air2sea(i), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  enddo

  ! Loop over coupling timesteps until stopTime
  do while ( .not. (ESMF_ClockIsStopTime(clock)))

    !> Directly manipulate the fields from import and export
    !> states.  In less basic applications, this should be
    !> handled by a mediator component.
    do i=1, schismCount
      call ESMF_FieldRegrid(field_in, fields_out(i), &
          routeHandle=routehandles_air2sea(i), rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    enddo

    call ESMF_GridCompRun(atmos_component, importState=atmos_import, &
      exportState=atmos_export, clock=clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    do i=1, schismCount
      call ESMF_GridCompRun(schism_components(i), importState=schism_imports(i), &
        exportState=schism_exports(i), clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    enddo

    call ESMF_ClockAdvance(clock, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end do

  ! Clean up

  deallocate(fields_out)
  deallocate(routehandles_air2sea)

  do i = 1, schismCount

    call ESMF_GridCompFinalize(schism_components(i), importState=schism_imports(i), &
      exportState=schism_exports(i), clock=clock, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateDestroy(schism_imports(i), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateDestroy(schism_exports(i), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_GridCompDestroy(schism_components(i), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end do
  deallocate(schism_exports)
  deallocate(schism_imports)
  deallocate(schism_components)

  call ESMF_ClockDestroy(clock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateDestroy(atmos_import, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateDestroy(atmos_export, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompDestroy(atmos_component, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_Finalize(rc=localrc)

end program main
