! a main program for concurrently running schism component in esmf
!
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister

program concurrent_esmf_test

use esmf
use schism_esmf_component, only: schismSetServices=>SetServices
use dummy_grid_component, only: atmosSetServices=>SetServices
implicit none

integer :: rc
type(ESMF_GridComp)     :: schism_component
type(ESMF_GridComp)     :: atmos_component
type(ESMF_State)        :: schism_import, schism_export
type(ESMF_State)        :: atmos_import, atmos_export
type(ESMF_TimeInterval) :: timestep
type(ESMF_Time)         :: start_time, stop_time
type(ESMF_Clock)        :: clock
type(ESMF_FIELD)        :: field,field_in,field_out
type(ESMF_RouteHandle)  :: routehandle_air2sea, routehandle_sea2air
type(ESMF_VM)           :: vm
integer                 :: petcount
integer,allocatable     :: petlist_schism(:),petlist_atmos(:)
integer                 :: i, inum
real(ESMF_KIND_R8), pointer :: ptr1d(:)

call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN)

call ESMF_VMGetGlobal(vm,rc=rc)
call ESMF_VMGet(vm,petCount=petcount,rc=rc)
allocate(petlist_schism(max(1,petcount-1)))
do i=1,max(1,petcount-1)
  petlist_schism(i)=i-1
end do
allocate(petlist_atmos(1))
petlist_atmos(1) = petcount-1


schism_component = ESMF_GridCompCreate(name='schism component',petList=petlist_schism,rc=rc)
if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)

atmos_component = ESMF_GridCompCreate(name='atmosphere component',petList=petlist_atmos,rc=rc)
if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)

call ESMF_GridCompSetServices(schism_component, schismSetServices, rc=rc)
if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)
call ESMF_GridCompSetServices(atmos_component, atmosSetServices, rc=rc)
if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)

schism_export = ESMF_StateCreate(name='schism export state',rc=rc)
schism_import = ESMF_StateCreate(name='schism import state',rc=rc)
atmos_export = ESMF_StateCreate(name='atmos export state',rc=rc)
atmos_import = ESMF_StateCreate(name='atmos import state',rc=rc)

! Create a clock

call ESMF_TimeIntervalSet(timestep, s=3600, rc=rc)
call ESMF_TimeSet(start_time, yy=2018, mm=4, dd=10, rc=rc)
stop_time = start_time + 2*timestep

clock = ESMF_ClockCreate(timestep, start_time, stopTime=stop_time, &
          name='main clock',rc=rc)
if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)

! Initialize SCHISM
call ESMF_GridCompInitialize(schism_component, importState=schism_import, &
         exportState=schism_export, clock=clock, rc=rc)
call ESMF_GridCompInitialize(atmos_component, importState=atmos_import, &
         exportState=atmos_export, clock=clock, rc=rc)

call ESMF_StateReconcile(schism_import,vm=vm,rc=rc)
if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)
call ESMF_StateReconcile(atmos_export,vm=vm,rc=rc)
if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)

call ESMF_StateGet(schism_import,'wind_x-velocity_in_10m_height',field=field_out,rc=rc)
if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)
call ESMF_StateGet(atmos_export,'wind_x-velocity',field=field_in,rc=rc)
if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)
call ESMF_FieldRegridStore(field_in,field_out,routehandle=routehandle_air2sea,rc=rc)
if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)

! Run SCHISM until StopTime
do while ( .not. (ESMF_ClockIsStopTime(clock)))

  call ESMF_FieldRegrid(field_in,field_out,routeHandle=routehandle_air2sea,rc=rc)
  if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)

#if 0
  ! schism should run in single-processor mode (so max.2 processes for this
  ! test); hgrid.gr3 has 63 nodes, schism runs in cartesian coordinates,
  ! coordinates of test nodes:
  !  [1]  = (0.    ,0.   )
  !  [3]  = (0.    ,1000.)
  !  [61] = (20000.,0.   )
  !  [63] = (20000.,1000.)
  call ESMF_FieldGet(field_out,farrayPtr=ptr1d, rc=rc)
  if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)
  write(0,*) 'interpolated wind field:'
  write(0,*) '  expected [1]=0.0, got: ',ptr1d(1)
  write(0,*) '  expected [3]=1.5, got: ',ptr1d(3)
  write(0,*) '  expected [61]=-3.0, got: ',ptr1d(61)
  write(0,*) '  expected [63]=-1.0, got: ',ptr1d(63)
#endif

  call ESMF_GridCompRun(atmos_component, importState=atmos_import, &
           exportState=atmos_export, clock=clock, rc=rc)
  if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)

  call ESMF_GridCompRun(schism_component, importState=schism_import, &
           exportState=schism_export, clock=clock, rc=rc)
  if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc)

  call ESMF_ClockAdvance(clock)
end do

! Finalize SCHISM
call ESMF_GridCompFinalize(schism_component, importState=schism_import, &
         exportState=schism_export, clock=clock, rc=rc)

call ESMF_ClockDestroy(clock)
call ESMF_StateDestroy(schism_import)
call ESMF_StateDestroy(schism_export)
call ESMF_StateDestroy(atmos_import)
call ESMF_StateDestroy(atmos_export)
call ESMF_GridCompDestroy(schism_component)
call ESMF_GridCompDestroy(atmos_component)

call ESMF_Finalize()

end program

