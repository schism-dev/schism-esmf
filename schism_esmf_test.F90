! a main program for coupled applications with esmf
!
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister

#define ESMF_CONTEXT  line=__LINE__,file=ESMF_FILENAME,method=ESMF_METHOD
#define ESMF_ERR_PASSTHRU msg="SCHISM subroutine call returned error"
#undef ESMF_FILENAME
#define ESMF_FILENAME "schism_esmf_test.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

#undef  ESMF_METHOD
#define ESMF_METHOD "Main"
program schism_esmf_test

use esmf
use schism_esmf_component, only: schismSetServices=>SetServices
implicit none

integer :: rc
type(ESMF_GridComp)     :: schism_component
type(ESMF_State)        :: schism_import, schism_export
type(ESMF_TimeInterval) :: timestep
type(ESMF_Time)         :: start_time, stop_time
type(ESMF_Clock)        :: clock
real(ESMF_KIND_R8), pointer :: windx(:)
type(ESMF_FIELD)        :: field
integer                 :: i, inum
integer(ESMF_KIND_I4)   :: localrc

call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

schism_component = ESMF_GridCompCreate(name='schism component', rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

call ESMF_GridCompSetServices(schism_component, schismSetServices, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

schism_export = ESMF_StateCreate(name='schism export state', rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

schism_import = ESMF_StateCreate(name='schism import state', rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

! Create a clock

call ESMF_TimeIntervalSet(timestep, s=3600, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

call ESMF_TimeSet(start_time, yy=2018, mm=4, dd=10, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

stop_time = start_time + timestep

clock = ESMF_ClockCreate(timestep, start_time, stopTime=stop_time, &
          name='main clock',rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

! Initialize SCHISM
call ESMF_GridCompInitialize(schism_component, importState=schism_import, &
         exportState=schism_export, clock=clock, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

#if 0
! set pointer for wind field
call ESMF_StateGet(schism_import, &
  'wind_x-velocity_in_10m_height',field=field,rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

call ESMF_FieldGet(field,farrayPtr=windx,rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

inum=size(windx)
do i=1,inum
  windx(i) = 10.0*i/inum
end do
#endif

! Run SCHISM until StopTime
call ESMF_GridCompRun(schism_component, importState=schism_import, &
         exportState=schism_export, clock=clock, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

! Finalize SCHISM
call ESMF_GridCompFinalize(schism_component, importState=schism_import, &
         exportState=schism_export, clock=clock, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

call ESMF_ClockDestroy(clock, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
call ESMF_StateDestroy(schism_import, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
call ESMF_StateDestroy(schism_export, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
call ESMF_GridCompDestroy(schism_component, rc=localrc)
_SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

call ESMF_Finalize()

end program
