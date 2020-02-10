! dummy grid component
!
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister

module dummy_grid_component
use esmf
use mpi

implicit none
integer               :: iths=0,ntime=0

public SetServices

contains

  subroutine SetServices(comp,rc)
  type(ESMF_GridComp) :: comp
  integer,intent(out) :: rc

  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
                            userRoutine=Init, rc=rc)
  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
                            userRoutine=Run, rc=rc)
  call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
                            userRoutine=Finalize, rc=rc)

  rc = ESMF_SUCCESS
  end subroutine SetServices


  subroutine Init(comp, importState, exportState, clock, rc)
  type(ESMF_GridComp)   :: comp
  type(ESMF_State)      :: importState
  type(ESMF_State)      :: exportState
  type(ESMF_Field)      :: field
  type(ESMF_Clock)      :: clock
  type(ESMF_VM)         :: vm
  type(ESMF_DistGrid)   :: distgrid
  type(ESMF_Grid)       :: grid
  type(ESMF_Array)      :: xarray,yarray,dataarray
  real(ESMF_KIND_R8), dimension(:,:),pointer :: ptr2d
  integer, intent(out)  :: rc
  integer               :: i,n,petcount

  ! get communicator
  call ESMF_GridCompGet(comp,vm=vm,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_VMGet(vm,petCount=petcount,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !Simple box grid
  grid  = ESMF_GridCreateNoPeriDim( &
         minIndex=(/1,1/),maxIndex=(/2,2/),  coordsys=ESMF_COORDSYS_CART, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_GridGet(grid,distgrid=distgrid,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  xarray = ESMF_ArrayCreate(distgrid,ESMF_TYPEKIND_R8,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_ArrayGet(xarray,farrayPtr=ptr2d,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  ptr2d(1,1)=0.0
  ptr2d(1,2)=0.0
  ptr2d(2,1)=20000.0
  ptr2d(2,2)=20000.0

  yarray = ESMF_ArrayCreate(distgrid,ESMF_TYPEKIND_R8,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_ArrayGet(yarray,farrayPtr=ptr2d,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  ptr2d(1,1)=0.0
  ptr2d(1,2)=2000.0
  ptr2d(2,1)=0.0
  ptr2d(2,2)=2000.0

  call ESMF_GridSetCoord(grid, &
    staggerLoc=ESMF_STAGGERLOC_CENTER, coordDim=1, array=xarray, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_GridSetCoord(grid, &
    staggerLoc=ESMF_STAGGERLOC_CENTER, coordDim=2,array=yarray,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! define fields for export
  field = ESMF_FieldCreate(grid, name='wind_x-velocity', &
                           typekind=ESMF_TYPEKIND_R8, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_FieldGet(field,farrayPtr=ptr2d,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  ptr2d(1,1)=0.0
  ptr2d(1,2)=3.0
  ptr2d(2,1)=-3.0
  ptr2d(2,2)=1.0

  call ESMF_StateAddReplace(exportstate, (/field/), rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  rc = ESMF_SUCCESS
  end subroutine Init



  subroutine Run(comp, importState, exportState, clock, rc)
  implicit none
  type(ESMF_GridComp)     :: comp
  type(ESMF_State)        :: importState
  type(ESMF_State)        :: exportState
  type(ESMF_Clock)        :: clock
  integer, intent(out)    :: rc


  rc = ESMF_SUCCESS
  end subroutine Run



  subroutine Finalize(comp, importState, exportState, clock, rc)
  implicit none
  type(ESMF_GridComp)   :: comp
  type(ESMF_State)      :: importState
  type(ESMF_State)      :: exportState
  type(ESMF_Clock)      :: clock
  type(ESMF_Grid)       :: grid
  type(ESMF_Field)      :: field
  integer, intent(out)  :: rc

  call ESMF_StateGet(importState,'wind_x-velocity',field=field,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_FieldGet(field,grid=grid,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_FieldDestroy(field,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
  call ESMF_GridDestroy(grid,rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  rc = ESMF_SUCCESS
  end subroutine Finalize


end module dummy_grid_component
