! This code is part of the SCHISM-ESMF interface.  It defines
! the schism component for a NUOPC coupled system
!
! @copyright (C) 2021-2023 Helmholtz-Zentrum Hereon
! @copyright (C) 2022-2023 Virginia Institute of Marine Science
! @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
!
! @author Carsten Lemmen <carsten.lemmen@hereon.de>
! @author Joseph Y. Zhang >jzhang@vims.edu>
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
#define ESMF_FILENAME "schism_nuopc_cap.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism_nuopc_cap

  use esmf
  use nuopc
  use NUOPC_Model, &
    model_routine_SS           => SetServices, &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_SetClock       => label_SetClock, &
    model_label_Advance        => label_Advance


  use schism_bmi
  use schism_esmf_util
  use schism_nuopc_util

  implicit none

  private
  public SetServices

contains

#undef ESMF_METHOD
#define ESMF_METHOD "SetServices"
subroutine SetServices(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4)             :: localrc

  character(len=ESMF_MAXSTR), parameter :: label_InternalState = 'InternalState'
  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()

  rc = ESMF_SUCCESS

  call NUOPC_CompDerive(comp, model_routine_SS, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(internalState%wrap, stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_UserCompSetInternalState(comp, label_InternalState, &
    internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! NUOPC automatically has an entry point for InitializeP0, so do not?
  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p0"/), userRoutine=InitializeP0, rc=localrc)

  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeAdvertise, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeRealize, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_SetClock, &
    specRoutine=SetClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_DataInitialize, &
    specRoutine=DataInitialize, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_Advance, &
    specRoutine=ModelAdvance, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Do we need a specialization of Finalize, by adding a label?
  !call NUOPC_CompSpecialize(comp, specLabel=model_label_Finalize, &
  !  specRoutine=Finalize, rc=localrc)
  ! Yes, we should release the haloHandle
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine SetServices

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeP0"
subroutine InitializeP0(comp, importState, exportState, parentClock, rc)

!  phase 0: Set Initialize Phase Definition Version (REQUIRED, NUOPC PROVIDED)
!  – Initialize the InitializePhaseMap Attribute according to the NUOPC Initialize Phase Definition (IPD) version 00
! ∗ IPDv00p1 = 1: (REQUIRED, IMPLEMENTOR PROVIDED) · Advertise Fields in import and export States.
! ∗ IPDv00p2 = 2: (REQUIRED, IMPLEMENTOR PROVIDED)
! · Realize the advertised Fields in import and export States.
! ∗ IPDv00p3 = 3: (REQUIRED, NUOPC PROVIDED)
! · Check compatibility of the Fields’ Connected status.
! ∗ IPDv00p4 = 4: (REQUIRED, NUOPC PROVIDED)
! · Handle Field data initialization. Time stamp the export Fields.

  type(ESMF_GridComp)   :: comp
  type(ESMF_State)      :: importState
  type(ESMF_State)      :: exportState
  type(ESMF_Clock)      :: parentClock
  integer, intent(out)  :: rc

  character(len=10)           :: InitializePhaseMap(4)
  integer                     :: localrc
  logical                     :: isPresent
  character(len=ESMF_MAXSTR)  :: configFileName, compName
  character(len=ESMF_MAXSTR)  :: message
  type(ESMF_Config)           :: config
  type(ESMF_Info)             :: info

  rc=ESMF_SUCCESS

  call NUOPC_CompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetClock(comp, externalClock=parentClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initializing (p=0) component ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  !> @todo add other IPD versions mappings, i.e. IPDv00p1, IPDv01p1, IPDv02p1,
  !> IPDv03p1, IPDv04p1, IPDv05p1; all map to 1 ?
  InitializePhaseMap = (/"IPDv00p1=1","IPDv00p2=2", &
    "IPDv00p3=3","IPDv00p4=4"/)

  !call ESMF_AttributeAdd(comp, convention="NUOPC", &
  !  purpose="General", &
  !  attrList=(/"InitializePhaseMap"/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_InfoGetFromHost(comp, info=info, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_InfoSet(info, key="NUOPC/General/InitializePhaseMap", &
    values=InitializePhaseMap, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

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

      write(message,'(A)')  trim(compName)//' read configuration from '// trim(configFileName)
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    else
      write(message,'(A)')  trim(compName)//' has no configuration; use global config'
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    endif

    call ESMF_GridCompSet(comp, config=config, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

end subroutine InitializeP0

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeAdvertise"
!> @description The purpose of this phase is for your model to advertise its import and export 
!> fields. This means that your model announces which model variables it is capable of exporting
!> (e.g., an atmosphere might export air pressure at sea level) and which model variables it
!> requires (e.g., an atmosphere might require sea surface temperature as a boundary condition).
!> The reason there is an explicit advertise phase is because NUOPC dynamically matches fields
!> among all the models participating in a coupled simulation during runtime. So, we need to
!> collect the list of possible input and output fields from all the models during their initialization.
!> Note that NUOPC does not allocate memory for fields during the advertise phase or when
!> NUOPC_Advertise is called. Instead, this is simply a way for models to communicate the standard
!> names of fields. During a later phase, only those fields that are connected (e.g., a field
!> exported from one model that is imported by another) need to have memory allocated. Also, since
!> ESMF will accept pointers to pre-allocated memory, it is usually not necessary to change
!> how memory is allocated for your model’s variables.
!> Should be mapped to IPDv00p1, IPDv01p1, IPDv02p1, IPDv03p1, IPDv04p1, IPDv05p1;
!> its implementation is required
subroutine InitializeAdvertise(comp, importState, exportState, clock, rc)

  implicit none

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4)       :: localrc, mpiCommunicator, mpiCommDuplicate
  integer(ESMF_KIND_I4)       :: ntime=0, iths=0, i
  character(len=ESMF_MAXSTR)  :: message, compName, cvalue
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
  logical                     :: isPresent, isSet

  type(ESMF_VM)                     :: vm
  type(type_InternalStateWrapper)   :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()

  rc = ESMF_SUCCESS
  localrc = ESMF_SUCCESS

  !> more possikeywordds lities for interface: verbosity, profiling, diagnostic
  call NUOPC_CompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initialize (p=1) component ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  allocate(internalState%wrap, stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  isDataPtr => internalState%wrap
  isDataPtr%numOwnedNodes = 0
  isDataPtr%numForeignNodes = 0

  call SCHISM_InitializePtrMap(comp, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (.not.ESMF_StateIsCreated(importState)) then
    importState=ESMF_StateCreate(name=trim(compName)//'Import', stateintent= &
      ESMF_STATEINTENT_IMPORT, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message,'(A)') trim(compName)//' created state "'//trim(compName)// &
      'Import" for import'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  endif

  if (.not.ESMF_StateIsCreated(exportState)) then
    exportState=ESMF_StateCreate(name=trim(compName)//'Export', stateintent= &
      ESMF_STATEINTENT_EXPORT, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message,'(A)') trim(compName)//' created state "'//trim(compName)// &
      'Export" for export'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  endif

  ! Get VM for this component
  call ESMF_GridCompGet(comp, vm=vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_VMGet(vm, mpiCommunicator=mpiCommunicator, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

#ifndef ESMF_MPIUNI
  call MPI_Comm_dup(mpiCommunicator, mpiCommDuplicate, rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initializing parallel environment ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  call schism_parallel_init(communicator=mpiCommDuplicate)

  write(message, '(A)') trim(compName)//' initialized parallel environment.'
#endif

  call ESMF_UtilIOMkDir ('./outputs',  relaxedFlag=.true., rc=localrc)
  write(message, '(A)') trim(compName)//' writes results to directory "./outputs".'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  write(message, '(A)') trim(compName)//' initializing science model ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  inquire(file='param.nml', exist=isPresent)
  if (.not.isPresent) then
    write(message, '(A)') trim(compName)//' cannot start without required file "param.nml".'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_FILE_OPEN
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  ! query attributes
  call NUOPC_CompAttributeGet(comp, name='meshloc', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  if (isPresent .and. isSet) then
    if (trim(cvalue) == 'node') then
      meshloc = ESMF_MESHLOC_NODE
    else
      meshloc = ESMF_MESHLOC_ELEMENT
    end if
  else
    cvalue = 'node'
    meshloc = ESMF_MESHLOC_NODE
  end if
  write(message, '(A)') trim(compName)//' meshloc is set to '//trim(cvalue)
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! debug option
  call NUOPC_CompAttributeGet(comp, name='debug_level', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  debug_level = 0
  if (isPresent .and. isSet) then
     read(cvalue,*) debug_level
  end if
  write(message, '(A,I1)') trim(compName)//' debug_level is set to ', debug_level
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! init schism
  call schism_init(0, './', iths, ntime)
  write(message, '(A)') trim(compName)//' initialized science model'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! Deal with setting up the Field dictionary here and advertising the
  ! variables.
  ! @todo this should be customized by a field dictionary-like configuration
  ! file

  call NUOPC_FieldDictionaryAddIfNeeded("air_pressure_at_sea_level", "N m-2", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("surface_downwelling_photosynthetic_radiative_flux", "W m-2 s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("sea_surface_temperature", "K", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("inst_temp_height2m", "K", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("sea_surface_salinity", "PSU", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("inst_merid_wind_height10m", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("inst_zonal_wind_height10m", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("eastward_wave_radiation_stress", "N m-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("eastward_northward_wave_radiation_stress", "N m-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("northward_wave_radiation_stress", "N m-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("depth-averaged_x-velocity", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("depth-averaged_y-velocity", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("ocn_current_zonal", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("ocn_current_merid", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("sea_surface_height_above_sea_level", "m", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("ocean_mask", "1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! for coupling to ATM/DATM
  call NUOPC_Advertise(importState, "air_pressure_at_sea_level", rc=localrc)
  call NUOPC_Advertise(importState, "inst_zonal_wind_height10m", rc=localrc)
  call NUOPC_Advertise(importState, "inst_merid_wind_height10m", rc=localrc)
  call NUOPC_Advertise(importState, "inst_temp_height2m", rc=localrc)
  call NUOPC_Advertise(importState, "inst_spec_humid_height2m", rc=localrc)
  call NUOPC_Advertise(importState, "inst_net_sw_flx", rc=localrc)
  call NUOPC_Advertise(importState, "inst_down_lw_flx", rc=localrc)
  call NUOPC_Advertise(importState, "inst_prec_rate", rc=localrc)

  ! for coupling to WW3/WDAT
  call NUOPC_Advertise(importState, "sea_surface_wave_significant_height", rc=localrc)
  call NUOPC_Advertise(importState, "sea_water_waves_effect_on_currents_bernoulli_head_adjustment", rc=localrc)
  call NUOPC_Advertise(importState, "sea_surface_x_stress_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_surface_y_stress_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_bottom_upward_x_stress_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_bottom_upward_y_stress_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_bed_orbital_x_velocity_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_bed_orbital_y_velocity_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_surface_wave_mean_direction", rc=localrc)
  call NUOPC_Advertise(importState, "sea_surface_wave_mean_period", rc=localrc)
  call NUOPC_Advertise(importState, "sea_surface_wave_mean_number", rc=localrc)
  call NUOPC_Advertise(importState, "eastward_surface_stokes_drift_current", rc=localrc)
  call NUOPC_Advertise(importState, "northward_surface_stokes_drift_current", rc=localrc)
  call NUOPC_Advertise(importState, "eastward_wave_radiation_stress", rc=localrc)
  call NUOPC_Advertise(importState, "eastward_northward_wave_radiation_stress", rc=localrc)
  call NUOPC_Advertise(importState, "northward_wave_radiation_stress", rc=localrc)

  !> The mesh information is usually not in CF standard and therefore needs
  !> to be added to the FieldDictionary before advertising
  allocate(itemNameList(4))
  itemNameList=(/ 'mesh_topology                 ', &
                  'mesh_global_node_id           ', &
                  'mesh_global_element_id        ', &
                  'mesh_element_node_connectivity'/)

  do i=1,4
    call NUOPC_FieldDictionaryAddIfNeeded(trim(itemNameList(i)), "1", localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_Advertise(exportState, StandardName=trim(itemNameList(i)), &
      SharePolicyField='share', SharePolicyGeomObject='share', &
      TransferOfferGeomObject='will provide', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  enddo

  ! for coupling to WW3
  call NUOPC_FieldAdvertise(exportState, "ocn_current_zonal", "m s-1", localrc)
  call NUOPC_FieldAdvertise(exportState, "ocn_current_merid", "m s-1", localrc)

  call NUOPC_FieldAdvertise(exportState, "sea_surface_temperature", "K", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "temperature", "K", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "sea_surface_salinity", "PSU", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "depth-averaged_x-velocity", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "depth-averaged_y-velocity", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "sea_surface_height_above_sea_level", "m", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! required for coupling through CMEPS mediator
  call NUOPC_FieldAdvertise(exportState, "ocean_mask", "1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeRealize"
!> @description During this phase, fields that were previously advertised should now be realized. Realizing a field means that an ESMF_Field object is created and it is added to the appropriate ESMF_State, either import or export. In order to create an ESMF_Field, you’ll first need to create one of the ESMF geometric types, ESMF_Grid, ESMF_Mesh, or ESMF_LocStream. For 2D and 3D logically rectangular grids (such as a lat-lon grid), the typical choice is ESMF_Grid. For unstructured grids, use an ESMF_Mesh. Fields are put into import or export States by calling NUOPC_Realize.
!> Should be mapped to IPDv00p2, IPDv01p3, IPDv02p3, IPDv03p3, IPDv04p3, IPDv05p4;
!> Required if providing any geometry.  For higher phases IPDv03p5, IPDv04p5, IPDv05p6
!> a separate routine should be provided for accepting fields.
subroutine InitializeRealize(comp, importState, exportState, clock, rc)

  use schism_esmf_util, only : SCHISM_MeshCreateNode
  use schism_esmf_util, only : SCHISM_MeshCreateElement
  
  !> @todo move all use statements of schism into schism_bmi
  use schism_glbl, only: np, pr2, windx2, windy2, srad, nws, rkind
  use schism_esmf_util, only: SCHISM_StateFieldCreateRealize
  implicit none

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  type(ESMF_TimeInterval) :: stabilityTimeStep
  type(ESMF_Field)        :: field
  integer(ESMF_KIND_I4)   :: localrc, i

  type(ESMF_CoordSys_Flag) :: coordsys
  type(ESMF_Mesh)          :: mesh2d

  character(len=ESMF_MAXSTR)              :: message, compName
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
  type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
  integer(ESMF_KIND_I4)                   :: itemCount

  type(type_InternalStateWrapper)    :: internalState
  type(type_InternalState), pointer  :: isDataPtr => null()
  type(ESMF_DistGrid)                :: nodalDistgrid
  type(ESMF_Array)                   :: array

  real(ESMF_KIND_R8), pointer :: farrayPtr1(:) => null()

  !> @todo move to internal state
  !> Maybe think more generally about how to handle these intermediate varialbes, 
  !> best of course dealt with in a mediator
  real(ESMF_KIND_R8), target, allocatable :: eastward_wave_radiation_stress(:)
  real(ESMF_KIND_R8), target, allocatable :: eastward_northward_wave_radiation_stress(:)
  real(ESMF_KIND_R8), target, allocatable :: northward_wave_radiation_stress(:)

  rc = ESMF_SUCCESS
  localrc= ESMF_SUCCESS

!<<<<<<< HEAD
!=======
  !> @todo move addSchismMesh back to schism_esmf_util to share with ESMF cap
  !> call addSchismMesh(comp, ownedNodes=ownedNodes, foreignNodes=foreignNodes, rc=localrc)
  !call addSchismMesh(comp, localrc)
!>>>>>>> a9a0ce0 (Use MeshCreate instead off add Mesh)
  if (meshloc == ESMF_MESHLOC_NODE) then
    call SCHISM_MeshCreateNode(comp, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  else
    call SCHISM_MeshCreateElement(comp, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end if

  call ESMF_GridCompGet(comp, mesh=mesh2d, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_MeshGet(mesh2d, nodalDistgrid=nodalDistgrid, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_zonal_wind_height10m", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_merid_wind_height10m", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="air_pressure_at_sea_level", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="downwelling_short_photosynthetic_radiation_at_water_surface", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_wave_significant_height", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_water_waves_effect_on_currents_bernoulli_head_adjustment", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_x_stress_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_y_stress_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_bottom_upward_x_stress_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_bottom_upward_y_stress_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_bed_orbital_x_velocity_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_bed_orbital_y_velocity_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_wave_mean_direction", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_wave_mean_period", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_wave_mean_number", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="eastward_surface_stokes_drift_current", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="northward_surface_stokes_drift_current", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo add more atmospheric fields (like humidity)
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_temp_height2m", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_spec_humid_height2m", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_net_sw_flx", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_down_lw_flx", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_prec_rate", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Wave parameters, for now we only have those from the WW3DATA cap in 
  !> NOAA's CoastalApp.  
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="eastward_wave_radiation_stress", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="eastward_northward_wave_radiation_stress", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="northward_wave_radiation_stress", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> The list of export states is declared in InitializeAdvertise
  call ESMF_StateGet(exportState, itemCount=itemCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(itemNameList(itemCount), itemTypeList(itemCount), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateGet(exportState, itemTypeList=itemTypeList,  &
    itemNameList=itemNameList, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i=1, itemCount

    if (itemTypeList(i) /= ESMF_STATEITEM_FIELD) cycle

    call SCHISM_FieldRealize(exportState, itemNameList(i), &
      mesh=mesh2d, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message,'(A)') trim(compName)//' realized field '//trim(itemNameList(i))// &
      ' in its export state'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  enddo

  if (allocated(itemNameList)) deallocate(itemNameList)
  if (allocated(itemTypeList)) deallocate(itemTypeList)

  call ESMF_StateGet(importState, itemCount=itemCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (itemCount > 0) then 
    allocate(itemNameList(itemCount), itemTypeList(itemCount), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateGet(importState, itemTypeList=itemTypeList,  &
      itemNameList=itemNameList, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  do i=1, itemCount

    if (itemTypeList(i) /= ESMF_STATEITEM_FIELD) cycle

    write(message,'(A)') trim(compName)//' realized field '//trim(itemNameList(i))// &
      ' in its import state'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  enddo

  
  !> As IPDv01p3 (NUOPC Provided) fails when fields are not connected, we here
  !> remove all unconnected fields from the import and export stabilityTimeStep
  call SCHISM_RemoveUnconnectedFields(importState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_RemoveUnconnectedFields(exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  !@TODO: should we destroy field & array vars?
end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "SetClock"
!> @description A model's clock copies start/stop time and timestep from its
!> parent's clock.  If a model has a timestep that is different (smaller) than
!> the parent's it needs to be set here.
subroutine SetClock(comp, rc)

  use schism_glbl, only : start_year, start_month, start_day, start_hour, rnday
  use schism_glbl, only : wtiminc, wtime2

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Clock)           :: driverClock, modelClock
  type(ESMF_TimeInterval)    :: runDur, timeStep
  type(ESMF_Time)            :: startTime, stopTime
  integer                    :: localrc, d, h, m
  real(ESMF_KIND_R8)         :: seconds
  character(len=ESMF_MAXSTR) :: message

  rc = ESMF_SUCCESS

  call NUOPC_ModelGet(comp, driverClock=driverClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Set start time
  call ESMF_TimeSet(startTime, yy=start_year, mm=start_month, dd=start_day, h=int(start_hour), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Set stop time
  d = int(rnday)
  h = int((rnday-d)*24)
  m = int((rnday-d)*24*60-h*60)
  call ESMF_TimeIntervalSet(runDur, d=d, h=h, m=m, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  stopTime = startTime+runDur

  !> Time step must be same with the driver
  call ESMF_ClockGet(driverClock, timeStep=timeStep, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Create component clock
  modelClock = ESMF_ClockCreate(timeStep, startTime, stopTime=stopTime, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Update component clock
  call ESMF_GridCompSet(comp, clock=modelClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Check that wtiminc, i.e. the time between two new atmospheric inputs
  !> corresponds to the parent (coupling) time step
  call ESMF_TimeIntervalGet(timeStep, s_r8=seconds, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (abs(wtiminc - seconds) > 1e-6) then 
    write(message, '(A,I7,A,I7)') 'Check setting of wtiminc = ', int(wtiminc), &
      ' in param.nml! Auto-resetting to wtiminc = ', int(seconds)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif

  wtiminc = seconds
  wtime2 = wtiminc !Also need to reset wtime2 to make it consistent

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "SetRunClock"
!> @description If the timestep of the parent is dynamic, then there might
!> be mismatches between intergral timesteps of a model and the parent
!> timestep.  This is reported as an error unless dealt with in this label
!> @todo (not registered and fully implemented yet, so not used)
subroutine SetRunClock(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Clock) :: driverClock, modelClock
  type(ESMF_Time)  :: driverCurrTime
  integer          :: localrc

  rc = ESMF_SUCCESS

  !> query driver and the component clocks
  call NUOPC_ModelGet(comp, driverClock=driverClock, modelClock=modelClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockGet(driverClock, currTime=driverCurrTime, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> set model clock to have the current start time as the driver clock
  call ESMF_ClockSet(modelClock, currTime=driverCurrTime, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> check the component clock against the driver clock
  call NUOPC_CompCheckSetClock(comp, driverClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "DataInitialize"
subroutine DataInitialize(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  ! local variables
  type(ESMF_Time) :: currTime
  type(ESMF_Clock) :: clock
  type(ESMF_State) :: exportState
  integer(ESMF_KIND_I4) :: localrc
  character(len=*), parameter :: subname = '(DataInitialize): '
  !--------------------------------

  rc = ESMF_SUCCESS
  call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

  !> Query component for its clock, import and export states
  call NUOPC_ModelGet(comp, modelClock=clock, exportState=exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Update fields on export state
  call SCHISM_Export(comp, exportState, clock, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

end subroutine DataInitialize

#undef ESMF_METHOD
#define ESMF_METHOD "ModelAdvance"
!> @description As described in the section 4.2, the subroutine ModelAdvance (shown below) has been registered to the special- ization point with the label model_label_Advance in the SetServices subroutine. This specialization point subroutine is called within the generic NUOPC_Model run phase in order to request that your model take a timestep forward. The code to do this is model dependent, so it does not appear in the subroutine below. Each NUOPC component maintains its own clock (an ESMF_Clock object). The clock is used to indicate the current model time and the timestep size. When the subroutine finishes, your model should be moved ahead in time from the current time by one timestep. NUOPC will automatically advance the clock for you, so there is no explicit call to do that here.
!> Because the import/export states and the clock do not come in through the parameter list, they must be accessed via a call to NUOPC_ModelGet
subroutine ModelAdvance(comp, rc)

  use schism_glbl, only: dt

  implicit none

  !> Input/output variables
  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  !> Local variables 
  type(ESMF_Clock)        :: clock
  type(ESMF_State)        :: importState, exportState
  type(ESMF_Time)         :: currTime
  type(ESMF_TimeInterval) :: timeStep
  character(len=160)      :: message
  integer(ESMF_KIND_I4)   :: localrc
  integer, save           :: it = 1
  integer                 :: i, num_schism_steps
  real(ESMF_KIND_R8)      :: seconds
  character(len=*), parameter :: subname = '(ModelAdvance): '
  !--------------------------------

  rc = ESMF_SUCCESS
  call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

  !> Query component for its clock, import and export states
  call NUOPC_ModelGet(comp, modelClock=clock, importState=importState, &
    exportState=exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Update fields on import state
  call SCHISM_Import(comp, importState, clock, rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Write log about advance
  call ESMF_ClockPrint(clock, options="currTime", &
      preString="--- advancing schism from ", unit=message, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_TimePrint(currTime + timeStep, &
      preString=trim(message)//" to ", unit=message, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogWrite(message, ESMF_LOGMSG_INFO, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Run SCHISM
  call ESMF_TimeIntervalGet(timeStep, s_r8=seconds, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  num_schism_steps=int(seconds/dt) 
  if (mod(seconds,dt) /= 0) then
          write(message, '(A)') 'Coupling step cannot be divided by dt, please adjust to avoid lack of steps! '
     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  end if

  do i = it, it+num_schism_steps-1
     call schism_step(i)
     it = it + 1
  end do

  !> Update fields on export state
  call SCHISM_Export(comp, exportState, clock, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

end subroutine ModelAdvance

#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_RemoveUnconnectedFields"
subroutine SCHISM_RemoveUnconnectedFields(state, rc)

  implicit none

  type(ESMF_State), intent(inout)      :: state
  integer(kind=ESMF_KIND_I4), optional :: rc

  integer(kind=ESMF_KIND_I4)              :: rc_, localrc, itemCount, i
  type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
  character(len=ESMF_MAXSTR)              :: message
  type(ESMF_Field)                        :: field
  logical                                 :: isPresent

  if (present(rc)) rc = ESMF_SUCCESS

  call ESMF_StateGet(state, itemCount=itemCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  allocate(itemTypeList(itemCount))
  allocate(itemNameList(itemCount))

  call ESMF_StateGet(state, itemTypeList=itemTypeList,  &
    itemNameList=itemNameList, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  do i=1, itemCount

    if (itemTypeList(i) /= ESMF_STATEITEM_FIELD) cycle

    call ESMF_StateGet(state, trim(itemNameList(i)), field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    isPresent=.true.
!    call ESMF_AttributeGet(field, name='Connected', isPresent=isPresent, rc=localrc)
!    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A)') '--- checking connection state of '//trim(itemNameList(i))
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    !if (isPresent) isPresent = NUOPC_IsConnected(state,trim(itemNameList(i)), rc=localrc)
    if (isPresent) isPresent = NUOPC_IsConnected(field, rc=localrc)
    if(localrc/=ESMF_SUCCESS) then
      isPresent=.false.
      localrc=ESMF_SUCCESS
    endif
    
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    if (.not.isPresent) then
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

      call ESMF_StateRemove(state, itemNameList(i:i), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

      write(message,'(A)') '--- removed unconnected field '//trim(itemNameList(i))
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    endif

  enddo
end subroutine SCHISM_RemoveUnconnectedFields

#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_Export"
subroutine SCHISM_Export(comp, exportState, clock, rc)

  use schism_glbl,      only: nvrt, eta2, dav, uu2, vv2, tr_nd, idry_e, npa
  use schism_esmf_util, only: SCHISM_StateUpdate

  implicit none

  !> Input/output variables
  type(ESMF_GridComp), intent(in)    :: comp
  type(ESMF_State)   , intent(inout) :: exportState
  type(ESMF_Clock)   , intent(in)    :: clock
  integer            , intent(inout) :: rc

  !> Local variables
  type(ESMF_Time) :: currTime
  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()
  integer(ESMF_KIND_I4) :: localrc
  real(ESMF_KIND_R8), allocatable, save, target :: idry_r8(:)
  real(ESMF_KIND_R8), allocatable, save, target :: sst_K(:)
  character(len=ESMF_MAXSTR) :: timeStr
  character(len=*), parameter :: subname = '(SCHISM_Export): '
  !--------------------------------

  rc = ESMF_SUCCESS
  call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

  !> Query internal state
  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  isDataPtr => internalState%wrap

  if(.not.associated(isDataPtr)) localrc = ESMF_RC_PTR_NULL
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Update fields on export state
  !> sea surface height
  call SCHISM_StateUpdate(exportState, 'sea_surface_height_above_sea_level', eta2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> depth average current in x direction
  call SCHISM_StateUpdate(exportState, 'depth-averaged_x-velocity', dav(1,:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> depth average current in y direction
  call SCHISM_StateUpdate(exportState, 'depth-averaged_y-velocity', dav(2,:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> surface current in x direction
  call SCHISM_StateUpdate(exportState, 'ocn_current_zonal', uu2(nvrt,:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> surface current in y direction
  call SCHISM_StateUpdate(exportState, 'ocn_current_merid', vv2(nvrt,:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> surface temperature
  if (.not. allocated(sst_K)) then
     allocate(sst_K(npa))
     sst_K(:) =0.0d0+273.15d0
  end if
  sst_K(:) = tr_nd(1,nvrt,:)+273.15d0 !Change unit to K
  call SCHISM_StateUpdate(exportState, 'sea_surface_temperature', sst_K(:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> surface salinity
  call SCHISM_StateUpdate(exportState, 'sea_surface_salinity', tr_nd(2,nvrt,:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> ocean mask
  !> mediator expects ocean mask in double rather then integer
  if (.not. allocated(idry_r8)) then
     allocate(idry_r8(size(idry_e)))
     idry_r8(:) =0.0d0
  end if
  idry_r8(:) = dble(idry_e(:))

  call SCHISM_StateUpdate(exportState, 'ocean_mask', idry_r8, &
    isPtr=isDataPtr, onElement=.true., rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Write fields on export state for debugging
  if (debug_level > 0) then
     call ESMF_ClockGet(clock, currTime=currTime, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

     call ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=localrc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

     call SCHISM_StateWriteVTK(exportState, 'export_'//trim(timeStr), rc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end if

  call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

end subroutine SCHISM_Export

#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_Import"
subroutine SCHISM_Import(comp, importState, clock, rc)

  use schism_glbl     , only: RADFLAG, windx2, windy2, pr2
  use schism_glbl     , only: airt2,shum2,srad,hradd,fluxprc
  use schism_esmf_util, only: SCHISM_StateImportWaveTensor
  use schism_esmf_util, only: SCHISM_StateImportWave3dVortex
  use schism_esmf_util, only: SCHISM_StateUpdate

  implicit none

  !> Input/output variables
  type(ESMF_GridComp), intent(in)    :: comp
  type(ESMF_State)   , intent(inout) :: importState
  type(ESMF_Clock)   , intent(in)    :: clock
  integer            , intent(inout) :: rc

  !> Local variables
  type(ESMF_Time) :: currTime
  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()
  integer(ESMF_KIND_I4) :: localrc
  character(len=ESMF_MAXSTR) :: timeStr
  character(len=*), parameter :: subname = '(SCHISM_Import): '
  !--------------------------------

  rc = ESMF_SUCCESS
  call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

  !> Query internal state
  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  isDataPtr => internalState%wrap

  if(.not.associated(isDataPtr)) localrc = ESMF_RC_PTR_NULL
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Update fields on import state
  if (RADFLAG == 'VOR') then
     !> Obtain required variables from wave component to do coupling with vortex formulation
     call SCHISM_StateImportWave3dVortex(importState, isDataPtr, rc=localrc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  else
     !> Obtain radiation tensor from wave component and calculate the wave stress
     call SCHISM_StateImportWaveTensor(importState, isDataPtr, rc=localrc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end if

  !> surface wind component in x direction
  call SCHISM_StateUpdate(importState, 'inst_zonal_wind_height10m', windx2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Surface wind component in y direction
  call SCHISM_StateUpdate(importState, 'inst_merid_wind_height10m', windy2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Surface air pressure
  call SCHISM_StateUpdate(importState, 'air_pressure_at_sea_level', pr2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Air temperature at sea level
  call SCHISM_StateUpdate(importState, 'inst_temp_height2m', airt2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Air humidity at sea level
  call SCHISM_StateUpdate(importState, 'inst_spec_humid_height2m', shum2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Downward shortwave at sea level
  call SCHISM_StateUpdate(importState, 'inst_net_sw_flx', srad, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Downward longwave at sea level
  call SCHISM_StateUpdate(importState, 'inst_down_lw_flx', hradd, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Precipatation rate
  call SCHISM_StateUpdate(importState, 'inst_prec_rate', fluxprc, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Write fields on import state for debugging
  if (debug_level > 0) then
     call ESMF_ClockGet(clock, currTime=currTime, rc=localrc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

     call ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=localrc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

     call SCHISM_StateWriteVTK(importState, 'import_'//trim(timeStr), rc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end if

  call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

end subroutine SCHISM_Import

#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateWriteVTK"
subroutine SCHISM_StateWriteVTK(state, prefix, rc)

  implicit none

  !> Input/output variables
  type(ESMF_State), intent(in) :: state
  character(len=*), intent(in) :: prefix
  integer, intent(out), optional :: rc

  !> local variables
  integer :: i, itemCount
  type(ESMF_Field) :: field
  character(ESMF_MAXSTR), allocatable :: itemNameList(:)
  character(len=*),parameter :: subname='(SCHISM_StateWriteVTK)'
  !--------------------------------

  rc = ESMF_SUCCESS
  call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

  !> Get number of fields in the state
  call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  !> Get item names
  if (.not. allocated(itemNameList)) allocate(itemNameList(itemCount))

  call ESMF_StateGet(state, itemNameList=itemNameList, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  !> Loop over fields and write them
  do i = 1, itemCount
     !> Get field
     call ESMF_StateGet(state, itemName=trim(itemNameList(i)), field=field, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

     !> Write field
     call ESMF_FieldWriteVTK(field, trim(prefix)//'_'//trim(itemNameList(i)), rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
  end do

  !> Clean temporary variables
  if (allocated(itemNameList)) deallocate(itemNameList)

  call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

end subroutine SCHISM_StateWriteVTK

end module schism_nuopc_cap
