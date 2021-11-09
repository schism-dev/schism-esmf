! This code is part of the SCHISM-ESMF interface, it defines utility
! functions used by the NUOPC cap
!
! @copyright 2021 Helmholtz-Zentrum Hereon
! @author Carsten Lemmen <carsten.lemmen@hereon.de>
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
#define ESMF_FILENAME "schism_esmf_util.F90"

#ifndef _SCHISM_LOG_AND_FINALIZE_ON_ERROR_
#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
#endif

module schism_nuopc_util

  use esmf
  use NUOPC

  implicit none

  public NUOPC_FieldAdvertise, NUOPC_FieldDictionaryAddIfNeeded, SCHISM_StateImportWaveTensor
  private

contains

#undef  ESMF_METHOD
#define ESMF_METHOD "NUOPC_FieldDictionaryAddIfNeeded"
subroutine NUOPC_FieldDictionaryAddIfNeeded(name, unit, rc)

  use NUOPC, only : NUOPC_FieldDictionaryHasEntry, NUOPC_FieldDictionaryAddEntry
  implicit none

  character(len=*), intent(in)           :: name
  character(len=*), intent(in)           :: unit
  integer(ESMF_KIND_I4), intent(out), optional     :: rc

  logical                 :: isPresent
  integer(ESMF_KIND_I4)   :: localrc, rc_
  character(len=ESMF_MAXSTR) :: message

  if (.not.NUOPC_FieldDictionaryHasEntry(trim(name), rc=localrc)) then
    call NUOPC_FieldDictionaryAddEntry(trim(name), trim(unit), rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)
  endif

  if (present(rc)) rc=localrc

end subroutine NUOPC_FieldDictionaryAddIfNeeded

#undef  ESMF_METHOD
#define ESMF_METHOD "NUOPC_FieldAdvertise"
subroutine NUOPC_FieldAdvertise(state, name, unit, rc)

  use NUOPC, only : NUOPC_Advertise
  implicit none

  type(ESMF_State), intent(inout)              :: state
  character(len=*), intent(in)                 :: name
  character(len=*), intent(in)                 :: unit
  integer(ESMF_KIND_I4), intent(out), optional :: rc

  logical                 :: isPresent
  integer(ESMF_KIND_I4)   :: localrc, rc_
  character(len=ESMF_MAXSTR) :: message

  call NUOPC_FieldDictionaryAddIfNeeded(trim(name), trim(unit), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  call NUOPC_Advertise(state, trim(name), trim(unit),  &
    TransferOfferGeomObject='will provide',  rc=localrc)
  !   SharePolicyField='share', SharePolicyGeomObject='not share', &
   _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (present(rc)) rc=localrc

end subroutine NUOPC_FieldAdvertise

#undef  ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateImportWaveTensor""
subroutine SCHISM_StateImportWaveTensor(state, rc)

  use schism_glbl, only: nsa

  implicit none

  type(ESMF_State), intent(in)                 :: state
  integer(ESMF_KIND_I4), intent(out), optional :: rc

  logical                    :: isPresent
  integer(ESMF_KIND_I4)      :: localrc, rc_, i
  character(len=ESMF_MAXSTR) :: message
  type(ESMF_Field)           :: field
  type(ESMF_StateItem_Flag)  :: itemType

  real(ESMF_KIND_R8), pointer :: farrayPtr1(:) => null()
  real(ESMF_KIND_R8), pointer :: radiation_stress_component_sxx(:) => null()
  real(ESMF_KIND_R8), pointer :: radiation_stress_component_sxy(:) => null()
  real(ESMF_KIND_R8), pointer :: radiation_stress_component_syy(:) => null()

  if (present(rc)) rc=localrc

  call ESMF_StateGet(state, "radiation_stress_component_sxx", itemType=itemType, rc=localrc)
  if (itemType /= ESMF_STATEITEM_FIELD) return

  allocate(radiation_stress_component_sxx(nsa))
  call ESMF_StateGet(state, "radiation_stress_component_sxx", field=field, rc=localrc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)

  do i=1,nsa
    radiation_stress_component_sxx(i) = farrayPtr1(i)
  enddo
  
  call ESMF_StateGet(state, "radiation_stress_component_sxy", itemType=itemType, rc=localrc)
  if (itemType /= ESMF_STATEITEM_FIELD) return

  allocate(radiation_stress_component_sxy(nsa))
  call ESMF_StateGet(state, "radiation_stress_component_sxy", field=field, rc=localrc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)

  do i=1,nsa
    radiation_stress_component_sxy(i) = farrayPtr1(i)
  enddo

  call ESMF_StateGet(state, "radiation_stress_component_syy", itemType=itemType, rc=localrc)
  if (itemType /= ESMF_STATEITEM_FIELD) return

  allocate(radiation_stress_component_syy(nsa))
  call ESMF_StateGet(state, "radiation_stress_component_syy", field=field, rc=localrc)

  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)

  do i=1,nsa
    radiation_stress_component_syy(i) = farrayPtr1(i)
  enddo
  nullify(farrayPtr1)

  ! @todo check RSXX etc must have dimnesion of m*m/s/s!
!  call compute_waveforce_from_stress(radiation_stress_component_sxx, & 
!    radiation_stress_component_sxy,radiation_stress_component_syy)
  call compute_wave_force_lon(radiation_stress_component_sxx, & 
    radiation_stress_component_sxy,radiation_stress_component_syy)

  deallocate(radiation_stress_component_sxx)
  deallocate(radiation_stress_component_sxy)
  deallocate(radiation_stress_component_syy)
  
end subroutine SCHISM_StateImportWaveTensor

end module schism_nuopc_util

