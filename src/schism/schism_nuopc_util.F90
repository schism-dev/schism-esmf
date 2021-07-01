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

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism_nuopc_util

  use esmf
  use NUOPC

  implicit none

  public NUOPC_FieldAdvertise
  private

contains

#undef  ESMF_METHOD
#define ESMF_METHOD "NUOPC_FieldDictionaryAddIfNeeded"
subroutine NUOPC_FieldDictionaryAddIfNeeded(name, unit, rc)

  use NUOPC, only : NUOPC_FieldDictionaryHasEntry, NUOPC_FieldDictionaryAddEntry
  implicit none

  character(len=ESMF_MAXSTR), intent(in)           :: name
  character(len=ESMF_MAXSTR), intent(in)           :: unit
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

  type(ESMF_State), intent(inout)                  :: state
  character(len=ESMF_MAXSTR), intent(in)           :: name
  character(len=ESMF_MAXSTR), intent(in)           :: unit
  integer(ESMF_KIND_I4), intent(out), optional     :: rc

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

end module schism_nuopc_util

