! This code is part of the SCHISM-ESMF interface and defines
! a generic BMI (basic model interface) to schism
!
! @copyright 2021-2022 Helmholtz-Zentrum Hereon
! @copyright 2018-2021 Helmholtz-Zentrum Geesthacht
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
#define ESMF_FILENAME "schism_bmi.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism_bmi

  ! We *should not* use ESMF in a BMI, but this is used here for logging only
  use esmf!, only:: ESMF_LogFoundError, ESMF_END_ABORT, ESMF_Finalize, ESMF_SUCCESS

  interface
    subroutine parallel_init(communicator)
      implicit none
      integer, optional :: communicator
    end subroutine parallel_init

    subroutine parallel_finalize
      implicit none
    end subroutine parallel_finalize

    subroutine schism_init(iorder,indir,iths,ntime)
      implicit none
      integer, intent(in)          :: iorder
      character(len=*), intent(in) :: indir
      integer, intent(out)         :: iths,ntime
    end subroutine schism_init

    subroutine schism_step(it)
      implicit none
      integer, intent(in) :: it
    end subroutine schism_step

    subroutine schism_finalize()
      implicit none
    end subroutine schism_finalize

  end interface

contains

#undef  ESMF_METHOD
#define ESMF_METHOD "schism_parallel_init"
subroutine schism_parallel_init(communicator)

  use schism_msgp, only: parallel_init
  use schism_msgp, only: schism_mpi_comm=>comm

  implicit none
  integer :: communicator

  call parallel_init(communicator)

end subroutine schism_parallel_init

#undef  ESMF_METHOD
#define ESMF_METHOD "schismTimeStep"
subroutine schismTimeStep(seconds)

  use schism_glbl, only: dt

  implicit none
  double precision, intent(out) :: seconds

  seconds = dt

end subroutine schismTimeStep

#undef ESMF_METHOD
#define ESMF_METHOD 'schismPtr1'
function schismPtr1(varname) result(farrayPtr)

  use schism_glbl, only: airt2
  character(len=*), intent(in) :: varname
  double precision, pointer :: farrayPtr(:)

  select case((trim(varname)))
  case ('airt2')
    farrayPtr => airt2
  case default
    nullify(farrayPtr)
  end select

end function schismPtr1

end module schism_bmi
