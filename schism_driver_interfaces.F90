! module with interfaces as used in schism_driver

module schism_driver_interfaces

  interface
    subroutine parallel_init(communicator)
      implicit none
      integer, optional :: communicator
    end subroutine parallel_init


    subroutine parallel_finalize
      implicit none
    end subroutine parallel_finalize


    subroutine schism_init(indir,iths,ntime)
      implicit none
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

end module schism_driver_interfaces
