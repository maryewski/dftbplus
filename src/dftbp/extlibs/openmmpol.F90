!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"

module dftbp_extlibs_openmmpol

#:if WITH_OPENMMPOL
    use ommp_interface, only : ommp_system, ommp_init_mmp, ommp_terminate
#:endif

    implicit none

! SECTION: type definitions

#:if WITH_OPENMMPOL
    type :: TOMMPInterface
        type(ommp_system), pointer :: pSystem
    end type
#:else
    type :: TOMMPInterface
    end type
#:endif

type :: TOMMPInput
    character(:), allocatable :: filename
end type

public TOMMPInterface, TOMMPInterface_init

contains 

! SECTION: subroutines
#:if WITH_OPENMMPOL
        subroutine TOMMPInterface_init(this, openmmpolInput)
            type(TOMMPInterface), intent(out) :: this
            type(TOMMPInput) :: openmmpolInput
            ! character(len=*), intent(in) :: filename

            call ommp_init_mmp(this%pSystem, openmmpolInput%filename)
        end subroutine
#:else
    subroutine TOMMPInterface_init()
    end subroutine
#:endif


end module
