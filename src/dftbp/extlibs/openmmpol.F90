!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"

module dftbp_extlibs_openmmpol

#:if WITH_OPENMMPOL
    use ommp_interface, only : ommp_system, ommp_init_mmp, ommp_terminate, ommp_print_summary, &
                               ommp_prepare_qm_ele_ene, ommp_set_external_field, ommp_potential_mmpol2ext, &
                               ommp_qm_helper, ommp_init_qm_helper, ommp_terminate_qm_helper, ommp_set_verbose, &
                               OMMP_VERBOSE_DEBUG, OMMP_VERBOSE_LOW, OMMP_VERBOSE_HIGH

    use dftbp_common_environment, only : TEnvironment
    use dftbp_dftb_periodic, only : TNeighbourList
    use dftbp_type_commontypes, only : TOrbitals
    use dftbp_dftb_charges, only : getSummedCharges
#:endif
    use dftbp_common_accuracy, only : dp
    implicit none

! SECTION: type definitions

#:if WITH_OPENMMPOL
    type :: TOMMPInterface
        type(ommp_system), pointer :: pSystem
        type(ommp_qm_helper), pointer :: pQMHelper
        integer :: solver
        real(dp), allocatable :: qmmmCouplingEnergyPerAtom(:)
        real(dp) :: forceFieldEnergy
    contains
        ! procedure :: updateQMCoords
        procedure :: updateQMCharges
        procedure :: getEnergies
        procedure :: addPotential
        ! procedure :: addGradients
    end type
#:else
    type :: TOMMPInterface
    end type
#:endif

type :: TOMMPInput
    character(:), allocatable :: filename
    integer :: solver
end type

public TOMMPInterface, TOMMPInterface_init

contains 

! SECTION: subroutines
#:if WITH_OPENMMPOL
        subroutine TOMMPInterface_init(this, openmmpolInput, nQMatoms, atomTypes, qmAtomCoords)
            type(TOMMPInterface), intent(out) :: this
            type(TOMMPInput), intent(in) :: openmmpolInput
            integer, intent(in) :: nQMatoms
            integer, dimension(nQMatoms), intent(in) :: atomTypes
            real(dp), dimension(3, nQMatoms), intent(in) :: qmAtomCoords
            real(dp), allocatable, dimension(:) :: netCharges

            this%solver = openmmpolInput%solver
            this%forceFieldEnergy = 0.0_dp

            call ommp_init_mmp(this%pSystem, openmmpolInput%filename)
            ! call ommp_set_verbose(OMMP_VERBOSE_DEBUG)
            allocate(netCharges(nQMatoms))
            call ommp_init_qm_helper(this%pQMHelper, nQMatoms, qmAtomCoords, netCharges, atomTypes)
            deallocate(netCharges)

        end subroutine

        subroutine TOMMPInterface_terminate(this)
            type(TOMMPInterface), intent(out) :: this
            call ommp_terminate(this%pSystem)
            call ommp_terminate_qm_helper(this%pQMHelper)
        end subroutine

        subroutine updateQMCoords(this)
            class(TOMMPInterface), intent(inout) :: this

            !> TODO: to be cleared
            this%pSystem%eel%ipd_done = .false.
            this%pSystem%eel%D2mgg_done = .false.
            this%pSystem%eel%D2dgg_done = .false.
            this%pQMHelper%E_n2p_done = .false.
        end subroutine

        subroutine updateQMCharges(this, env, species, neighList, qq, q0, img2CentCell, orb)
            class(TOMMPInterface), intent(inout) :: this

            !> Computational environment settings
            type(TEnvironment), intent(in) :: env

            !> Species, shape: [nAtom]
            integer, intent(in) :: species(:)
      
            !> Neighbour list.
            type(TNeighbourList), intent(in) :: neighList
      
            !> Orbital charges.
            real(dp), intent(in) :: qq(:,:,:)
      
            !> Reference orbital charges.
            real(dp), intent(in) :: q0(:,:,:)
      
            !> Mapping on atoms in central cell.
            integer, intent(in) :: img2CentCell(:)
      
            !> Orbital information
            type(TOrbitals), intent(in) :: orb

            ! Charge per atom (internal variable)
            real(dp), allocatable :: qPerAtom(:)

            write(*, *) "Call to QM charges"
            allocate(qPerAtom(size(species)))

            !> getSummedCharges computes population, necessary to multiply by -1 to get charge
            call getSummedCharges(species, orb, qq, q0=q0, dQatom=qPerAtom)
            qPerAtom = -qPerAtom

            !> Set charges
            this%pQMHelper%qqm = qPerAtom 
            deallocate(qPerAtom)

            !> Charges updated, re-evaluation of quantities is requested
            this%pSystem%eel%ipd_done = .false.
            this%pSystem%eel%D2mgg_done = .false.
            this%pSystem%eel%D2dgg_done = .false.
            this%pQMHelper%E_n2p_done = .false.

            ! > Compute electric field produced by QM part of the system on MM atoms
            call ommp_prepare_qm_ele_ene(this%pSystem, this%pQMHelper)

            !> Set external field for MM, solve the polarization equations
            call ommp_set_external_field(this%pSystem, this%pQMHelper%E_n2p, this%solver, .true.)

            write(*, "(A,F12.6)") "E_QMMM: ", dot_product(this%pQMHelper%qqm, this%pQMHelper%V_m2n)
        end subroutine

        subroutine addPotential(this, shiftPerAtom)
            class(TOMMPInterface) :: this

            !> Computed potential is added to this vector
            real(dp), intent(inout) :: shiftPerAtom(:)
            
            !> Computed shift
            real(dp), allocatable :: shiftFromOpenmmpol(:)

            !> TODO: debug
            write(*, *) "Call to addPotential"
            write(*, *) this%pQMHelper%qqm

            !> Add potential from openmmpol to the vector of potentials
            shiftPerAtom = shiftPerAtom + this%pQMHelper%V_m2n

        end subroutine

        subroutine getEnergies(this)
            class(TOMMPInterface) :: this

        end subroutine

#:else
    subroutine TOMMPInterface_init()
    end subroutine

    subroutine TOMMPInterface_terminate()
    end subroutine
#:endif


end module
