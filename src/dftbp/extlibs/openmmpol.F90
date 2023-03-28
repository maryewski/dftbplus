!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"
#:include "error.fypp"

module dftbp_extlibs_openmmpol

#:if WITH_OPENMMPOL
    use ommp_interface, only : ommp_system, ommp_init_mmp, ommp_init_xyz, ommp_terminate, ommp_print_summary, &
                               ommp_prepare_qm_ele_ene, ommp_set_external_field, ommp_potential_mmpol2ext, &
                               ommp_qm_helper, ommp_init_qm_helper, ommp_terminate_qm_helper, ommp_set_verbose, &
                               OMMP_VERBOSE_DEBUG, OMMP_VERBOSE_LOW, OMMP_VERBOSE_HIGH
#:endif
    use dftbp_io_message, only : error, warning
    use dftbp_common_environment, only : TEnvironment
    use dftbp_dftb_periodic, only : TNeighbourList
    use dftbp_type_commontypes, only : TOrbitals
    use dftbp_dftb_charges, only : getSummedCharges
    use dftbp_common_accuracy, only : dp
    implicit none

! SECTION: type definitions
    type :: TOMMPInterface

    #:if WITH_OPENMMPOL
        type(ommp_system), pointer :: pSystem
        type(ommp_qm_helper), pointer :: pQMHelper
    #:endif

        integer :: solver
        real(dp), allocatable :: qmmmCouplingEnergyPerAtom(:)
        real(dp), allocatable :: qmAtomsPotential(:)
        real(dp) :: forceFieldEnergy
    contains
        ! procedure :: updateQMCoords
        procedure :: updateQMCharges
        procedure :: addAtomEnergies
        procedure :: addPotential
        ! procedure :: addGradients
    end type

type :: TOMMPInput
    character(:), allocatable :: inputFormat
    character(:), allocatable :: geomFilename
    character(:), allocatable :: paramsFilename
    integer :: solver
end type

public TOMMPInterface, TOMMPInterface_init

contains 

! SECTION: subroutines
        subroutine TOMMPInterface_init(this, openmmpolInput, nQMatoms, atomTypes, qmAtomCoords)
            type(TOMMPInterface), intent(out) :: this
            type(TOMMPInput), intent(in) :: openmmpolInput
            integer, intent(in) :: nQMatoms
            integer, dimension(nQMatoms), intent(in) :: atomTypes
            real(dp), dimension(3, nQMatoms), intent(in) :: qmAtomCoords
            real(dp), allocatable, dimension(:) :: netCharges
            
            #:if WITH_OPENMMPOL
            this%solver = openmmpolInput%solver
            this%forceFieldEnergy = 0.0_dp
            allocate(this%qmAtomsPotential(nQMatoms))
            this%qmAtomsPotential(:) = 0.0_dp
            allocate(this%qmmmCouplingEnergyPerAtom(nQMatoms))
            this%qmmmCouplingEnergyPerAtom(:) = 0.0_dp

            if (openmmpolInput%inputFormat == "Tinker") then
                call ommp_init_xyz(this%pSystem, openmmpolInput%geomFilename, openmmpolInput%paramsFilename)
            else if (openmmpolInput%inputFormat == "mmp") then
                call ommp_init_mmp(this%pSystem, openmmpolInput%paramsFilename)
            else
                call error("Bad openmmpol input format supplied to initializer!")
            end if

            ! TODO: add verbosity control
            call ommp_set_verbose(OMMP_VERBOSE_DEBUG)
            allocate(netCharges(nQMatoms))
            call ommp_init_qm_helper(this%pQMHelper, nQMatoms, qmAtomCoords, netCharges, atomTypes)
            deallocate(netCharges)
            #:else 
            call notImplementedError
            #:endif

        end subroutine

        subroutine TOMMPInterface_terminate(this)
            type(TOMMPInterface), intent(out) :: this
            #:if WITH_OPENMMPOL
            call ommp_terminate(this%pSystem)
            call ommp_terminate_qm_helper(this%pQMHelper)
            !> TODO: why do those deallocations cause segfault?
            ! deallocate(this%qmAtomsPotential)
            ! deallocate(this%qmmmCouplingEnergyPerAtom)
            #:else 
            call notImplementedError
            #:endif
        end subroutine

        subroutine updateQMCoords(this)
            class(TOMMPInterface), intent(inout) :: this

            #:if WITH_OPENMMPOL
            !> TODO: to be cleared
            this%pSystem%eel%ipd_done = .false.
            this%pSystem%eel%D2mgg_done = .false.
            this%pSystem%eel%D2dgg_done = .false.
            this%pQMHelper%E_n2p_done = .false.

            #:else 
            call notImplementedError
            #:endif

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

            #:if WITH_OPENMMPOL
            !> Debug: print a message every time charges are updated
            write(*, *) "Call to QM charges"

            allocate(qPerAtom(size(species)))

            !> getSummedCharges computes population, necessary to
            !  multiply by -1 to get correct charge sign
            call getSummedCharges(species, orb, qq, q0=q0, dQatom=qPerAtom)
            qPerAtom = -qPerAtom

            !> Set charges in the helper object
            this%pQMHelper%qqm = qPerAtom 
            deallocate(qPerAtom)

            !> Since charges are now updated, re-evaluation of
            !  charge-related quantities is requested
            this%pSystem%eel%ipd_done = .false.
            this%pSystem%eel%D2mgg_done = .false.
            this%pSystem%eel%D2dgg_done = .false.
            this%pQMHelper%E_n2p_done = .false.

            ! > Compute electric field produced by QM part of the system on MM atoms
            call ommp_prepare_qm_ele_ene(this%pSystem, this%pQMHelper)

            !> Set external field for MM, solve the polarization equations
            call ommp_set_external_field(this%pSystem, this%pQMHelper%E_n2p, this%solver, .true.)

            !> Store external potential for later access
            this%qmAtomsPotential(:) = this%pQMHelper%V_m2n
            this%qmmmCouplingEnergyPerAtom(:) = this%pQMHelper%qqm * this%qmAtomsPotential 
            
            !> Debug: prints total QM/MM coupling energy on every step
            ! write(*, "(A,F12.6, A)") "E_QMMM: ", sum(this%qmmmCouplingEnergyPerAtom) * 627.5, " kJ/mol"

            #:else 
            call notImplementedError
            #:endif
        end subroutine

        subroutine addPotential(this, shiftPerAtom)
            class(TOMMPInterface) :: this

            !> Computed potential is added to this vector
            real(dp), intent(inout) :: shiftPerAtom(:)
            
            !> Computed shift
            real(dp), allocatable :: shiftFromOpenmmpol(:)

            #:if WITH_OPENMMPOL
            !> Debug: print charges every time the potential is computed
            write(*, *) "Call to addPotential"
            write(*, *) this%pQMHelper%qqm

            !> Add potential from openmmpol to the vector of potentials
            shiftPerAtom = shiftPerAtom + this%pQMHelper%V_m2n

            #:else 
            call notImplementedError
            #:endif

        end subroutine

        subroutine addAtomEnergies(this, energyArray)
            class(TOMMPInterface) :: this

            !> Array to write energy into
            real(dp), intent(inout), allocatable :: energyArray

            #:if WITH_OPENMMPOL
            !> Write energy per atom (charge Q * potential V)
            ! @:ASSERT (size(energyArray) == size(this%qmAtomsPotential))
            ! energyArray = energyArray + this%pQMHelper%qqm * this%pQMHelper%V_m2n
            #:else 
            call notImplementedError
            #:endif

        end subroutine

        subroutine notImplementedError
            call error("DFTB+ compiled without support for openmmpol library")
        end subroutine notImplementedError
end module

