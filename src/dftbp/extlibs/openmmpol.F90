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
    !> debug
    use dftbp_type_integral, only : TIntegral
    use dftbp_dftb_sparse2dense, only : unpackHS, blockSymmetrizeHS
    use dftbp_type_multipole, only : TMultipole
    use dftbp_dftb_populations, only : getChargePerShell, denseSubtractDensityOfAtoms, mulliken,&
    & denseMulliken, denseBlockMulliken, skewMulliken, getOnsitePopulation, &
    & getAtomicMultipolePopulation
    !> end debug
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
        procedure :: testNumericalMatrixElementsDebug
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
      
            !>> Orbital information
            type(TOrbitals), intent(in) :: orb

            !> Charge per atom (internal variable)
            real(dp), allocatable :: qPerAtom(:)

            !> Sum of charges supplied to the solver (internal variable)
            real(dp) :: qSum

            #:if WITH_OPENMMPOL

            allocate(qPerAtom(size(species)))

            !> getSummedCharges computes population, necessary to
            !  multiply by -1 to get correct charge sign
            call getSummedCharges(species, orb, qq, q0=q0, dQatom=qPerAtom)
            qPerAtom = -qPerAtom
            qSum = sum(qPerAtom)

            ! write(*, *) "Charges in calculation:"
            ! write(*, *) qPerAtom
            !> Set charges in the helper object
            this%pQMHelper%qqm = qPerAtom 
            deallocate(qPerAtom)

            !> Since charges are now updated, re-evaluation of
            !  charge-related quantities is requested
            ! > Compute electric field produced by QM part of the system on MM atoms
            this%pQMHelper%E_n2p_done = .false. ! Only computes E_n2p (and V_m2n/V_p2n at first call)
            call ommp_prepare_qm_ele_ene(this%pSystem, this%pQMHelper)
            
            !> Set external field for MM, solve the polarization equations
            this%pSystem%eel%ipd_done = .false.
            this%pSystem%eel%D2mgg_done = .false.
            this%pSystem%eel%D2dgg_done = .false.

            call ommp_set_external_field(this%pSystem, this%pQMHelper%E_n2p, this%solver, .true.)

            ! > Compute electostatic potential produced by MM+Pol part on QM nucleai
            this%pQMHelper%V_p2n_done = .false. ! Only computes V_p2n, after having updated the external field / IPDs
            call ommp_prepare_qm_ele_ene(this%pSystem, this%pQMHelper)

            !> Store external potential for later access
            this%qmAtomsPotential = this%pQMHelper%V_m2n + this%pQMHelper%V_p2n
            this%qmmmCouplingEnergyPerAtom = this%pQMHelper%qqm * this%qmAtomsPotential
            
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
            ! write(*, *) "Call to addPotential"
            ! write(*, *) this%pQMHelper%qqm

            !> Add potential from openmmpol to the vector of potentials
            shiftPerAtom = shiftPerAtom + this%qmAtomsPotential

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

        subroutine testNumericalMatrixElementsDebug(this, env, rhoPrim, ints, orb, species, q0, neighbourList, nNeighbourSK,&
                                                 & iAtomStart, iSparseStart, img2CentCell)
            class(TOMMPInterface) :: this

            ! ARGUMENTS:
            type(TEnvironment), intent(in) :: env

            !> Non-symmetrized density matrix in sparse form
            real(dp), intent(in), allocatable :: rhoPrim(:, :)

            !> Integral container
            type(TIntegral), intent(in) :: ints

            !> Atomic orbital information
            type(TOrbitals), intent(in) :: orb

            !> Reference atomic populations
            real(dp), intent(in) :: q0(:, :, :)

            !> Species of each atoms
            integer, intent(in) :: species(:)
        
            !> Atomic neighbours
            type(TNeighbourList), intent(in) :: neighbourList

            !> Misc
            integer, intent(in) :: nNeighbourSK(:)
            integer, intent(in) :: iAtomStart(:)
            integer, intent(in) :: iSparseStart(:, :)
            integer, intent(in) :: img2CentCell(:)
            ! integer, intent(in) :: sqrHamSize

            ! INTERNAL VARIABLES:
            !> Imaginary component for the density matrix (do not allocate!)
            real(dp), allocatable :: iRhoPrim(:, :)

            !> Unperturbed density matrix
            real(dp), allocatable :: dm_0(:, :)

            !> Perturbed density matrix
            real(dp), allocatable :: dm_perturbed(:, :)

            !> Analytical hamiltonian contributions
            real(dp), allocatable :: h_exact(:)

            !> Numerical hamiltonian contributions
            real(dp), allocatable :: h_numeric(:)

            !> Charges:
            !> perturbed orbital (work)
            real(dp), allocatable :: qOrbPerturbed(:, :, :)
            !> unperturbed orbital (work)
            real(dp), allocatable :: qOrb0(:, :, :)
            !> dual atomic (do not allocate)
            real(dp), allocatable :: qBlock(:, :, :, :)
            !> imaginary dual atomic (do not allocate)
            real(dp), allocatable :: qiBlock(:, :, :, :)
            !> net atom (do not allocate)
            real(dp), allocatable :: qNetAtom(:)
            !> unperturbed net atom (work)
            real(dp), allocatable :: qPerAtom0(:)
            !> perturbed net atom (work)
            real(dp), allocatable :: qPerAtomPerturbed(:)
            !> end of charges
            
            !> Perturbed polarization energy
            real(dp) :: E_pol_perturbed

            !> Polarization energy on the unperturbed density matrix
            real(dp) :: E_pol_0

            !> Density matrix perturbation step
            real(dp) :: P_step = 1e-05

            !> Iterator index for the sparse matrix
            integer :: i

            !> Allocate necessary arrays
            allocate(dm_0, mold=rhoPrim)
            allocate(dm_perturbed, mold=rhoPrim)
            allocate(h_exact(size(dm_0, dim=1)))
            allocate(h_numeric(size(dm_0, dim=1)))
            allocate(qOrbPerturbed(orb%mOrb, size(orb%nOrbAtom), 1))
            allocate(qOrb0, mold=qOrbPerturbed)
            allocate(qPerAtom0(size(species)))
            allocate(qPerAtomPerturbed, mold=qPerAtom0)

            !> Init arrays
            dm_0 = rhoPrim(:, :) + 0.0_dp
            dm_perturbed = 0.0_dp
            h_exact = 0.0_dp
            h_numeric = 0.0_dp
            qOrbPerturbed = 0.0_dp
            qOrb0 = 0.0_dp
            qPerAtom0 = 0.0_dp
            qPerAtomPerturbed = 0.0_dp
            
            !> Compute the unperturbed population from density matrix
            call getMullikenPopulation(env, dm_0, ints, orb, neighbourList, nNeighbourSK,&
                                     & img2CentCell, iSparseStart, qOrb0, iRhoPrim, qBlock, qiBlock, qNetAtom)

            !> Get QM/MM energy at unperturbed density
            call this%updateQMCharges(env, species, neighbourList, qOrb0, q0, img2CentCell, orb)
            write(*, *) this%qmmmCouplingEnergyPerAtom
            E_pol_0 = sum(this%qmmmCouplingEnergyPerAtom)

            !> Unpack and block symmetrize the overlap matrix
            ! call unpackHS(overlap_dense, overlap, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart, img2CentCell)
            ! call blockSymmetrizeHS(overlap_dense, iAtomStart)
            
            !> Computational procedure for analytical matrix elements:

            !> Computational procedure for numerical matrix elements:
            ! For each pair of basis functions \mu \nu, vary density matrix element with a step dStep in [1e-4, 1e-7],
            ! then evaluate charges, then optimize QM/MM coupling energy for each set of charges,
            ! then compute (E_pol_1 - E_pol_0)/dStep - this is the Fock matrix contribution

            do i = 1, size(dm_0, 1)

                !> Construct perturbed density matrix
                dm_perturbed = dm_0 + 0.0_dp
                dm_perturbed(i, 1) = dm_perturbed(i, 1) + P_step

                !> Evaluate perturbed charges from perturbed density matrix
                call getMullikenPopulation(env, dm_perturbed, ints, orb, neighbourList, nNeighbourSK,&
                                         & img2CentCell, iSparseStart, qOrbPerturbed, iRhoPrim, qBlock, qiBlock, qNetAtom)

                !> Get variational energy for perturbed charges
                call this%updateQMCharges(env, species, neighbourList, qOrbPerturbed, q0, img2CentCell, orb)
                !
                write(*, *) "QM/MM coupling energy per atom:"
                write(*, *) this%qmmmCouplingEnergyPerAtom
                !
                E_pol_perturbed = sum(this%qmmmCouplingEnergyPerAtom)

                !> Evaluate numerical Fock matrix element
                h_numeric(i) = (E_pol_perturbed - E_pol_0) / P_step

                qOrbPerturbed = 0.0_dp

            end do
            
            !> Print numeric hamiltonian
            ! write(*, *) "Numeric Fock matrix contribution:"
            ! write(*, *) h_numeric
            !> Set charges back to initial
            call this%updateQMCharges(env, species, neighbourList, qOrb0, q0, img2CentCell, orb)

            !> Deallocate all temporary arrays
            deallocate(dm_0)
            deallocate(dm_perturbed)
            deallocate(h_exact)
            deallocate(h_numeric)
            deallocate(qOrb0)
            deallocate(qOrbPerturbed)
            deallocate(qPerAtom0)
            deallocate(qPerAtomPerturbed)

            #:if WITH_OPENMMPOL
            #:else 
            call notImplementedError
            #:endif

        end subroutine

        subroutine notImplementedError
            call error("DFTB+ compiled without support for openmmpol library")
        end subroutine notImplementedError

    !> Calculate Mulliken population from sparse density matrix.
        subroutine getMullikenPopulation(env, rhoPrim, ints, orb, neighbourList, nNeighbourSK,&
            & img2CentCell, iSparseStart, qOrb, iRhoPrim, qBlock, qiBlock, qNetAtom, multipoles)
      
          !> Environment settings
          type(TEnvironment), intent(in) :: env
      
          !> sparse density matrix
          real(dp), intent(in) :: rhoPrim(:,:)
      
          !> Integral container
          type(TIntegral), intent(in) :: ints
      
          !> Atomic orbital information
          type(TOrbitals), intent(in) :: orb
      
          !> Atomic neighbours
          type(TNeighbourList), intent(in) :: neighbourList
      
          !> Number of neighbours for each atom within overlap distance
          integer, intent(in) :: nNeighbourSK(:)
      
          !> image to actual atom indexing
          integer, intent(in) :: img2CentCell(:)
      
          !> sparse matrix indexing array
          integer, intent(in) :: iSparseStart(:,:)
      
          !> orbital charges
          real(dp), intent(out) :: qOrb(:,:,:)
      
          !> imaginary part of density matrix
          real(dp), intent(in), allocatable :: iRhoPrim(:,:)
      
          !> Dual atomic charges
          real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)
      
          !> Imaginary part of dual atomic charges
          real(dp), intent(inout), allocatable :: qiBlock(:,:,:,:)
      
          !> Onsite Mulliken charges per atom
          real(dp), intent(inout), allocatable :: qNetAtom(:)
      
          !> Multipole moments
          type(TMultipole), intent(inout), optional :: multipoles
      
          integer :: iSpin
      
          qOrb(:,:,:) = 0.0_dp
          do iSpin = 1, size(rhoPrim, dim=2)
            call mulliken(env, qOrb(:,:,iSpin), ints%overlap, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
                & nNeighbourSK, img2CentCell, iSparseStart)
          end do
      
          if (allocated(qBlock)) then
            qBlock(:,:,:,:) = 0.0_dp
            do iSpin = 1, size(rhoPrim, dim=2)
              call mulliken(env, qBlock(:,:,:,iSpin), ints%overlap, rhoPrim(:,iSpin), orb,&
                  & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
            end do
          end if
      
          if (allocated(qiBlock)) then
            qiBlock(:,:,:,:) = 0.0_dp
            do iSpin = 1, size(iRhoPrim, dim=2)
              call skewMulliken(env, qiBlock(:,:,:,iSpin), ints%overlap, iRhoPrim(:,iSpin), orb,&
                  & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
            end do
          end if
      
          if (allocated(qNetAtom)) then
            call getOnsitePopulation(rhoPrim(:,1), orb, iSparseStart, qNetAtom)
          end if
      
          if (present(multipoles)) then
      
            if (allocated(multipoles%dipoleAtom)) then
              call getAtomicMultipolePopulation(multipoles%dipoleAtom, ints%dipoleBra, ints%dipoleKet, &
                  & rhoPrim, orb, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, &
                  & iSparseStart)
            end if
      
            if (allocated(multipoles%quadrupoleAtom)) then
              call getAtomicMultipolePopulation(multipoles%quadrupoleAtom, ints%quadrupoleBra,&
                  & ints%quadrupoleKet, rhoPrim, orb, neighbourList%iNeighbour, nNeighbourSK,&
                  & img2CentCell, iSparseStart)
            end if
      
          end if
      
        end subroutine getMullikenPopulation
end module

