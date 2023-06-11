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
    use dftbp_dftb_sparse2dense, only : unpackHS, blockSymmetrizeHS, symmetrizeHS
    use dftbp_type_multipole, only : TMultipole
    use dftbp_dftb_populations, only : getChargePerShell, denseSubtractDensityOfAtoms, mulliken,&
    & denseMulliken, denseBlockMulliken, skewMulliken, getOnsitePopulation, &
    & getAtomicMultipolePopulation
    use dftbp_dftb_shift, only : addShift
    use dftbp_type_densedescr, only : TDenseDescr
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
        procedure :: bigMatrixElementDebugTest
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
            this%qmAtomsPotential = 0.0_dp
            allocate(this%qmmmCouplingEnergyPerAtom(nQMatoms))
            this%qmmmCouplingEnergyPerAtom = 0.0_dp

            if (openmmpolInput%inputFormat == "Tinker") then
                call ommp_init_xyz(this%pSystem, openmmpolInput%geomFilename, openmmpolInput%paramsFilename)
            else if (openmmpolInput%inputFormat == "mmp") then
                call ommp_init_mmp(this%pSystem, openmmpolInput%paramsFilename)
            else
                call error("Bad openmmpol input format supplied to initializer!")
            end if

            ! TODO: add verbosity control
            ! call ommp_set_verbose(OMMP_VERBOSE_DEBUG)
            call ommp_set_verbose(OMMP_VERBOSE_LOW)
            allocate(netCharges(nQMatoms))
            netCharges = 0.0_dp
            call ommp_init_qm_helper(this%pQMHelper, nQMatoms, qmAtomCoords, netCharges, atomTypes)
            deallocate(netCharges)
            #:else 
            call notImplementedError
            #:endif

        end subroutine

        subroutine TOMMPInterface_terminate(this)
            type(TOMMPInterface), intent(inout) :: this
            #:if WITH_OPENMMPOL
            call ommp_terminate_qm_helper(this%pQMHelper)
            call ommp_terminate(this%pSystem)
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

            #:if WITH_OPENMMPOL

            !> DFTB+ internally uses charges and potential that are opposite to the chemical
            !  convention, therefore we need a sign inversion when passing them to openmmpol
            call getSummedCharges(species, orb, qq, q0=q0, dQatom=this%pQMHelper%qqm)
            this%pQMHelper%qqm = -this%pQMHelper%qqm

            !> Compute electric field produced by QM part of the system on MM atoms
            this%pQMHelper%E_n2p_done = .false. 

            !> Only computes E_n2p (and V_m2n/V_p2n at first call)
            call ommp_prepare_qm_ele_ene(this%pSystem, this%pQMHelper)
            
            !> Set external field for MM, solve the polarization equations
            this%pSystem%eel%ipd_done = .false.
            this%pSystem%eel%D2mgg_done = .false.
            this%pSystem%eel%D2dgg_done = .false.

            call ommp_set_external_field(this%pSystem, &
                                         this%pQMHelper%E_n2p, &
                                         this%solver, .true.)

            !> Compute electostatic potential produced by MM+Pol part on QM nuc
            this%pQMHelper%V_p2n_done = .false. 

            !> Only computes V_p2n, after having updated the external field/IPDs
            call ommp_prepare_qm_ele_ene(this%pSystem, this%pQMHelper)

            !> Store external potential for later access
            this%qmAtomsPotential = -(this%pQMHelper%V_m2n + 2*this%pQMHelper%V_p2n)
            this%qmmmCouplingEnergyPerAtom = this%pQMHelper%qqm * (this%pQMHelper%V_m2n + this%pQMHelper%V_p2n)
            
            write(*, *) "Constant potential:"
            write(*, *) this%pQMHelper%V_m2n
            write(*, *) "Polarization potential:"
            write(*, *) this%pQMHelper%V_p2n
            !> Debug: prints total QM/MM coupling energy on every step
            write(*, "(A,F12.6, A)") "E_QMMM: ", sum(this%qmmmCouplingEnergyPerAtom) * 627.5, " kJ/mol"

            #:else 
            call notImplementedError
            #:endif
        end subroutine

        subroutine addPotential(this, shiftPerAtom)
            class(TOMMPInterface) :: this

            !> Computed potential is added to this vector
            real(dp), intent(inout) :: shiftPerAtom(:)

            #:if WITH_OPENMMPOL

            !> Add potential from openmmpol to the vector of potentials
            shiftPerAtom = shiftPerAtom + this%qmAtomsPotential

            #:else 
            call notImplementedError
            #:endif

        end subroutine

        subroutine addAtomEnergies(this, energyArray)
            class(TOMMPInterface) :: this

            !> Array to write energy into
            real(dp), intent(inout), allocatable :: energyArray(:)

            #:if WITH_OPENMMPOL
            !> Write energy per atom (charge Q * potential V)
            ! @:ASSERT (size(energyArray) == size(this%qmAtomsPotential))
            ! energyArray = energyArray + this%pQMHelper%qqm * this%pQMHelper%V_m2n
            #:else 
            call notImplementedError
            #:endif

        end subroutine

        subroutine bigMatrixElementDebugTest(this, env, rhoPrim, ints, orb, species, q0, neighbourList, nNeighbourSK,&
                                                  & iSparseStart, img2CentCell, denseDesc)
            class(TOMMPInterface) :: this

            ! ARGUMENTS:
            !> Computational environment
            type(TEnvironment), intent(in) :: env

            !> Non-symmetrized density matrix in sparse form
            real(dp), intent(in) :: rhoPrim(:, :)

            !> Integral container
            type(TIntegral), intent(in) :: ints

            !> Atomic orbital information
            type(TOrbitals), intent(in) :: orb

            !> Reference atomic populations
            real(dp), intent(in) :: q0(:, :, :)

            !> Species for all atoms
            integer, intent(in) :: species(:)
        
            !> Atomic neighbours
            type(TNeighbourList), intent(in) :: neighbourList

            !> Misc
            integer, intent(in) :: nNeighbourSK(:)
            integer, intent(in) :: iSparseStart(:, :)
            integer, intent(in) :: img2CentCell(:)
            type(TDenseDescr), intent(in) :: denseDesc

            ! INTERNAL VARIABLES:
            !> Unperturbed density matrix
            real(dp), allocatable :: dm_0(:, :, :)

            !> Perturbed density matrix
            real(dp), allocatable :: dm_perturbed(:, :, :)

            !> Dense overlap matrix
            real(dp), allocatable :: overlap_dense(:, :)

            !> Analytical hamiltonian contributions (dense)
            real(dp), allocatable :: h_exact(:, :, :)

            !> Analytical hamiltonian contributions (sparse)
            real(dp), allocatable :: h_exact_sparse(:, :)

            !> Numerical hamiltonian contributions
            real(dp), allocatable :: h_numeric(:, :, :)

            !> Perturbed orbital charges (work)
            real(dp), allocatable :: qOrbPerturbed(:, :, :)

            !> Unperturbed orbital (work)
            real(dp), allocatable :: qOrb0(:, :, :)

            !> Potential for exact matrix element evaluation
            real(dp), allocatable :: potential(:, :)

            !> Perturbed polarization energy (work)
            real(dp) :: E_pol_perturbed

            !> Polarization energy on the unperturbed density matrix
            real(dp) :: E_pol_0

            !> Numerical differention step
            real(dp), parameter :: diff_step = 1e-7_dp

            !> Iterator indices for rank 2 martrices
            integer :: i, j

            !> Number of spins in the system
            integer :: nSpin

            !> Iterator index for spins
            integer :: iSpin

            !> Number of atoms in the system
            integer :: nAtoms

            !> Derivative factor
            real(dp) :: derFactor = 1.0_dp

            nSpin = size(rhoPrim, dim=2)
            nAtoms = size(species)

            !> Allocate and initialize unperturbed density matrix
            allocate(dm_0(denseDesc%fullSize, denseDesc%fullSize, nSpin))
            dm_0 = 0.0_dp
            do iSpin = 1, nSpin
                call unpackHS(dm_0(:, :, iSpin), rhoPrim(:, iSpin), neighbourList%iNeighbour, nNeighbourSK,&
                            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
                call blockSymmetrizeHS(dm_0(:, :, iSpin), denseDesc%iAtomStart)
            end do

            !> Allocate and initialize other matrices
            allocate(dm_perturbed, mold=dm_0)
            dm_perturbed = 0.0_dp
            allocate(h_exact, mold=dm_0)
            allocate(h_numeric, mold=dm_0)
            allocate(h_exact_sparse, mold=rhoPrim)
            h_exact = 0.0_dp
            h_numeric = 0.0_dp
            h_exact_sparse = 0.0_dp
            allocate(potential(nAtoms, 1))
            potential = 0.0_dp

            !> Allocate and initialize charges and potentials
            allocate(qOrbPerturbed(orb%mOrb, size(orb%nOrbAtom), 1))
            allocate(qOrb0, mold=qOrbPerturbed)
            qOrb0 = 0.0_dp
            qOrbPerturbed = 0.0_dp

            !> Make dense overlap matrix
            allocate(overlap_dense(denseDesc%fullSize, denseDesc%fullSize))
            overlap_dense = 0.0_dp
            call unpackHS(overlap_dense, ints%overlap, neighbourList%iNeighbour, nNeighbourSK,&
                        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
            call blockSymmetrizeHS(overlap_dense, denseDesc%iAtomStart)

            !> Debug: write overlap and unperturbed density matrices
            call write_array_newfile(overlap_dense, "overlap.txt")
            call write_array_newfile(dm_0(:, :, 1), "dm.txt")

            !> Compute the unperturbed population from density matrix
            call denseBlockMullikenFullMatrix(dm_0, overlap_dense, denseDesc%iAtomStart, qOrb0)

            ! !> Get QM/MM energy at unperturbed density
            ! call this%updateQMCharges(env, species, neighbourList, qOrb0, q0, img2CentCell, orb)

            !> Compute QM/MM coupling energy at unpertubed density matrix
            E_pol_0 = sum(this%qmmmCouplingEnergyPerAtom)

            !> Compute Fock matrix QM/MM terms by numeric differentiation in a loop
            do iSpin = 1, nSpin
                do j = 1, denseDesc%fullSize
                    do i = 1, denseDesc%fullSize
                        !> Construct perturbed density matrix
                        dm_perturbed = dm_0
                        dm_perturbed(i, j, iSpin) = dm_perturbed(i, j, iSpin) + diff_step
                        if (i /= j) then
                            dm_perturbed(j, i, iSpin) = dm_perturbed(j, i, iSpin) + diff_step 
                            derFactor = 0.5_dp
                        else
                            derFactor = 1.0_dp
                        end if
                        
                        !> Compute charges from perturbed density matrix
                        call denseBlockMullikenFullMatrix(dm_perturbed, overlap_dense, denseDesc%iAtomStart, qOrbPerturbed)
        
                        !> Get variational energy for perturbed charges
                        call this%updateQMCharges(env, species, neighbourList, qOrbPerturbed, q0, img2CentCell, orb)
                        E_pol_perturbed = sum(this%qmmmCouplingEnergyPerAtom)
        
                        !> Evaluate numerical Fock matrix element
                        h_numeric(i, j, iSpin) = derFactor * (E_pol_perturbed - E_pol_0) / diff_step

                    end do
                end do
            end do

            !> Set charges back to initial
            call this%updateQMCharges(env, species, neighbourList, qOrb0, q0, img2CentCell, orb)

            !> Get analytical potential
            call this%addPotential(potential(:, 1))
            ! call this%addKVectorNumeric(potential(:, 1), qOrb0, q0, env, species, &
                                    !   & neighbourList, img2CentCell, orb)

            !> Compute analytic Fock matrix element
            call addShift(env, h_exact_sparse, ints%overlap, neighbourList%nNeighbour, neighbourList%iNeighbour,&
                              & species, orb, iSparseStart, size(orb%nOrbAtom), img2CentCell, potential, .true.)

            !> Build and symmetrize analytic Fock matrix in dense form
            do iSpin = 1, nSpin
                call unpackHS(h_exact(:, :, iSpin), h_exact_sparse(:, iSpin), neighbourList%iNeighbour,&
                            & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell)
                call blockSymmetrizeHS(h_exact(:, :, iSpin), denseDesc%iAtomStart)
            end do

            !> Write analytical dense Fock matrix (spin channel 1)
            call write_array_newfile(h_exact(:, :, 1), "qmmm_fock_analytic.txt")

            !> Write numerical dense Fock matrix (spin channel 1)
            call write_array_newfile(h_numeric(:, :, 1), "qmmm_fock_numeric.txt")

            !> Total difference:
            ! write(*, *) "Total difference between analytic and numeric Hamiltonians:"
            ! write(*, *) sum(abs(h_exact(:, :, 1) - h_numeric(:, :, 1)))
            
            !> Write difference matrix
            call write_array_newfile(abs(h_exact(:, :, 1) - h_numeric(:, :, 1)), "qmmm_fock_difference.txt")
            
            !> Set small matrix elements to zero
            ! h_numeric = h_numeric * merge(1.0_dp, 0.0_dp, abs(h_numeric) >= 1e-10)

            
            !> Deallocate all temporary arrays
            deallocate(dm_0)
            deallocate(dm_perturbed)
            deallocate(h_exact)
            deallocate(h_exact_sparse)
            deallocate(h_numeric)
            deallocate(qOrb0)
            deallocate(qOrbPerturbed)
            deallocate(overlap_dense)
            deallocate(potential)

            #:if WITH_OPENMMPOL
            #:else 
            call notImplementedError
            #:endif

        end subroutine

        subroutine notImplementedError
            call error("DFTB+ compiled without support for openmmpol library")
        end subroutine notImplementedError

        subroutine denseBlockMullikenFullMatrix(rhoSqr, overSqr, iSquare, qq)
            !> Square spin polarized density matrix in dense format
            real(dp), intent(in) :: rhoSqr(:, :, :)

            !> Square overlap matrix in dense format
            real(dp), intent(in) :: overSqr(:, :)
        
            !> Atom positions in the row/column of square matrices
            integer, intent(in) :: iSquare(:)
        
            !> Mulliken charges on output (mOrb, nAtom, nSpin)
            real(dp), intent(out) :: qq(:,:,:)

            !> Element-wise product of P and S
            real(dp), allocatable :: A(:, :, :)

            integer :: iSpin, iAtom, nAtom, nSpin, nOrb, ii, jj

            qq(:, :, :) = 0.0_dp

            nAtom = size(qq, dim=2)
            nSpin = size(qq, dim=3)

            allocate(A, mold=rhoSqr)
            A = 0.0_dp

            !> Loop over spins
            do iSpin = 1, nSpin
                A(:, :, iSpin) = rhoSqr(:, :, iSpin) * overSqr
                !> Loop over atoms in cell
                do iAtom = 1, nAtom
                    ii = iSquare(iAtom)
                    jj = iSquare(iAtom+1)
                    nOrb = jj - ii
                    qq(:nOrb, iAtom, iSpin) = sum(A(ii:jj-1, :, iSpin), dim=2)
                end do
            end do

            deallocate(A)

        end subroutine denseBlockMullikenFullMatrix

        subroutine write_array_newfile(array, name)
            real(dp), intent(in) :: array(:, :)
            character(len=*), intent(in) :: name

            !> IO identifier
            integer :: io

            !> Array iterator
            integer :: i

            !> Open file for write
            open(newunit=io, file=name, status="replace", action="write")

            do i = 1, size(array, dim=1)
                write(io, *) array(i, :), ''
            end do

            !> Close file
            close(io)

        end subroutine write_array_newfile

end module

