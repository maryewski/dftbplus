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
      ommp_get_full_ele_energy, OMMP_VERBOSE_DEBUG, OMMP_VERBOSE_LOW, OMMP_VERBOSE_HIGH, &
      ommp_get_fixedelec_energy, ommp_get_polelec_energy, ommp_get_full_bnd_energy, &
      ommp_get_vdw_energy, ommp_qm_helper_init_vdw_prm, ommp_qm_helper_vdw_energy,&
      ommp_qm_helper_set_attype, ommp_qm_helper_vdw_energy, ommp_print_summary_to_file, &
      ommp_ignore_duplicated_angle_prm, OMMP_SOLVER_NONE, OMMP_FF_AMOEBA, OMMP_FF_WANG_AL,&
      OMMP_FF_WANG_DL

#:endif
   use dftbp_io_message, only : error, warning
   use dftbp_common_environment, only : TEnvironment
   use dftbp_dftb_periodic, only : TNeighbourList
   use dftbp_type_commontypes, only : TOrbitals
   use dftbp_dftb_charges, only : getSummedCharges
   use dftbp_common_accuracy, only : dp, mc
   !> debug
   ! use dftbp_type_integral, only : TIntegral
   ! use dftbp_dftb_sparse2dense, only : unpackHS
   ! use dftbp_type_multipole, only : TMultipole
   ! use dftbp_dftb_populations, only : mulliken
   ! use dftbp_dftb_shift, only : addShift
   ! use dftbp_type_densedescr, only : TDenseDescr
   !> debug end
   implicit none

   !> TODO: import from openmmpol
   integer, protected :: OMMP_FF_AMBER = 0

   type :: TOMMPInterface

#:if WITH_OPENMMPOL
      type(ommp_system), pointer :: pSystem
      type(ommp_qm_helper), pointer :: pQMHelper
#:endif

      integer :: solver
      real(dp), allocatable :: sitePotentialStatic(:)
      real(dp), allocatable :: sitePotentialPolarizationEnergy(:)
      real(dp), allocatable :: sitePotentialPolarizationFock(:)
      real(dp) :: bondedEnergy
      real(dp) :: nonBondedEnergy
      real(dp) :: electrostaticEnergy
   contains
      ! procedure :: updateQMCoords
      procedure :: updateQMCharges
      procedure :: addPotential
      procedure :: writeOutput
      ! procedure :: FockMatrixDebugTest
   end type

   type :: TOMMPInput
      !> May be "Tinker" or "mmp", determines which
      !  init routine is used
      character(:), allocatable :: inputFormat
      !> Path to the MM geometry file
      character(:), allocatable :: mmGeomFilename
      !> Path to a separate .prm or .key file with
      ! Tinker-style parameters (if present)
      character(:), allocatable :: mmParamsFilename
      !> Determines which solver is utilized by openmmpol
      !  to compute polarization dipoles
      integer :: solver
      !> Provides "atom types" for QM atoms
      !  (used to get vdW parameters for QM/MM interaction)
      integer, allocatable :: qmAtomTypes(:)
      !> Path to .prm or .key file that contains
      ! "atom type" specifications that are to be used
      ! for QM atoms
      ! (used to get vdW parameters for QM/MM interaction)
      character(:), allocatable :: qmParamsFilename

   end type

   public TOMMPInterface, TOMMPInterface_init

contains

   subroutine TOMMPInterface_init(this, openmmpolInput, atomNames, atomTypes, qmAtomCoords)
      ! TODO: redo %speciesNames (tblite has something)
      type(TOMMPInterface), intent(out) :: this
      type(TOMMPInput), intent(in) :: openmmpolInput
      character(len=*), intent(in) :: atomNames(:)
      integer, intent(in) :: atomTypes(:)
      real(dp), intent(in) :: qmAtomCoords(:, :)
      integer :: nQMatoms
      real(dp), allocatable, dimension(:) :: chargesDummy
      integer, allocatable :: atomNumbersVector(:)

#:if WITH_OPENMMPOL
      nQMatoms = size(atomTypes)

      this%solver = openmmpolInput%solver
      this%electrostaticEnergy = 0.0_dp
      this%bondedEnergy = 0.0_dp
      this%nonBondedEnergy = 0.0_dp
      allocate(this%sitePotentialStatic(nQMatoms))
      this%sitePotentialStatic = 0.0_dp
      allocate(this%sitePotentialPolarizationEnergy(nQMatoms))
      this%sitePotentialPolarizationEnergy = 0.0_dp
      allocate(this%sitePotentialPolarizationFock(nQMatoms))
      this%sitePotentialPolarizationFock = 0.0_dp

      ! Tinker compatibility by default
      call ommp_ignore_duplicated_angle_prm()
      
      select case (openmmpolInput%inputFormat)
       case ('Tinker')
         call ommp_init_xyz(this%pSystem, openmmpolInput%mmGeomFilename, openmmpolInput%mmParamsFilename)
       case ('mmp')
         call ommp_init_mmp(this%pSystem, openmmpolInput%mmParamsFilename)
       case default
         call error("Bad openmmpol input format supplied to initializer!")
      end select

      ! TODO: add verbosity control
      call ommp_set_verbose(OMMP_VERBOSE_DEBUG)
      ! call ommp_set_verbose(OMMP_VERBOSE_LOW)
      

      allocate(chargesDummy(nQMatoms))
      allocate(atomNumbersVector(nQMatoms))
      chargesDummy = 0.0_dp
      call getAtomNumbersVector(atomNumbersVector, atomNames, atomTypes)
      call ommp_init_qm_helper(this%pQMHelper, nQMatoms, qmAtomCoords, chargesDummy, atomNumbersVector)
      deallocate(chargesDummy)
      deallocate(atomNumbersVector)

      ! If using AMOEBA, request potentials from the P set of dipoles
      if(this%pSystem%amoeba) then
         this%pQMHelper%V_pp2n_req = .true.
      end if

      ! Init vdW part of the QM/MM interaction
      call ommp_qm_helper_set_attype(this%pQMHelper, openmmpolInput%qmAtomTypes)
      call ommp_qm_helper_init_vdw_prm(this%pQMHelper, openmmpolInput%qmParamsFilename)

      !> Evaluate bonded and vdW terms at initial geometry
      this%bondedEnergy = ommp_get_full_bnd_energy(this%pSystem)
      this%nonBondedEnergy = ommp_get_vdw_energy(this%pSystem) + &
      &ommp_qm_helper_vdw_energy(this%pQMHelper, this%pSystem)


#:else
      call notImplementedError
#:endif

   end subroutine

   subroutine TOMMPInterface_terminate(this)
      type(TOMMPInterface), intent(inout) :: this
#:if WITH_OPENMMPOL
      call ommp_terminate_qm_helper(this%pQMHelper)
      call ommp_terminate(this%pSystem)
      deallocate(this%sitePotentialStatic)
      deallocate(this%sitePotentialPolarizationEnergy)
      deallocate(this%sitePotentialPolarizationFock)
#:else
      call notImplementedError
#:endif
   end subroutine

   subroutine updateQMCoords(this)
      class(TOMMPInterface) :: this

#:if WITH_OPENMMPOL
      !> TODO: to be cleared
      this%pSystem%eel%ipd_done = .false.
      this%pSystem%eel%D2mgg_done = .false.
      this%pSystem%eel%D2dgg_done = .false.
      this%pQMHelper%E_n2p_done = .false.

      this%bondedEnergy = ommp_get_full_bnd_energy(this%pSystem)
      this%nonBondedEnergy = ommp_get_vdw_energy(this%pSystem) + &
      &ommp_qm_helper_vdw_energy(this%pQMHelper, this%pSystem)
#:else
      call notImplementedError
#:endif

   end subroutine

   subroutine updateQMCharges(this, env, species, neighList, qq, q0, img2CentCell, orb)
      class(TOMMPInterface):: this

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

      if (this%pSystem%ff_type .ne. OMMP_FF_AMBER) then
         !> Set external field for MM, solve the polarization equations
         this%pSystem%eel%D2Mgg_done = .false.
         this%pSystem%eel%D2Dgg_done = .false.

         call ommp_set_external_field(this%pSystem, &
            this%pQMHelper%E_n2p, &
            this%solver, &
            OMMP_SOLVER_NONE, &
            .true.)

         !> Compute electostatic potential produced by MM+Pol part on QM nuclei
         this%pQMHelper%V_p2n_done = .false.
         if(this%pSystem%amoeba) then
            this%pQMHelper%V_pp2n_done = .false.
         end if

         !> Only computes V_p2n, after having updated the external field/IPDs
         call ommp_prepare_qm_ele_ene(this%pSystem, this%pQMHelper)

         this%sitePotentialPolarizationEnergy = -this%pQMHelper%V_p2n
         !> If using AMOEBA, the potential is mean for two sets of dipoles
         if (this%pSystem%amoeba) then
            this%sitePotentialPolarizationFock = -(this%pQMHelper%V_pp2n + this%pQMHelper%V_p2n) * 0.5_dp
         else
            this%sitePotentialPolarizationFock = -this%pQMHelper%V_p2n
         end if
      end if

      this%sitePotentialStatic = -this%pQMHelper%V_m2n

      !> Store total QM/MM electrostatic energy
      this%electrostaticEnergy = ommp_get_polelec_energy(this%pSystem)   &
         + ommp_get_fixedelec_energy(this%pSystem) &
         + dot_product(this%pQMHelper%qqm,         &
      & -(this%sitePotentialStatic + 0.5_dp * this%sitePotentialPolarizationEnergy))

#:else

      call notImplementedError

#:endif
   end subroutine

   subroutine addPotential(this, shiftPerAtom, qq, q0, env, species,&
   & neighbourList, img2CentCell, orb)

      class(TOMMPInterface):: this

      !> Computed potential is added to this vector
      real(dp), intent(inout) :: shiftPerAtom(:)

      !> Orbital-resolved Mulliken population
      real(dp), intent(in) :: qq(:, :, :)

      !> Reference orbital-resolved Mulliken population
      real(dp), intent(in) :: q0(:, :, :)

      !> Environment class
      type(TEnvironment), intent(in) :: env

      !> Species vector
      integer, intent(in) :: species(:)

      !> Neighbour list class
      type(TNeighbourList), intent(in) :: neighbourList

      !> Image-to-central-cell mapping
      integer, intent(in) :: img2CentCell(:)

      !> Orbital information class
      type(TOrbitals), intent(in) :: orb

#:if WITH_OPENMMPOL
      !> Add potential from openmmpol to the vector of potentials
      shiftPerAtom = shiftPerAtom + this%sitePotentialStatic + this%sitePotentialPolarizationFock

#:else

      call notImplementedError

#:endif

   end subroutine

   ! Write output
   subroutine writeOutput(this)
      class(TOMMPInterface) :: this

#:if WITH_OPENMMPOL
      call ommp_print_summary_to_file(this%pSystem, "openmmpol.out")
#:else
      call notImplementedError
#:endif

   end subroutine

   ! subroutine FockMatrixDebugTest(this, env, rhoPrim, ints, orb, species, q0, neighbourList, nNeighbourSK,&
   !                                           & iSparseStart, img2CentCell, denseDesc)
   !     class(TOMMPInterface) :: this

   !     ! ARGUMENTS:
   !     !> Computational environment
   !     type(TEnvironment), intent(in) :: env

   !     !> Non-symmetrized density matrix in sparse form
   !     real(dp), intent(in) :: rhoPrim(:, :)

   !     !> Integral container
   !     type(TIntegral), intent(in) :: ints

   !     !> Atomic orbital information
   !     type(TOrbitals), intent(in) :: orb

   !     !> Reference atomic populations
   !     real(dp), intent(in) :: q0(:, :, :)

   !     !> Species for all atoms
   !     integer, intent(in) :: species(:)

   !     !> Atomic neighbours
   !     type(TNeighbourList), intent(in) :: neighbourList

   !     !> Misc
   !     integer, intent(in) :: nNeighbourSK(:)
   !     integer, intent(in) :: iSparseStart(:, :)
   !     integer, intent(in) :: img2CentCell(:)
   !     type(TDenseDescr), intent(in) :: denseDesc

   !     ! INTERNAL VARIABLES:

   !     !> Perturbed density matrix
   !     real(dp), allocatable :: dm_perturbed(:, :)

   !     !> Analytical hamiltonian contributions (dense)
   !     real(dp), allocatable :: h_exact(:, :, :)

   !     !> Analytical hamiltonian contributions (sparse)
   !     real(dp), allocatable :: h_exact_sparse(:, :)

   !     !> Numerical hamiltonian contributions (dense)
   !     real(dp), allocatable :: h_numeric(:, :, :)

   !     !> Numeric hamiltonian contributions (sparse)
   !     real(dp), allocatable :: h_numeric_sparse(:, :)

   !     !> Perturbed orbital charges (work)
   !     real(dp), allocatable :: qOrbPerturbed(:, :, :)

   !     !> Unperturbed orbital (work)
   !     real(dp), allocatable :: qOrb0(:, :, :)

   !     !> Potential for exact matrix element evaluation
   !     real(dp), allocatable :: potential(:, :)

   !     !> Perturbed polarization energy (work)
   !     real(dp) :: E_pol_perturbed

   !     !> Polarization energy on the unperturbed density matrix
   !     real(dp) :: E_pol_0

   !     !> Numerical differention step
   !     real(dp), parameter :: diff_step = 1e-7_dp

   !     !> Iterator indices for rank 2 martrices
   !     integer :: i, j

   !     !> Number of spins in the system
   !     integer :: nSpin

   !     !> Iterator index for spins
   !     integer :: iSpin

   !     !> Number of atoms in the system
   !     integer :: nAtoms

   !     nSpin = size(rhoPrim, dim=2)
   !     nAtoms = size(species)

   !     !> Allocate and initialize other matrices
   !     allocate(dm_perturbed, mold=rhoPrim)
   !     dm_perturbed = 0.0_dp
   !     allocate(h_exact(denseDesc%fullSize, denseDesc%fullSize, nSpin))
   !     allocate(h_exact_sparse, mold=rhoPrim)
   !     allocate(h_numeric, mold=h_exact)
   !     allocate(h_numeric_sparse, mold=rhoPrim)
   !     h_exact = 0.0_dp
   !     h_numeric = 0.0_dp
   !     h_numeric_sparse = 0.0_dp
   !     h_exact_sparse = 0.0_dp
   !     allocate(potential(nAtoms, nSpin))
   !     potential = 0.0_dp

   !     !> Allocate and initialize charges and potentials
   !     allocate(qOrbPerturbed(orb%mOrb, size(orb%nOrbAtom), nSpin))
   !     allocate(qOrb0, mold=qOrbPerturbed)
   !     qOrb0 = 0.0_dp
   !     qOrbPerturbed = 0.0_dp

   !     !> Compute the unperturbed population from density matrix
   !     do iSpin = 1, nSpin
   !         call mulliken(env, qOrb0(:, :, iSpin), ints%overlap, rhoPrim(:, iSpin), &
   !                     & orb, neighbourList%iNeighbour, nNeighbourSK,&
   !                     & img2CentCell, iSparseStart)
   !     end do

   !     !> Compute QM/MM coupling energy at unpertubed density matrix
   !     E_pol_0 = this%electrostaticEnergy

   !     !> Get analytical potential
   !     call this%addPotential(potential(:, 1), qOrb0, q0, env, species,&
   !                                    & neighbourList, img2CentCell, orb)

   !     !> Compute analytic Fock matrix element
   !     call addShift(env, h_exact_sparse, ints%overlap, neighbourList%nNeighbour, neighbourList%iNeighbour,&
   !                       & species, orb, iSparseStart, nAtoms, img2CentCell, potential, .true.)

   !     !> Build and symmetrize analytic Fock matrix in dense form
   !     do iSpin = 1, nSpin
   !         call unpackHS(h_exact(:, :, iSpin), h_exact_sparse(:, iSpin), neighbourList%iNeighbour,&
   !                     & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell)
   !         call blockSymmetrizeHS(h_exact(:, :, iSpin), denseDesc%iAtomStart)
   !     end do

   !     !> Compute Fock matrix QM/MM terms by numeric differentiation in a loop
   !     do iSpin = 1, nSpin
   !         do i = 1, size(dm_perturbed, dim=1)
   !             qOrbPerturbed = 0.0_dp
   !             dm_perturbed = rhoPrim
   !             dm_perturbed(i, iSpin) = dm_perturbed(i, iSpin) + diff_step
   !             call mulliken(env, qOrbPerturbed(:, :, iSpin), ints%overlap, dm_perturbed(:, iSpin), &
   !                         & orb, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, &
   !                         & iSparseStart)
   !             call this%updateQMCharges(env, species, neighbourList, qOrbPerturbed, q0,&
   !                                     & img2CentCell, orb)
   !             E_pol_perturbed = this%electrostaticEnergy
   !             h_numeric_sparse(i, iSpin) = (E_pol_perturbed - E_pol_0) / diff_step
   !         end do
   !     end do

   !     !> Unpack sparse numerical Hamiltonian and correct off-diagonal values
   !     !> by multiplying them by 1/2
   !     do iSpin = 1, nSpin
   !         call unpackHS(h_numeric(:, :, iSpin), h_numeric_sparse(:, iSpin), &
   !                     & neighbourList%iNeighbour, nNeighbourSK, &
   !                     & denseDesc%iAtomStart, iSparseStart, img2CentCell)
   !         call blockSymmetrizeHS(h_numeric(:, :, iSpin), denseDesc%iAtomStart)
   !         do j = 1, denseDesc%fullSize
   !             do i = 1, denseDesc%fullSize
   !                 if (i /= j) then
   !                     h_numeric(i, j, iSpin) = h_numeric(i, j, iSpin) * 0.5_dp
   !                 end if
   !             end do
   !         end do
   !     end do

   !     !> Set charges back to initial
   !     call this%updateQMCharges(env, species, neighbourList, qOrb0, q0, img2CentCell, orb)

   !     !> Write analytical dense Fock matrix (spin channel 1)
   !     call write_array_newfile(h_exact(:, :, 1), "qmmm_fock_analytic.txt")

   !     !> Write numerical dense Fock matrix (spin channel 1)
   !     call write_array_newfile(h_numeric(:, :, 1), "qmmm_fock_numeric.txt")

   !     !> Write difference matrix
   !     call write_array_newfile(abs(h_exact(:, :, 1) - h_numeric(:, :, 1)), "qmmm_fock_difference.txt")

   !     !> Set small matrix elements to zero
   !     ! h_numeric = h_numeric * merge(1.0_dp, 0.0_dp, abs(h_numeric) >= 1e-10)

   !     !> Deallocate all temporary arrays
   !     deallocate(dm_perturbed)
   !     deallocate(h_exact)
   !     deallocate(h_exact_sparse)
   !     deallocate(h_numeric)
   !     deallocate(h_numeric_sparse)
   !     deallocate(qOrb0)
   !     deallocate(qOrbPerturbed)
   !     deallocate(potential)

   !     #:if WITH_OPENMMPOL
   !     #:else
   !     call notImplementedError
   !     #:endif

   ! end subroutine

   subroutine notImplementedError
      call error("DFTB+ compiled without support for openmmpol library")
   end subroutine notImplementedError

   ! subroutine denseBlockMullikenFullMatrix(rhoSqr, overSqr, iSquare, qq)
   !     !> Square spin polarized density matrix in dense format
   !     real(dp), intent(in) :: rhoSqr(:, :, :)

   !     !> Square overlap matrix in dense format
   !     real(dp), intent(in) :: overSqr(:, :)

   !     !> Atom positions in the row/column of square matrices
   !     integer, intent(in) :: iSquare(:)

   !     !> Mulliken charges on output (mOrb, nAtom, nSpin)
   !     real(dp), intent(out) :: qq(:,:,:)

   !     !> Element-wise product of P and S
   !     real(dp), allocatable :: A(:, :, :)

   !     integer :: iSpin, iAtom, nAtom, nSpin, nOrb, ii, jj

   !     qq(:, :, :) = 0.0_dp

   !     nAtom = size(qq, dim=2)
   !     nSpin = size(qq, dim=3)

   !     allocate(A, mold=rhoSqr)
   !     A = 0.0_dp

   !     !> Loop over spins
   !     do iSpin = 1, nSpin
   !         A(:, :, iSpin) = rhoSqr(:, :, iSpin) * overSqr
   !         !> Loop over atoms in cell
   !         do iAtom = 1, nAtom
   !             ii = iSquare(iAtom)
   !             jj = iSquare(iAtom+1)
   !             nOrb = jj - ii
   !             qq(:nOrb, iAtom, iSpin) = sum(A(ii:jj-1, :, iSpin), dim=2)
   !         end do
   !     end do

   !     deallocate(A)

   ! end subroutine denseBlockMullikenFullMatrix

   ! subroutine write_array_newfile(array, name)
   !     real(dp), intent(in) :: array(:, :)
   !     character(len=*), intent(in) :: name

   !     !> IO identifier
   !     integer :: io

   !     !> Array iterator
   !     integer :: i

   !     !> Open file for write
   !     open(newunit=io, file=name, status="replace", action="write")

   !     do i = 1, size(array, dim=1)
   !         write(io, *) array(i, :), ''
   !     end do

   !     !> Close file
   !     close(io)

   ! end subroutine write_array_newfile

   ! Given a species vector and a list of atom names,
   ! construct a vector of nuclear charges
   subroutine getAtomNumbersVector(numberVector, speciesNames, species)
      !> Nuclear charge vector
      integer, intent(inout) :: numberVector(:)
      !> Atom type names
      character(len=*), intent(in) :: speciesNames(:)
      !> Species in internal format
      integer, intent(in) :: species(:)

      !> Iterator index
      integer :: i
      !> Current species ordinal number
      integer :: currentSpecies

      do i = 1, size(species)
         currentSpecies = species(i)
         select case (speciesNames(currentSpecies))
          case ('H')
            numberVector(i) = 1
          case ('He')
            numberVector(i) = 2
          case ('Li')
            numberVector(i) = 3
          case ('Be')
            numberVector(i) = 4
          case ('B')
            numberVector(i) = 5
          case ('C')
            numberVector(i) = 6
          case ('N')
            numberVector(i) = 7
          case ('O')
            numberVector(i) = 8
          case ('F')
            numberVector(i) = 9
          case ('Ne')
            numberVector(i) = 10
          case ('Na')
            numberVector(i) = 11
          case ('Mg')
            numberVector(i) = 12
          case ('Al')
            numberVector(i) = 13
          case ('Si')
            numberVector(i) = 14
          case ('P')
            numberVector(i) = 15
          case ('S')
            numberVector(i) = 16
          case ('Cl')
            numberVector(i) = 17
          case ('Ar')
            numberVector(i) = 18
          case ('K')
            numberVector(i) = 19
          case ('Ca')
            numberVector(i) = 20
          case ('Sc')
            numberVector(i) = 21
          case ('Ti')
            numberVector(i) = 22
          case ('V')
            numberVector(i) = 23
          case ('Cr')
            numberVector(i) = 24
          case ('Mn')
            numberVector(i) = 25
          case ('Fe')
            numberVector(i) = 26
          case ('Co')
            numberVector(i) = 27
          case ('Ni')
            numberVector(i) = 28
          case ('Cu')
            numberVector(i) = 29
          case ('Zn')
            numberVector(i) = 30
          case ('Ga')
            numberVector(i) = 31
          case ('Ge')
            numberVector(i) = 32
          case ('As')
            numberVector(i) = 33
          case ('Se')
            numberVector(i) = 34
          case ('Br')
            numberVector(i) = 35
          case ('Kr')
            numberVector(i) = 36
          case ('Rb')
            numberVector(i) = 37
          case ('Sr')
            numberVector(i) = 38
          case ('Y')
            numberVector(i) = 39
          case ('Zr')
            numberVector(i) = 40
          case ('Nb')
            numberVector(i) = 41
          case ('Mo')
            numberVector(i) = 42
          case ('Tc')
            numberVector(i) = 43
          case ('Ru')
            numberVector(i) = 44
          case ('Rh')
            numberVector(i) = 45
          case ('Pd')
            numberVector(i) = 46
          case ('Ag')
            numberVector(i) = 47
          case ('Cd')
            numberVector(i) = 48
          case ('In')
            numberVector(i) = 49
          case ('Sn')
            numberVector(i) = 50
          case ('Sb')
            numberVector(i) = 51
          case ('Te')
            numberVector(i) = 52
          case ('I')
            numberVector(i) = 53
          case ('Xe')
            numberVector(i) = 54
          case ('Cs')
            numberVector(i) = 55
          case ('Ba')
            numberVector(i) = 56
          case ('La')
            numberVector(i) = 57
          case ('Ce')
            numberVector(i) = 58
          case ('Pr')
            numberVector(i) = 59
          case ('Nd')
            numberVector(i) = 60
          case ('Pm')
            numberVector(i) = 61
          case ('Sm')
            numberVector(i) = 62
          case ('Eu')
            numberVector(i) = 63
          case ('Gd')
            numberVector(i) = 64
          case ('Tb')
            numberVector(i) = 65
          case ('Dy')
            numberVector(i) = 66
          case ('Ho')
            numberVector(i) = 67
          case ('Er')
            numberVector(i) = 68
          case ('Tm')
            numberVector(i) = 69
          case ('Yb')
            numberVector(i) = 70
          case ('Lu')
            numberVector(i) = 71
          case ('Hf')
            numberVector(i) = 72
          case ('Ta')
            numberVector(i) = 73
          case ('W')
            numberVector(i) = 74
          case ('Re')
            numberVector(i) = 75
          case ('Os')
            numberVector(i) = 76
          case ('Ir')
            numberVector(i) = 77
          case ('Pt')
            numberVector(i) = 78
          case ('Au')
            numberVector(i) = 79
          case ('Hg')
            numberVector(i) = 80
          case ('Tl')
            numberVector(i) = 81
          case ('Pb')
            numberVector(i) = 82
          case ('Bi')
            numberVector(i) = 83
          case ('Po')
            numberVector(i) = 84
          case ('At')
            numberVector(i) = 85
          case ('Rn')
            numberVector(i) = 86
          case ('Fr')
            numberVector(i) = 87
          case ('Ra')
            numberVector(i) = 88
          case default
            call error("Unrecognized atom name")
         end select
      end do

   end subroutine getAtomNumbersVector

end module

