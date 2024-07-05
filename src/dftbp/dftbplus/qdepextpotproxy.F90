!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a proxy communicating with external generators of population dependent potentials.
module dftbp_dftbplus_qdepextpotproxy
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_shift, only : totalShift
  use dftbp_dftbplus_qdepextpotgen, only : TQDepExtPotGen, TQDepExtPotGenLLNode
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: TQDepExtPotProxy, TQDepExtPotProxy_init


  !> Collection of external q-dependent potentials queried during the SCC-cycle
  type :: TQDepExtPotProxy
    private
    !> collection of external potentials
    type(TQDepExtPotGenLLNode), pointer :: generatorsEntryNode
    !> energy contributions to DFTB atoms due to potentials
    real(dp), allocatable :: energyAtom(:)
    !> Total internal energy from all providers
    real(dp) :: energyInternal
  contains
    !> add potential contribution
    procedure :: addPotential => TQDepExtPotProxy_addPotential
    !> add energy contribution
    procedure :: addEnergy => TQDepExtPotProxy_addEnergy
    !> add internal energy contribtuion
    procedure :: addEnergyInternal => TQDepExtPotProxy_addEnergyInternal
    !> add force contribution
    procedure :: addGradientDc => TQDepExtPotProxy_addGradientDc
  end type TQDepExtPotProxy


contains

  !> Initializes proxy for querying a collection of q-dependent external potential generators.
  subroutine TQDepExtPotProxy_init(this, entryNode)

    !> Instance.
    type(TQDepExtPotProxy), intent(out) :: this

    !> External potential generators to consider.
    type(TQDepExtPotGenLLNode), intent(in), pointer :: entryNode

    this%generatorsEntryNode => entryNode

  end subroutine TQDepExtPotProxy_init


  !> Adds the external potential from the external potential generators.
  subroutine TQDepExtPotProxy_addPotential(this, deltaQAtom, deltaQShell, orb, species, potential)

    !> Instance.
    class(TQDepExtPotProxy), intent(inout) :: this

    !> Spin unpolarised net population per atom.
    real(dp), intent(in) :: deltaQAtom(:)

    !> Spin unpolarised net population per shell.
    real(dp), intent(in) :: deltaQShell(:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species
    integer, intent(in) :: species(:)

    !> Shell resolved potential to update.
    real(dp), intent(inout) :: potential(:,:,:,:)

    !> Current potential generator
    type(TQDepExtPotGenLLNode), pointer :: qmmmProviderCurrent

    real(dp), allocatable :: potAtom(:,:), potShell(:,:,:), potAtomTmp(:), potShellTmp(:,:)
    real(dp) :: energyInternalTmp
    integer :: mShell, nAtom, nSpin

    qmmmProviderCurrent => this%generatorsEntryNode  

    mShell = size(deltaQShell, dim=1)
    nAtom = size(potential, dim=3)
    nSpin = size(potential, dim=4)
    allocate(potAtom(nAtom, nSpin))
    allocate(potShell(mShell, nAtom, nSpin))
    allocate(potAtomTmp(nAtom))
    allocate(potShellTmp(mShell, nAtom))
    potAtom(:,:) = 0.0_dp
    potShell(:,:,:) = 0.0_dp
    if (.not. allocated(this%energyAtom)) then
      allocate(this%energyAtom(nAtom))
    end if
    this%energyAtom(:) = 0.0_dp
    this%energyInternal = 0.0_dp

    do while (associated(qmmmProviderCurrent))
      call qmmmProviderCurrent%instance%getExternalPot(deltaQAtom, deltaQShell, potAtomTmp,&
          & potShellTmp)
      
      ! Compute QM/MM interaction energy contribution
      this%energyAtom(:) = this%energyAtom + deltaQAtom * potAtomTmp
      this%energyAtom(:) = this%energyAtom + sum(deltaQShell * potShellTmp, dim=1)

      ! If Fock potential is screened, replace potential with the screened one
      ! before adding the potential that will be later used to compute the Fock
      ! matrix elements.
      if (qmmmProviderCurrent%instance%hasScreenedFock) then
        call qmmmProviderCurrent%instance%getExternalPotFock(deltaQAtom, deltaQShell, potAtomTmp,&
            & potShellTmp)
      end if

      potAtom(:,1) = potAtom(:,1) + potAtomTmp
      potShell(:,:,1) = potShell(:,:,1) + potShellTmp

      ! Add internal energy contributions if present
      if (qmmmProviderCurrent%instance%hasInternalEnergy) then
        call qmmmProviderCurrent%instance%getInternalEnergy(energyInternalTmp)
        this%energyInternal = this%energyInternal + energyInternalTmp
      end if

      qmmmProviderCurrent => qmmmProviderCurrent%next

    end do
    call totalShift(potShell, potAtom, orb, species)
    call totalShift(potential, potShell, orb, species)

  end subroutine TQDepExtPotProxy_addPotential


  !> Adds the energy contribution of the q-dependent external potentials.
  !>
  !> Note: This should only be called after the potential had been queried via the
  !> addPotential() procedure.
  !>
  subroutine TQDepExtPotProxy_addEnergy(this, energies)

    !> Instance.
    class(TQDepExtPotProxy), intent(inout) :: this

    !> Energy per atoms, to which the energies from the external potentials should be added to.
    real(dp), intent(inout) :: energies(:)

    energies(:) = energies + this%energyAtom

  end subroutine TQDepExtPotProxy_addEnergy


  !> Adds the internal MM contribution to the total QM/MM energy.
  !!
  !! Note: should only be called after calling addPotential().
  subroutine TQDepExtPotProxy_addEnergyInternal(this, energy)

    !> Instance.
    class(TQDepExtPotProxy), intent(inout) :: this

    !> Energy to add to
    real(dp), intent(inout) :: energy

    energy = energy + this%energyInternal

  end subroutine TQDepExtPotProxy_addEnergyInternal


  !> Adds the "double counting" gradient contribution of the q-dependent external potentials.
  !>
  !> The double counting part of the gradients is the one, which is not obtained by the derivatives
  !> of the shift-vectors.
  !>
  subroutine TQDepExtPotProxy_addGradientDc(this, deltaQAtom, deltaQShell, gradients)

    !> Instance.
    class(TQDepExtPotProxy), intent(inout) :: this

    !> Net population per atom.
    real(dp), intent(in) :: deltaQAtom(:)

    !> Net population per shell.
    real(dp), intent(in) :: deltaQShell(:,:)

    !> Gradients to upgrade.
    real(dp), intent(inout) :: gradients(:,:)

    real(dp), allocatable :: extPotGrad(:,:), deltaQAtomSpread(:,:)

    type(TQDepExtPotGenLLNode), pointer :: qmmmProviderCurrent
    qmmmProviderCurrent => this%generatorsEntryNode

    allocate(extPotGrad(size(gradients, dim=1), size(gradients, dim=2)))
    deltaQAtomSpread = spread(deltaQAtom, 1, 3)

    do while (associated(qmmmProviderCurrent))
      call qmmmProviderCurrent%instance%getExternalPotGrad(deltaQAtom, deltaQShell, extPotGrad)
      gradients(:,:) = gradients + deltaQAtomSpread * extPotGrad

      qmmmProviderCurrent => qmmmProviderCurrent%next
    end do

  end subroutine TQDepExtPotProxy_addGradientDc


end module dftbp_dftbplus_qdepextpotproxy
