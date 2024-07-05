!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains the interface for an external population dependent potential.
module dftbp_dftbplus_qdepextpotgen
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TQDepExtPotGen, TQDepExtPotGenLLNode

  !> Base class for generating external population dependent potentials.
  !>
  !> It should be extended whenever DFTB+ needs to be coupled with an external potential provider.
  !> Additional attributes and methods can be added in order to ease the conversion between the
  !> external tool and DFTB+.
  !>
  type, abstract ::   TQDepExtPotGen
    !> specifies whether the generator has screened
    !! Fock potential
    logical :: hasScreenedFock = .false.
    !> specifies whether the generator has internal
    !! energy that should be added to the total energy
    logical :: hasInternalEnergy = .false.
  contains
    !> get the external potential
    procedure(getExtPotIface), deferred :: getExternalPot
    !> get the external screened potential for the Fock matrix
    procedure(getExtPotFockIface), deferred :: getExternalPotFock
    !> get the gradient of the potential wrt DFTB atom positions
    procedure(getExtPotGradIface), deferred :: getExternalPotGrad
    !> get internal energy of the model
    procedure(getInternalEnergyIface), deferred :: getInternalEnergy
  end type TQDepExtPotGen


  !>  Linked list wrapper around external potential generator instance
  type :: TQDepExtPotGenLLNode

    !> TQDepExtPotGen instance to wrap.
    class(TQDepExtPotGen), allocatable :: instance

    !> Index of LL node
    integer :: index

    !> Pointer to the next LL node
    class(TQDepExtPotGenLLNode), pointer :: next => null()

  end type TQDepExtPotGenLLNode


  abstract interface

    !> Called, when DFTB+ needs the value of the external potential at the position of the atoms.
    !>
    !> The routine is called in every SCC iteration so that the external potential can be updated
    !> according to the atom charges. The actual implementation should return both potential
    !> contributions (or zero them out, if not needed). The atom and shell-resolved contribution
    !> will be then added in DFTB+.
    !>
    !> Note: External potential is defined as the external potential the electrons feel, so in case
    !> of an electrostatic potential you would have to invert its sign.
    !>
    subroutine getExtPotIface(this, chargePerAtom, chargePerShell, extPotAtom, extPotShell)
      import :: TQDepExtPotGen, dp

      !> Instance.
      class(TQDepExtPotGen), intent(inout) :: this

      !> Net number of electrons on the atom with respect to its reference configuration, the
      !> neutral atom.  Shape: [nAtom].
      real(dp), intent(in) :: chargePerAtom(:)

      !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
      real(dp), intent(in) :: chargePerShell(:,:)

      !> Atom dependent external potential contribution. Shape: [nAtom]
      real(dp), intent(out) :: extPotAtom(:)

      !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
      real(dp), intent(out) :: extPotShell(:,:)

    end subroutine getExtPotIface


    !> Called in cases where the potential that enters the Fock matrix
    !! differs from the electrostatic potential on the partial charges.
    !! This happens, for example, in AMOEBA QM/MM scenario, where
    !! the dipole-dipole interaction between polarizable MM atoms
    !! is screened, therefore a new contribution to the Fock matrix
    !! appears.
    subroutine getExtPotFockIface(this, chargePerAtom, chargePerShell, extPotAtom, extPotShell)
      import :: TQDepExtPotGen, dp

      !> Instance.
      class(TQDepExtPotGen), intent(inout) :: this

      !> Net number of electrons on the atom with respect to its reference configuration, the
      !> neutral atom.  Shape: [nAtom].
      real(dp), intent(in) :: chargePerAtom(:)

      !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
      real(dp), intent(in) :: chargePerShell(:,:)

      !> Atom dependent external potential contribution. Shape: [nAtom]
      real(dp), intent(out) :: extPotAtom(:)

      !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
      real(dp), intent(out) :: extPotShell(:,:)

    end subroutine getExtPotFockIface


    !> Called when DFTB needs the gradient of the external potential at the position of the atoms.
    !>
    !>
    !> The routine is only called once after finishing the SCC iteration, provides the calculator has
    !> been set up to calculate forces.
    !> Note: External potential is defined as the external potential the electrons feel, so in case
    !> of an electrostatic potential you would have to invert its sign.
    !>
    subroutine getExtPotGradIface(this, chargePerAtom, chargePerShell, extPotGrad)
      import :: TQDepExtPotGen, dp

      !> Class instance.
      class(TQDepExtPotGen), intent(inout) :: this

      !> Net number of electrons on the atom with respect to its reference configuration, the
      !> neutral atom.  Shape: [nAtom].
      real(dp), intent(in) :: chargePerAtom(:)

      !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
      real(dp), intent(in) :: chargePerShell(:,:)

      !> Gradient of the potential at each atom. Shape: [3, nAtom].
      real(dp), intent(out) :: extPotGrad(:,:)

    end subroutine getExtPotGradIface


    subroutine getInternalEnergyIface(this, energyInternal)
      import :: TQDepExtPotGen, dp

      !> Class instance.
      class(TQDepExtPotGen), intent(inout) :: this

      !> External energy contribution
      real(dp), intent(out) :: energyInternal
      
    end subroutine getInternalEnergyIface


  end interface


end module dftbp_dftbplus_qdepextpotgen
