!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Simplified C-interface with callbacks for population dependant external potential generators.
module dftbp_dftbplus_qdepextpotgenc
  use, intrinsic :: iso_c_binding
  use dftbp_common_accuracy, only : dp
  use dftbp_dftbplus_qdepextpotgen, only : TQDepExtPotGen
  implicit none
  private

  public :: getExtPotIfaceC, getExtPotGradIfaceC, getExtPotFockIfaceC, getInternalEnergyIfaceC
  public :: TQDepExtPotGenC, TQDepExtPotGenC_init

  !> Interface to the routine which calculates the external potential due to charges
  abstract interface

    !> Interface to set up external potential
    subroutine getExtPotIfaceC(refPtr, dQAtom, extPotAtom) bind(C)
      import :: c_double, c_ptr

      !> Reference pointer
      type(c_ptr), value, intent(in) :: refPtr

      !> Net number of electrons on each atom (note: positive number means electron excess)
      real(c_double), intent(in)  :: dQAtom(*)

      !> Potential on each atom (note: positive number means electron repulsion)
      real(c_double), intent(out) :: extPotAtom(*)

    end subroutine getExtPotIfaceC

    !> Interface to set up screened Fock external potential

    subroutine getExtPotFockIfaceC(refPtr, dQAtom, extPotAtom) bind(C)
      import :: c_double, c_ptr

      !> Reference pointer
      type(c_ptr), value, intent(in) :: refPtr

      !> Net number of electrons on each atom (note: positive number means electron excess)
      real(c_double), intent(in)  :: dQAtom(*)

      !> Potential on each atom (note: positive number means electron repulsion)
      real(c_double), intent(out) :: extPotAtom(*)

    end subroutine getExtPotFockIfaceC


    !> Interface to set up gradient of external potential
    subroutine getExtPotGradIfaceC(refPtr, dQAtom, extPotAtomGrad) bind(C)
      import :: c_double, c_ptr

      !> Reference pointer
      type(c_ptr), value, intent(in) :: refPtr

      !> Net number of electrons on each atom (note: positive number means electron excess)
      real(c_double), intent(in) :: dQAtom(*)

      !> Gradient of the potential on each atom (note: positive number means electron repulsion)
      real(c_double), intent(out) :: extPotAtomGrad(3, *)

    end subroutine getExtPotGradIfaceC

    !> Interface to set up external potential
    subroutine getInternalEnergyIfaceC(refPtr, energyInternal) bind(C)
      import :: c_double, c_ptr

      !> Reference pointer
      type(c_ptr), value, intent(in) :: refPtr

      !> Potential on each atom (note: positive number means electron repulsion)
      real(c_double), intent(out) :: energyInternal

    end subroutine getInternalEnergyIfaceC

  end interface

  !> Builds on the charge dependent external interface type from dftbp_dftbplus_qdepextpotgen for
  !> the C API
  type, extends(TQDepExtPotGen) :: TQDepExtPotGenC
    private
    type(c_ptr) :: refPtr 
    procedure(getExtPotIfaceC), nopass, pointer :: getExtPot
    procedure(getExtPotGradIfaceC), nopass, pointer :: getExtPotGrad
    procedure(getExtPotFockIfaceC), nopass, pointer :: getExtPotFock
    procedure(getInternalEnergyIfaceC), nopass, pointer :: getIntEnergy
  contains
    procedure :: getExternalPot => TDepExtPotGenC_getExternalPot
    procedure :: getExternalPotGrad => TQDepExtPotGenC_getExternalPotGrad
    procedure :: getExternalPotFock => TDepExtPotGenC_getExternalPotFock
    procedure :: getInternalEnergy => TDepExtPotGenC_getInternalEnergy
  end type TQDepExtPotGenC


contains


  !> Initialise an external charge dependent external potential within this type
  subroutine TQDepExtPotGenC_init(this, refPtr, extPotFunc, extPotGradFunc, extPotFockFunc,&
        & getInternalEnergyFunc)

    !> Instance
    type(TQDepExtPotGenC), intent(out) :: this

    !> Pointer to the C routine for the external potential
    type(c_ptr), intent(in) :: refPtr

    !> Function for the potential evaluation
    procedure(getExtPotIfaceC), pointer, intent(in) :: extPotFunc

    !> Function for the gradient of the potential
    procedure(getExtPotGradIfaceC), pointer, intent(in) :: extPotGradFunc

    !> Function for the gradient of the potential
    procedure(getExtPotFockIfaceC), pointer, intent(in) :: extPotFockFunc

    !> Function for getting internal energy
    procedure(getInternalEnergyIfaceC), pointer, intent(in) :: getInternalEnergyFunc

    this%getExtPot => extPotFunc
    this%getExtPotGrad => extPotGradFunc
    this%getExtPotFock => extPotFockFunc
    this%getIntEnergy => getInternalEnergyFunc
    this%refPtr = refPtr

  end subroutine TQDepExtPotGenC_init



  !> Extra routine a charge dependent external potential
  subroutine TDepExtPotGenC_getExternalPot(this, chargePerAtom, chargePerShell, extPotAtom,&
      & extPotShell)

    !> Instance.
    class(TQDepExtPotGenC), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Atom dependent external potential contribution. Shape: [nAtom]
    real(dp), intent(out) :: extPotAtom(:)

    !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
    real(dp), intent(out) :: extPotShell(:,:)

    call this%getExtPot(this%refPtr, chargePerAtom, extPotAtom)
    ! currently only atom resolved, no non-local l-dependent part
    extPotShell(:,:) = 0.0_dp

  end subroutine TDepExtPotGenC_getExternalPot


  !> Extra routine for interfacing gradients from a charge dependent external potential
  subroutine TQDepExtPotGenC_getExternalPotGrad(this, chargePerAtom, chargePerShell, extPotGrad)

    !> Class instance.
    class(TQDepExtPotGenC), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Gradient of the potential at each atom. Shape: [3, nAtom].
    real(dp), intent(out) :: extPotGrad(:,:)

    call this%getExtPotGrad(this%refPtr, chargePerAtom, extPotGrad)

  end subroutine TQDepExtPotGenC_getExternalPotGrad


  !> Extra routine for interfacing screened Fock external potential
  subroutine TDepExtPotGenC_getExternalPotFock(this, chargePerAtom, chargePerShell, extPotAtom, extPotShell)

    !> Class instance.
    class(TQDepExtPotGenC), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Screened Fock external potential
    real(dp), intent(out) :: extPotAtom(:)

    !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
    real(dp), intent(out) :: extPotShell(:,:)

    call this%getExtPotFock(this%refPtr, chargePerAtom, extPotAtom)

  end subroutine TDepExtPotGenC_getExternalPotFock


  !> Extra routine for interfacing screened Fock external potential
  subroutine TDepExtPotGenC_getInternalEnergy(this, energyInternal)

    !> Class instance.
    class(TQDepExtPotGenC), intent(inout) :: this

    !> QM/MM internal energy
    real(dp), intent(out) :: energyInternal

    call this%getIntEnergy(this%refPtr, energyInternal)

  end subroutine TDepExtPotGenC_getInternalEnergy


end module dftbp_dftbplus_qdepextpotgenc
