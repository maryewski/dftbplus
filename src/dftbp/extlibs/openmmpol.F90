!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"
#:include "error.fypp"

module dftbp_extlibs_openmmpol
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_dftbplus_qdepextpotgen, only : TQDepExtPotGen
  use dftbp_io_message, only : error, warning
#:if WITH_OPENMMPOL
  use ommp_interface
#:endif

  implicit none

  private
  public TOpenmmpolPotGen, TOpenmmpolPotGen_init, TOpenmmpolInput

  ! Constant for AMBER identification inside openmmpol;
  ! TODO: import from openmmpol
  integer, protected :: OMMP_FF_AMBER = 0

  !> Interface to control polarizable QM/MM
  !! through openmmpol external library
  type, extends(TQDepExtPotGen) :: TOpenmmpolPotGen
    private
#:if WITH_OPENMMPOL
    !> Pointer to openmmpol system object
    type(ommp_system), pointer :: pSystem
    !> Pointer to openmmpol QM helper object
    type(ommp_qm_helper), pointer :: pQMHelper
#:endif
    !> Potential to calculate total energy
    real(dp), allocatable :: potential(:)
    !> Effective potential to calculate Fock matrix elements;
    !! normally is equal to the ordinary potential,
    !! but in cases like AMOEBA it includes extra contributions
    !! from screened dipole-dipole interaction, which happen
    !! do not contribute to QM/MM interaction energy.
    real(dp), allocatable :: potentialFock(:)
    !> Total energy of all openmmpol-specific
    !! contributions
    real(dp) :: Einternal

  contains
    procedure :: init => TOpenmmpolPotGen_init
    !> method for evaluating external potential
    procedure :: getExternalPot => TOpenmmpolPotGen_getExternalPot
    !> method for evaluating screened Fock potential
    procedure :: getExternalPotFock => TOpenmmpolPotGen_getExternalPotFock
    !> method for evaluating gradient of external potential
    procedure :: getExternalPotGrad => TOpenmmpolPotGen_getExternalPotGrad
    !> method for adding internal energy contribution
    procedure :: getInternalEnergy => TOpenmmpolPotGen_getInternalEnergy
    !> method for writing openmmpol-specific data
    procedure :: writeOutput => TOpenmmpolPotGen_writeOutput
    !> method for updating coordinates in the QM zone
    ! procedure :: updateQMCoords => TOpenmmpolPotGen_updateQMCoords
    ! TODO: add methods for controlling dynamics and excited states contributions
  end type

  !> Data type for storing openmmpol-specific input parameters
  type :: TOpenmmpolInput
    !> Used input format; allowed values:
    !! "Tinker", "mmp"
    character(:), allocatable :: inputFormat
    !> Index of linear solver
    integer :: solver
    !> Path to MM geometry file
    character(:), allocatable :: mmGeomFilename
    !> Path to a separate parameter file, if present
    character(:), allocatable :: mmParamsFilename
    !> MM atom types for atoms in the QM zone;
    !! used for vdW parameter access
    integer, allocatable :: qmAtomTypes(:)
    !> Path to parameter file containing MM atom
    !! types for atoms in the QM zone
    character(:), allocatable :: qmParamsFilename
  end type

contains

  subroutine TOpenmmpolPotGen_writeOutput(this)
    class(TOpenmmpolPotGen) :: this

#:if WITH_OPENMMPOL
    call notImplementedError
#:else
    call notImplementedError
#:endif
  end subroutine TOpenmmpolPotGen_writeOutput


  subroutine TOpenmmpolPotGen_init(this, input)
    class(TOpenmmpolPotGen), intent(inout) :: this
    type(TOpenmmpolInput), intent(in) ::  input
    
#:if WITH_OPENMMPOL
    call ommp_print_summary_to_file(this%pSystem, "openmmpol.out")
#:else
    call notImplementedError
#:endif
  end subroutine TOpenmmpolPotGen_init


  subroutine TOpenmmpolPotGen_getExternalPot(this, chargePerAtom, chargePerShell, extPotAtom, extPotShell)
    !> Instance.
    class(TOpenmmpolPotGen), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Atom dependent external potential contribution. Shape: [nAtom]
    real(dp), intent(out) :: extPotAtom(:)

    !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
    real(dp), intent(out) :: extPotShell(:,:)

#:if WITH_OPENMMPOL
    ! TODO: implement adding potential
#:else
    call notImplementedError
#:endif
  end subroutine TOpenmmpolPotGen_getExternalPot


  subroutine TOpenmmpolPotGen_getExternalPotGrad(this, chargePerAtom, chargePerShell, extPotGrad)
    !> Class instance.
    class(TOpenmmpolPotGen), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Gradient of the potential at each atom. Shape: [3, nAtom].
    real(dp), intent(out) :: extPotGrad(:,:)

#:if WITH_OPENMMPOL
    call notImplementedError
#:else
    call notImplementedError
#:endif
  end subroutine TOpenmmpolPotGen_getExternalPotGrad


  subroutine TOpenmmpolPotGen_getExternalPotFock(this, chargePerAtom, chargePerShell, extPotAtom, extPotShell)
    !> Instance.
    class(TOpenmmpolPotGen), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Atom dependent external potential contribution. Shape: [nAtom]
    real(dp), intent(out) :: extPotAtom(:)

    !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
    real(dp), intent(out) :: extPotShell(:,:)

#:if WITH_OPENMMPOL
    ! TODO: implement adding potential
#:else
    call notImplementedError
#:endif
  end subroutine TOpenmmpolPotGen_getExternalPotFock


  subroutine TOpenmmpolPotGen_getInternalEnergy(this, energyInternal)
    !> Class instance.
    class(TOpenmmpolPotGen), intent(inout) :: this

    !> External energy contribution
    real(dp), intent(out) :: energyInternal

#:if WITH_OPENMMPOL
    energyInternal = energyInternal + this%Einternal
#:else
    call notImplementedError
#:endif
  end subroutine TOpenmmpolPotGen_getInternalEnergy

  
  subroutine notImplementedError
    call error("DFTB+ compiled without support for openmmpol library")
  end subroutine notImplementedError
end module