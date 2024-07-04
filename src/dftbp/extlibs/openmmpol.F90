!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"
#:include "error.fypp"

module dftbp_extlibs_openmmpol
#:if WITH_OPENMMPOL
  use ommp_interface, only : ommp_system, ommp_init_mmp, ommp_init_xyz, ommp_terminate,&
      & ommp_print_summary, ommp_prepare_qm_ele_ene, ommp_set_external_field,&
      & ommp_potential_mmpol2ext, ommp_qm_helper, ommp_init_qm_helper, ommp_terminate_qm_helper,&
      & ommp_set_verbose, ommp_get_full_ele_energy, OMMP_VERBOSE_DEBUG, OMMP_VERBOSE_LOW,&
      & OMMP_VERBOSE_HIGH, ommp_get_fixedelec_energy, ommp_get_polelec_energy,&
      & ommp_get_full_bnd_energy, ommp_get_vdw_energy, ommp_qm_helper_init_vdw_prm,&
      & ommp_qm_helper_vdw_energy, ommp_qm_helper_set_attype, ommp_qm_helper_vdw_energy,&
      & ommp_print_summary_to_file, OMMP_SOLVER_NONE, OMMP_FF_AMOEBA, OMMP_FF_WANG_AL,&
      & OMMP_FF_WANG_DL
#:endif
  use dftbp_io_message, only : error, warning
  use dftbp_common_environment, only : TEnvironment
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_dftb_charges, only : getSummedCharges
  use dftbp_common_accuracy, only : dp, mc

  public TOMMPInterface, TOMMPInterface_init

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


end module