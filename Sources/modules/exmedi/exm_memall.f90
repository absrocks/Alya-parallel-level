subroutine exm_memall
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_memall
  ! NAME 
  !    exm_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    equations of EXMEDI
  !    Allocation is done for:
  !      - Potentials (intra- and extracellular, and transmembrane)
  !      - Subcell and no-subcell models of ionic currents
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !    exm_turnon
  !***
  !-----------------------------------------------------------------------
  !
  ! DESCRIPTION OF VARIABLES' MEANING 
  !
  !-----------------------------------------------------------------------
  !
  ! Variable: elmag
  ! 
  !     Potentials (i.e. the unknowns) are stored in elmag(npoin,ncomp_exm):
  !        - internal action potential (v) ---> elmag(:,:)
  !        - external action potential (u) ---> extac_exm(:,:) TO BE DONE
  !        - FHN recovery potential        ---> refhn_exm(:,:)
  !
  !-----------------------------------------------------------------------
  !
  ! Variable: vauxi_exm
  !
  !     vauxi_exm(nauxi_exm,npoin,ncomp_exm)
  !        It is used in exm_ceauxi.f90 subroutine for computation of activation/
  !        inactivation states corresponding to Hodgkin-Huxley formalism
  !
  !-----------------------------------------------------------------------
  !
  ! Variable: 

  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      def_solver
  use      def_exmedi
  use      mod_exm_arrays, only : exm_arrays
  use mod_memory

  implicit none
  !integer(ip)    :: nauxi
  !integer(4)     :: istat,imate


  call exm_arrays('ALLOCATE')

end subroutine exm_memall
