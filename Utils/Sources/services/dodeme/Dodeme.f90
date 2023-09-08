subroutine Dodeme(order)
  !-----------------------------------------------------------------------
  !****f* dodeme/Dodeme
  ! NAME 
  !    Dodeme
  ! DESCRIPTION
  !    This routine is the bridge for Dodeme service
  ! USES
  !
  ! USED BY
  !
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_dodeme
  use def_elmtyp
  use mod_graphs
  use mod_outfor, only : outfor
  use mod_messages, only : livinf
  implicit none
  integer(ip), intent(in) :: order
  real(rp)                :: time0,time3,time1,time2

  servi = ID_DODEME

  if( ISLAVE ) return

  call cputim(time0)

  if( order == 0 ) then

     !-------------------------------------------------------------------
     !
     ! Read data
     !
     !-------------------------------------------------------------------

     call dod_reapro() 
 
  else if( kfl_servi(servi) == 1 ) then 

     !-------------------------------------------------------------------
     !
     ! Create new element: Extension, Dirichlet, Neumann, Hole
     !
     !-------------------------------------------------------------------

     if( order == 1 ) then

        call livinf(-4_ip,'DODEME',0_ip)

        !----------------------------------------------------------------
        !
        ! Prepare arrays
        !
        !----------------------------------------------------------------
        ! 
        ! Read data
        !
        call dod_turnon() 
        !
        ! Copy original mesh to subdomain wise mesh
        !
        call cputim(time1)
        call dod_copy_mesh()
        call cputim(time2) ; cpu_dod_copy_mesh = time2 - time1
        !
        ! Compute some auxiliary variables
        !
        call dod_inivar(2_ip)
        !
        ! Compute some graphs
        !
        call cputim(time1)
        call dod_graphs(1_ip)
        call cputim(time2) ; cpu_dod_graphs = time2 - time1
        !
        ! Prepare std-kdtree structure
        !
        call cputim(time1)
        call dod_kdtree(1_ip)
        call cputim(time2) ; cpu_dod_kdtree = time2 - time1
        !
        ! Prepare elsest structure
        !
        call cputim(time1)
        call dod_elsest()
        call cputim(time2) ; cpu_dod_elsest = time2 - time1

        !----------------------------------------------------------------
        !
        ! Treat the 3 types of topologies:
        ! 1. Chimera: identify SUBDOMAIN(ISUBD) % LBOCH(iboun) = BOEXT  
        ! 2. Patch: identify SUBDOMAIN(ISUBD) % LBOCH(iboun) = BOEXT  
        ! 3. Prescribed boundaries: compute PRESCRIBED_BOUNDARIES(IBOUN)
        !
        !    For these three topologies we end up with SUBDOMAIN(ISUBD) % LSUBD_NPOIN(IPOIN):
        !
        !                      Free   hole   fringe
        !          Chimera      0   -JSUBD   JSUBD
        !          Patch        0      x     JSUBD
        !          Prescribed   0      x     JSUBD
        !
        !----------------------------------------------------------------        
        !
        ! Bakcground: Hole cutting if IHOLE_DOD(ISUBD) = 1
        !
        call cputim(time1)
        call dod_holcut()
        call cputim(time2) ; cpu_dod_holcut = time2 - time1
        !
        ! Patch: detect boundary automatically if IPATC_DOD(ISUBD) = 1
        !
        call cputim(time1)
        call dod_patch_boundaries()
        call cputim(time2) ; cpu_dod_patch_boundaries = time2 - time1
        !
        ! Prescribed interfaces: detect boundary automatically if IPRES_DOD(ISUBD) = 1
        !
        call cputim(time1)
        call dod_prescribed_boundaries()
        call cputim(time2) ; cpu_dod_prescribed_boundaries = time2 - time1

        !----------------------------------------------------------------
        !
        ! Allocate memory for fringe node extension elements and
        ! look for candidates
        !
        !----------------------------------------------------------------        

        call cputim(time1)
        call dod_candidate_nodes()
        call cputim(time2) ; cpu_dod_candidate_nodes = time2 - time1

        !----------------------------------------------------------------
        !
        ! Create extensions
        !
        !----------------------------------------------------------------        

        call cputim(time1)
        call dod_extension()
        call cputim(time2) ; cpu_dod_extension = time2 - time1

        !----------------------------------------------------------------
        !
        ! Merge elements
        !
        !----------------------------------------------------------------        

        call cputim(time1)
        call dod_merge_new_elements()
        call cputim(time2) ; cpu_dod_merge_new_elements = time2 - time1

        !----------------------------------------------------------------
        !
        ! Deallocate memory and finish dodeme
        !
        !----------------------------------------------------------------        

        call dod_memall(-2_ip)
        !
        ! Output cpu times
        !
        call dod_outcpu()
        !
        ! Write tail for formatted files
        !
        call outfor(26_ip,lun_outpu_dod,' ')
        !
        ! Close output file
        !
        call dod_openfi(3_ip)

        call livinf(-5_ip,'END DODEME',0_ip)

     end if

  end if

  call cputim(time3)
  cpu_servi(1,servi) =  cpu_servi(1,servi) + (time3 - time0)


end subroutine Dodeme
