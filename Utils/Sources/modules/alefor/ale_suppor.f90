subroutine ale_suppor()
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_suppor
  ! NAME
  !    ale_suppor
  ! DESCRIPTION
  !    This routine applied boundary conditions
  ! OUTPUT 
  ! USES
  ! USED BY
  !    ale_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_alefor
  use mod_opebcs
  use mod_kdtree
  use def_kermod
  use mod_postpr
  use mod_messages, only : livinf
  implicit none
  integer(ip)  :: iboun,ipoin,idime
  real(rp)     :: dummr(3),dista,chkdi,proje(3),zeror1
  character(5) :: wopos(3)

  if( kfl_suppo_ale == 1 ) then

     if( ndivi > 0 ) then
        call livinf(0_ip,'PROJECT MM BOUNDARY MESH ONTO SUPPORT SURFACE',0_ip)
     else
        call livinf(0_ip,'PROJECT PART OF THE BOUNDARY MESH ONTO SUPPORT SURFACE',0_ip)
     end if

     if( INOTMASTER ) then
        !
        ! Allocate and construct KD-Tree
        !
        zeror1 = 1.0_rp + zeror
        chkdi  = 1.0e9_rp
        call kdtree(&
             1_ip,mnodb_ad,npoin_ad,nboun_ad,&
             coord_ad,lnodb_ad,ltypb_ad)
        !
        ! Look for minimum distance to the surface for all boundary nodes with fixity condition 2
        ! Displacement is BVESS_ALE on these nodes
        ! 
        do ipoin = 1,npoin
           if( kfl_fixno_ale(1,ipoin) == 2 ) then
              call dpopar(&
                   1_ip,coord(1:ndime,ipoin),npoin_ad,mnodb_ad,&
                   nboun_ad,chkdi,ltypb_ad,lnodb_ad,&
                   coord_ad,dista,dummr,proje,iboun) 
              do idime = 1,ndime
                 bvess_ale(idime,ipoin) = proje(idime)-coord(idime,ipoin)
              end do
           end if
        end do
        !
        ! Deallocate memory as we will no longer project onto the support surface
        !
        call kdtree(&
             2_ip,mnodb_ad,npoin_ad,nboun_ad,&
             coord_ad,lnodb_ad,ltypb_ad)
        call ale_membcs(-5_ip) 
     end if
     kfl_suppo_ale = -1
     !
     ! Postprocess
     !
     !wopos(1) = 'XXXXX' 
     !wopos(2) = 'VECTO'
     !wopos(3) = 'NPOIN'
     !call postpr(bvess_ale,wopos,ittim,cutim)   
     !do ipoin = 1,npoin_ad
     !   write(100,'(i7,3(1x,e12.6))') ipoin+npoin,coord_ad(1,ipoin),coord_ad(2,ipoin),coord_ad(3,ipoin)
     !end do
     !do iboun=1,nboun_ad
     !   write(101,'(10(1x,i7))') iboun+nelem,lnodb_ad(1:nnode(abs(ltypb_ad(iboun))),iboun)+npoin
     !end do
     !call runend('POPO')

  end if

end subroutine ale_suppor
