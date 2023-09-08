subroutine got_updtss()
!-----------------------------------------------------------------------
!****f* Gotita/got_updtss
! NAME 
!    got_updtss
! DESCRIPTION
!    This routine computes the time step size 
!    equation. 
! USED BY
!    got_begste
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  use mod_communications, only : PAR_MIN
  implicit none 
  integer(ip)      :: ielem,idime,inode,ipoin
  integer(ip)      :: pnode,pelty      
  real(rp)         :: dtcri

  if(kfl_timei_got/=0) then
     !
     ! Loop over elements
     !
     if(kfl_paral/=0) then

        timin_got = 1e30
        timax_got =-1e30
        dtmin(1)  = 1e6
        do ielem = 1,nelem
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           call got_elmgat(&
                2_ip,1_ip,pnode,lnods(1,ielem),elvdr_got,elcdr_got,&
                elcod_got,elvel_got,eldif_got)
           call got_elmtss(&
                pnode,lorde(pelty),elcod_got,elvel_got,elvdr_got,elcdr_got,&
                eldif_got,elmar(pelty)%shacg,elmar(pelty)%dercg,&
                elmar(pelty)%weicg,dtcri)           
           dtmin(1)=min(dtmin(1),dtcri)
           if(dtmin(1)>timax_got) timax_got=dtmin(1)
           if(dtmin(1)<timin_got) timin_got=dtmin(1)
        end do
        
     end if
     !
     ! Parall: Look for minimum over all subdomains (dtmin)
     !
     if(kfl_paral>=0) call PAR_MIN(dtmin)
     !
     ! Assign 1/dt
     !
     dtcri_got = dtmin(1)
     if(dtcri_got/=0.0_rp) dtinv_got = 1.0_rp/(dtcri_got*safet_got)
     if(kfl_timco==1)      dtinv=max(dtinv,dtinv_got)
        
  end if

end subroutine got_updtss
