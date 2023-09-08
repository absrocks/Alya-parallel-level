subroutine got_diffus()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_diffus
  ! NAME 
  !    got_diffus
  ! DESCRIPTION
  !    This routine smoothes the diffusion
  ! USES
  ! USED BY
  !    got_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_domain
  use def_master
  use def_gotita
  implicit none
  integer(ip) :: ipoin,inode,ielem,idime,pnode,pelty,jpoin
  real(rp)    :: gpdet,gpvol,gpcar(ndime,mnode)
  real(rp)    :: xjaci(9),xjacm(9),fact1,dummr
  real(rp)    :: tragl(9),chale(2),hleng(3)

  !integer(ip), save :: ipass=0
  !if(ipass==0) then
  !   ipass=1
  !   do ipoin=1,npoin
  !      gpdet=1e6
  !      do jpoin=1,npoin
  !         if(lpoty(jpoin)/=0) then
  !            if(  .not.(coord(1,jpoin)<xmini_got(1).or.coord(1,jpoin)>xmaxi_got(1).or.&
  !                 &     coord(2,jpoin)<xmini_got(2).or.coord(2,jpoin)>xmaxi_got(2))) then  
  !               gpvol=0.0_rp
  !               do idime=1,ndime
  !                  gpvol=gpvol+(coord(idime,ipoin)-coord(idime,jpoin))**2
  !               end do
  !               if(gpvol<gpdet) gpdet=gpvol
  !            end if
  !         end if
  !      end do
  !      gpdet=sqrt(gpdet)
  !      if(gpdet<0.1_rp) then
  !         diffm_got(ipoin)=1.0e-6
  !      else
  !         diffm_got(ipoin)=1.0e-8
  !      end if
  !   end do
  !end if
  !return

  if(kfl_diffu_got==2.and.kfl_probl_got/=2) then
     !
     ! Initialization
     !
     do ipoin=1,npoin
        diffm_got(ipoin)=0.0_rp
     end do
     !
     ! Compute LHS
     !
     do ielem=1,nelem
        pelty=ltype(ielem) 
        pnode=nnode(pelty)
        call got_elmgat(&
             3_ip,1_ip,pnode,lnods(1,ielem),elvdr_got,&
             elcdr_got,elcod_got,elvel_got,eldif_got)
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           call elmder(&
                pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                elcod_got,gpcar,gpdet,xjacm,xjaci)
           gpvol=elmar(pelty)%weigc(inode)*gpdet
           call elmlen(&
                ndime,pnode,elmar(pelty)%dercg,tragl,elcod_got,&
                hnatu(pelty),hleng)
           call elmchl(&
                dummr,dummr,elcod_got,elvdr_got,dummr,chale,pnode,&
                lorde(pelty),1.0_rp,1.0_ip,5_ip)
           call got_diffun(elcdr_got,chale,fact1)
           diffm_got(ipoin)=diffm_got(ipoin)+gpvol*fact1
        end do
     end do
     !
     ! Solve diagonal system
     !
     do ipoin=1,npoin
        diffm_got(ipoin)=diffm_got(ipoin)/vmasc(ipoin)
     end do

  end if

end subroutine got_diffus
