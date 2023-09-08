subroutine wav_updtss
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_updtss
  ! NAME 
  !    wav_updtss
  ! DESCRIPTION
  !    This routine computes the time step size for the temperature
  !    equation.
  ! USED BY
  !    wav_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_wavequ
  use mod_communications, only : PAR_MIN
  implicit none 
  integer(ip)      :: ielem,idime,inode,ipoin
  integer(ip)      :: pnode,pmate,pelty      
  real(rp)         :: dtcri,dtmin,rdinv,gpdet,gpden,gpkap,h,c
  real(rp)         :: elcod(ndime,mnode)
  real(rp)         :: gpcar(ndime,mnode),xjacm(9),xjaci(9)
  real(rp), target :: dtpar(1)

  dtmin = 1e6

  if( INOTMASTER ) then
     !
     ! Compute minimum element time step
     !
     rdinv = 1.0_rp/real(ndime)
     do ielem = 1,nelem
        pelty = ltype(ielem)
        pmate = lmate_wav(ielem)
        pnode = nnode(pelty)
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do idime = 1,ndime
              elcod(idime,inode) = coord(idime,ipoin)
           end do
        end do
        rdinv = 1.0_rp/real(ndime)
        call elmder(pnode,ndime,elmar(pelty)%dercg,elcod,gpcar,gpdet,xjacm,xjaci)
        h     = (elmar(pelty)%weicg*gpdet)**rdinv ! h
        gpden = densi_wav(1,pmate)                ! rho
        gpkap = kappa_wav(1,pmate)                ! kappa
        c     = sqrt(gpkap/gpden)                 ! c=sqrt(kappa/rho)
        dtcri = h/c                               ! dt=h/c
        dtmin = min(dtmin,dtcri)
     end do
  end if
  !
  ! Look for minimum over whole mesh
  !
  if( IPARALL) call PAR_MIN(dtmin)

  dtcri_wav = dtmin
  dtinv_wav = 1.0_rp/(dtcri_wav*safet_wav)
  if(kfl_timco==1) dtinv=max(dtinv,dtinv_wav) 

end subroutine wav_updtss
