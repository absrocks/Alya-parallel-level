subroutine nsa_elerhs_cdr(elrhs,inode,igaus,ndaux,dvolu,& 
     xvelo,xunkn,xconv,xdiff,dconv,ddiff,gunkn,xshai,cartd,hesma,xsube,shote)
  !-----------------------------------------------------------------------
  !****f* nastal/nsa_elerhs
  ! NAME 
  !    nsa_elerhs
  ! DESCRIPTION
  !    Computes the elrhs vector
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_nastal
  implicit none
  integer(ip) :: ndaux,jdofn,idofn,idime,jdime,kdime,ievat,inode,igaus
  real(rp)    :: xshai,dvolu,elrhs(nevat_nsa),state(ndofn_nsa),galte(ndofn_nsa),dtsub(ndofn_nsa), &
       xconv(ndofn_nsa,ndofn_nsa,ndime),gunkn(ndofn_nsa,ndime),cartd(ndime,mnode),ddiff(ndofn_nsa,ndofn_nsa,2,ndime), &
       xsube(ndofn_nsa,mgaus,3),xdiff(ndofn_nsa,ndofn_nsa,ndime,ndime),hesma(ndime,ndime,mnode), &
       xunkn(ndofn_nsa,mgaus,3),xvelo(ndime),shote(ndofn_nsa),dconv(ndofn_nsa,ndofn_nsa)
  real(rp) :: xx_rtemp,xx_rtempC,xx_rtempG,xx_rtempH,xx_rtempX
  
  state= 0.0_rp
  galte= 0.0_rp
  dtsub= 0.0_rp
  
  do idofn=1,ndofn_nsa
     xx_rtempX=xshai*xsube(idofn,igaus,1)
     do jdofn=1,ndaux  
        state(jdofn) = state(jdofn) + dconv(jdofn,idofn)*xx_rtempX
     enddo
  end do
  
  do idime=1,ndime
     xx_rtempC=cartd(idime,inode)
     do idofn=1,ndofn_nsa
        xx_rtempG=gunkn(idofn,idime)  * xshai
        xx_rtempX= xsube(idofn,igaus,1) * xx_rtempC
        do jdofn=1,ndaux  
           galte(jdofn) = galte(jdofn) + xconv(jdofn,idofn,idime) * xx_rtempG
           state(jdofn) = state(jdofn) + (xconv(jdofn,idofn,idime) + ddiff(jdofn,idofn,2,idime)) * xx_rtempX
        end do
        
     end do
  end do
  
  do jdime=1,ndime
     do idime=1,ndime
        xx_rtempC = cartd(idime,inode) 
        xx_rtempH = hesma(jdime,idime,inode)
        do idofn=1,ndofn_nsa
           xx_rtempG = gunkn(idofn,jdime) * xx_rtempC
           xx_rtempX=  xsube(idofn,igaus,1) * xx_rtempH
           do jdofn=1,ndaux  
              galte(jdofn) = galte(jdofn) + xdiff(jdofn,idofn,idime,jdime) * xx_rtempG
              state(jdofn) = state(jdofn) + xdiff(jdofn,idofn,idime,jdime) * xx_rtempX 
           end do
        end do
     end do
  enddo
  
  galte = - galte
  
  if (kfl_stabi_nsa == 0) state= 0.0_rp


!!!!!!!!!!!!!!!!!!!!!!!
  galte = 0.0_rp
  state = 0.0_rp
  state(ndime+1) = - react_nsa * xshai * xsube(ndime+1,igaus,1)
  do idime=1,ndime
     galte(ndime+1) = galte(ndime+1) - (xshai * (conve_nsa(idime) + react_nsa) &
          + diffu_nsa * cartd(idime,inode)) * gunkn(ndime+1,idime)
     state(ndime+1) = state(ndime+1) + (conve_nsa(idime) * cartd(idime,inode) &
          + diffu_nsa * hesma(idime,idime,inode)) * xsube(ndime+1,igaus,1)
  end do

!!$print*,'galte gunkn cartd dvolu',galte(ndime+1),gunkn(ndime+1,1),gunkn(ndime+1,2),cartd(1,inode),cartd(2,inode),dvolu
!!!!!!!!!!!!!!!!!!!!!!!!
    
  do jdofn=1,ndaux  
     ievat= (inode-1) * ndofn_nsa + jdofn
     elrhs(ievat) = elrhs(ievat) + dvolu * (galte(jdofn) + state(jdofn) + shote(jdofn) + dtsub(jdofn))
  end do
!!$print*,'elrhs',elrhs((inode-1) * ndofn_nsa + 1),elrhs((inode-1) * ndofn_nsa + 2), &
!!$     elrhs((inode-1) * ndofn_nsa + 3),elrhs((inode-1) * ndofn_nsa + 4) 
end subroutine nsa_elerhs_cdr
