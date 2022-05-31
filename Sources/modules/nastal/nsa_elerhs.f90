!------------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_elerhs.f90
!> @date    01/08/2012
!> @author  Mariano Vazquez
!> @brief   Computes the right-hand-side vector
!> @details Computes the right-hand-side vector on nodes, elrhs(ievat),\n
!!          using the values of the previous time-step.\n
!!          The right-hand-side vector is formed by the Galerkin (galer),\n
!!          stabilization (state) and shock capturing (shote) terms.!!\n
!> @}
!------------------------------------------------------------------------
subroutine nsa_elerhs(ielem,elrhs,inode,igaus,ndaux,dvolu,& 
     xvelo,xunkn,xconv,xdiff,dconv,ddiff,gunkn,dflux_conv,xshai,cartd,&
     hesma,xsube,shote,xvofo,heats,xtide)
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
  integer(ip) :: ndaux,jdofn,idofn,idime,jdime,kdime,ievat,inode,igaus,ielem
  real(rp)    :: xshai,dvolu,elrhs(nevat_nsa),state(ndofn_nsa),galte(ndofn_nsa),dtsub(ndofn_nsa), &
       xconv(ndofn_nsa,ndofn_nsa,ndime),gunkn(ndofn_nsa,ndime),cartd(ndime,mnode),ddiff(ndofn_nsa,ndofn_nsa,2,ndime), &
       xsube(ndofn_nsa,mgaus,3),xdiff(ndofn_nsa,ndofn_nsa,ndime,ndime),hesma(ndime,ndime,mnode), &
       xunkn(ndofn_nsa,mgaus,3),xvelo(ndime),shote(ndofn_nsa),dconv(ndofn_nsa,ndofn_nsa), dflux_conv(ndofn_nsa),&
       xtide(ndofn_nsa) ! Preconditioner * ( Unknown(last_iteration) - Unknown(last_time_step) ) / physical_time_step

  real(rp) :: xx_rtemp,xx_rtempC,xx_rtempG,xx_rtempH,xx_rtempX

  real(rp),      intent(in) :: xvofo(ndofn_nsa,ndofn_nsa),heats

  !-----------------------------------------------------------------------!
  ! LODI
  !-----------------------------------------------------------------------! 
  integer(ip) :: ipoin 
  real(rp)    :: xsour(ndofn_nsa), xaux(ndofn_nsa)
  xsour(1:ndofn_nsa) = 0.0_rp
  !-----------------------------------------------------------------------! 

  state= 0.0_rp
  galte= 0.0_rp
  dtsub= 0.0_rp

  do idime=1,ndime
     do idofn=1,ndofn_nsa
        xx_rtempG= gunkn(idofn,idime)  * xshai 
        do jdofn=1,ndaux  
           galte(jdofn) = galte(jdofn) + xconv(jdofn,idofn,idime) * xx_rtempG
        end do
     end do
  enddo
  
  if (kfl_skews_nsa == 1) then
     ! skew symmetric
     galte = galte + xshai * dflux_conv
  end if

  if (kfl_stabi_nsa >= 1) then
     do idime=1,ndime
        xx_rtempC=cartd(idime,inode)
        do idofn=1,ndofn_nsa
           xx_rtempX = xsube(idofn,igaus,1) * xx_rtempC
           do jdofn=1,ndaux  
              if (kfl_lopre_nsa < 2) then
                 state(jdofn) = state(jdofn) + (xconv(jdofn,idofn,idime) + ddiff(jdofn,idofn,2,idime)) * xx_rtempX
              else
                 ! (MM) This is an approximation. ddiff should appear here as in the 'kfl_lopre_nsa = 0' case
                 state(jdofn) = state(jdofn) + (xconv(jdofn,idofn,idime)) * xx_rtempX 
              end if
           end do
        end do
     enddo
  end if


  !-----------------------------------------------------------------------!
  ! LODI
  !-----------------------------------------------------------------------! 
  if(lodi_nsa==1) then
     galte = 0.0_rp 
     ipoin = lnods(inode,ielem)
     if(normal_nsa(ipoin)%id>0) then
        do idime = 1,ndime
           xaux( 1:ndofn_nsa) = matmul( Sx_nsa(1:ndofn_nsa,1:ndofn_nsa,idime,igaus,ielem), xchrc_nsa(1:ndofn_nsa,ielem,igaus,idime) )
           xaux( normal_nsa(ipoin)%ichrc ) = kfact_lodi_nsa * ( press(ipoin,1) - prefe_lodi_nsa ) * sound_nsa(ipoin) 
           xsour(1:ndofn_nsa) = xsour(1:ndofn_nsa) - (/xaux(2:ndime+1),xaux(1),xaux(ndime+2)/)
        enddo
     else
        do idime = 1,ndime
           galte(1:ndofn_nsa) = galte(1:ndofn_nsa) + matmul( xconv(1:ndofn_nsa,1:ndofn_nsa,idime), gunkn(1:ndofn_nsa,idime) * xshai )
        enddo
     endif
  endif
  !-----------------------------------------------------------------------! 

  do idofn=1,ndofn_nsa
     xx_rtempX=xshai*xsube(idofn,igaus,1)
     do jdofn=1,ndaux
        if (kfl_lopre_nsa < 2) then
           state(jdofn) = state(jdofn) + (dconv(jdofn,idofn) + xvofo(jdofn,idofn))*xx_rtempX
        else
           ! (MM) This is an approximation. dconv should appear here as in the 'kfl_lopre_nsa = 0' case
           state(jdofn) = state(jdofn) + xvofo(jdofn,idofn)*xx_rtempX 
        end if
        galte(idofn) = galte(idofn) - xvofo(idofn,jdofn)*xshai*xunkn(jdofn,igaus,1)
     enddo
  end do

  do jdime=1,ndime
     do idime=1,ndime
        xx_rtempC = cartd(idime,inode) 
        xx_rtempH = hesma(jdime,idime,inode)
        do idofn=1,ndofn_nsa
           xx_rtempG = gunkn(idofn,jdime)  * xx_rtempC
           xx_rtempX =  xsube(idofn,igaus,1) * xx_rtempH
           do jdofn=1,ndaux  
              galte(jdofn) = galte(jdofn) + xdiff(jdofn,idofn,idime,jdime) * xx_rtempG
              state(jdofn) = state(jdofn) + xdiff(jdofn,idofn,idime,jdime) * xx_rtempX 
           end do
        end do
     end do
  enddo

  galte = - galte

  galte(ndime+2) = galte(ndime+2) + heats * xshai  ! Addition of heat source from chemical reactions

  if (kfl_stabi_nsa == 0) state= 0.0_rp     


  do jdofn=1,ndaux  
     ievat= (inode-1) * ndofn_nsa + jdofn
     ! shote is negative because how is it made in nsa_shocap subroutine
     elrhs(ievat) = elrhs(ievat) + dvolu * &
          (-xtide(jdofn)*xshai + galte(jdofn) + state(jdofn) - shote(jdofn) + dtsub(jdofn))
     elrhs(ievat) = elrhs(ievat) + dvolu * xshai * xsour(jdofn)
  end do

!!$if (ielem == 1 .and. igaus == 1 .and. inode == 1) then
!!$ievat= (inode-1) * ndofn_nsa
!!$print*
!!$print*, 'dtinv',itinn(modul),1.0_rp/dtinv_nsa,1.0_rp/dtinv
!!$print*, 'elrhs',elrhs(ievat+1),elrhs(ievat+2),elrhs(ievat+3),elrhs(ievat+4)
!!$print*, 'galte',galte(1),galte(2),galte(3),galte(4)
!!$print*, 'state',state(1),state(2),state(3),state(4)
!!$end if


end subroutine nsa_elerhs
