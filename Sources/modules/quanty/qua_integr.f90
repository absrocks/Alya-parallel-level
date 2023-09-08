subroutine qua_integr(itask)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_integr
  ! NAME 
  !    qua_integr
  ! DESCRIPTION
  !    This routine make all the integrals over the entire domain
  !    for calculate total energy in DFT.
  ! USES
  ! 
  ! USED BY
  !    qua_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_quanty
  use def_domain
  implicit none
  real(rp)    :: Ecin,Ecoul,Ehartree,Exc,ah,ah1,ah2,intrho
  integer(ip) :: kpos,kk,jj,nk,nderiv,itask

  nk     = 0
  nderiv = 0
  
  call qua_3dinte_R(NK,NDERIV,RHOON,intrho)
  
  !call qua_3dinte(NK,NDERIV,rhoon(:,1),intrho)  ! integro la densidad final
  
  if( INOTSLAVE ) write(6,*) 'densidad: ',intrho

  if( INOTSLAVE ) write(eigen_sol(1)%lun_solei ,*) 'integral de Rho  ',intrho

  if( itask == 1 ) return

  call integraGauss(Ecin,Ecoul,Ehartree,Exc)   
   
  if( INOTSLAVE ) then
     write(6,*) 'Carga  =   ',intrho
     WRITE(6,*) 'Ecin   =   ',Ecin,'  Ecoul= ',Ecoul,' Ehartree= ',Ehartree*0.5_rp,' Exc= ',Exc
     WRITE(6,*) 'Etot   =   ',Ecin + Ecoul + Ehartree*0.5_rp + Exc
     !WRITE(6,*) 'Etot(2)=   ',noccupa(1)*eigva(1)- EHARTREE*0.5 - EXC*1./3.
     
     write(eigen_sol(1)%lun_solei ,*) 'Carga  =   ',intrho
     WRITE(eigen_sol(1)%lun_solei ,*) 'Ecin   =   ',Ecin,'  Ecoul= ',Ecoul,' Ehartree= ',Ehartree*0.5_rp,' Exc= ',Exc
     WRITE(eigen_sol(1)%lun_solei ,*) 'Etot   =   ',Ecin + Ecoul + Ehartree*0.5_rp + Exc
     !WRITE(eigen_sol(1)%lun_solei ,*) 'Etot(2)=   ',noccupa(1)*eigva(1)- EHARTREE*0.5 - EXC*1./3.
  end if

  contains

    subroutine integraGauss(Ecin,Ecoul,Ehartree,Exc)
      implicit none

      real(rp), intent(inout) :: Ecin,Ecoul,Ehartree,Exc
      real(rp)    :: thrd_1,aa_1,bb_1,ax_1,xcoef_1,ycoef_1,zcoef_1,xyln
      integer(ip) :: ielem,igaus,inode,ipoin,idime,ki          ! Indices and dimensions
      integer(ip) :: pelty,pmate,pnode
      integer(ip) :: pgaus,plapl,porde,ptopo

      real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
      real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj

      real(rp)    :: elcod(ndime,mnode),gpcod(ndime,mgaus)
      real(rp)    :: gpvol(mgaus)                          ! |J|*w
      real(rp)    :: ecgaus(mgaus), ehgaus(mgaus), exgaus(mgaus), akgaus(mgaus)
      real(rp)    :: elrho(mnode),  elcoul(mnode), elvhar(mnode),akcin(mnode),elvxc(mnode)
      real(rp)    :: gpsha(mnode,mgaus),phi_cin(mnode,nestates),auxi,xderi(3)
      real(rp)    :: dummr(4)

      thrd_1 = 0.333333333333333333_rp
      aa_1   = 0.93222_rp
      bb_1   = 0.947362E-2_rp
      ax_1   = 0.738558766382022406_rp

      Ecoul    = 0.0_rp
      Ehartree = 0.0_rp
      Exc      = 0.0_rp
      Ecin     = 0.0_rp

      if( INOTMASTER ) then

         do ielem = 1,nelem 


            pelty = ltype(ielem)
            pnode = nnode(pelty)
            pgaus = ngaus(pelty)
            plapl = 0
            porde = lorde(pelty)
            ptopo = ltopo(pelty)
            !
            ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
            !

            do inode=1,pnode
               ipoin=lnods(inode,ielem)
               do idime=1,ndime
                  elcod(idime,inode)=coord(idime,ipoin)
               end do
            end do

            call elmcar(&
                 pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                 elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
                 gphes,ielem)

            do inode = 1,pnode

               ipoin         = lnods(inode,ielem)
               elrho(inode)  = rhoon(ipoin,1)
               elcoul(inode) = v_pot_ps(ipoin)
               elvhar(inode) = v_hartree(ipoin)

               if( kfl_potxc_qua == 0 ) then 

                  elvxc(inode)= 0.0_rp

               elseif( kfl_potxc_qua == 1 ) then


                  xcoef_1 =  bb_1*elrho(inode)**thrd_1
                  ycoef_1 =  xcoef_1/(xcoef_1 + 1.0_rp)
                  zcoef_1 = -aa_1*xcoef_1/bb_1
                  xyln    =  xcoef_1*log(ycoef_1)

                  elvxc(inode)=zcoef_1 * (1.0_rp + xyln)
               endif

               do kk = 1,nestates
                  phi_cin(inode,kk) = eigen( (kk-1)*npoin + ipoin )
               enddo

            end do


            do igaus = 1,pgaus

               ecgaus(igaus) = 0.0_rp
               ehgaus(igaus) = 0.0_rp
               exgaus(igaus) = 0.0_rp
               akgaus(igaus) = 0.0_rp

               do inode=1,pnode
                  ecgaus(igaus) = ecgaus(igaus) + elcoul(inode) * elrho(inode) * elmar(pelty) % shape(inode,igaus)
                  ehgaus(igaus) = ehgaus(igaus) + elvhar(inode) * elrho(inode) * elmar(pelty) % shape(inode,igaus)
                  exgaus(igaus) = exgaus(igaus) + elvxc(inode)  * elrho(inode) * elmar(pelty) % shape(inode,igaus)
               end do

               do ki = 1,nestates
                  xderi(1) = 0.0_rp
                  xderi(2) = 0.0_rp
                  xderi(3) = 0.0_rp
                  do inode = 1,pnode
                     xderi(1) = xderi(1) + phi_cin(inode,ki) * gpcar(1,inode,igaus)
                     xderi(2) = xderi(2) + phi_cin(inode,ki) * gpcar(2,inode,igaus)
                     xderi(3) = xderi(3) + phi_cin(inode,ki) * gpcar(3,inode,igaus)
                  end do
                  auxi = xderi(1)*xderi(1) + xderi(2)*xderi(2) + xderi(3)*xderi(3)
                  akgaus(igaus) = akgaus(igaus) + 0.5_rp * noccupa(ki) * auxi
               end do

               Ecoul   = Ecoul    + gpvol(igaus) * ecgaus(igaus)
               Ehartree= Ehartree + gpvol(igaus) * ehgaus(igaus)
               Exc     = Exc      + gpvol(igaus) * exgaus(igaus)
               Ecin    = Ecin     + gpvol(igaus) * akgaus(igaus)

            end do

         end do

      end if

      dummr(1) = Ecoul 
      dummr(2) = Ehartree
      dummr(3) = Exc     
      dummr(4) = Ecin    
      call pararr('SUM',0_ip,4_ip,dummr)
      Ecoul    = dummr(1)
      Ehartree = dummr(2)
      Exc      = dummr(3)
      Ecin     = dummr(4)

    end subroutine integraGauss
  
end subroutine qua_integr

