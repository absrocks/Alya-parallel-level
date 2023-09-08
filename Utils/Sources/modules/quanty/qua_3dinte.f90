subroutine qua_3dinte(NK,NDERIV,funcion,rta)
  !-----------------------------------------------------------------------
  !****f* 3dinte
  ! NAME 
  !    3dinte
  ! DESCRIPTION
  !    This routine make an three-dimensional integral over the domain
  ! INPUT
  !   funcion (*) ?
  !  NK  number of phionda to be integrated
  !  NK==0  we integrate Rho
  !
  !  NDERIV = 1 make the integration of the kinetic energy 
  ! 
  ! OUTPUT
  !    RTA = is the final integral value 
  ! USES
  ! USED BY
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_quanty
  implicit none
  integer(ip) :: NK,NDERIV 
  complex(rp) :: funcion(*)
  real(rp)    :: rta,promed,phi(mnode)
  real(rp)    :: gpvol(mgaus)                          
  integer(ip) :: ielem,igaus,inode,ipoin                          
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj


  rta=0.0_rp

  if( INOTMASTER ) then

     do ielem=1,nelem ! DO JEL=1,NE
        
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        plapl = 0
        porde = lorde(pelty)
        ptopo = ltopo(pelty)

        IF( nk == 0 ) THEN
           do inode = 1,pnode
              ipoin         =  lnods(inode,ielem)
              elcod(1,inode) = coord(1,ipoin)
              elcod(2,inode) = coord(2,ipoin)
              elcod(3,inode) = coord(3,ipoin)
              PHI(inode)     = real(funcion(ipoin))  ! rho(ipoin)
           end do
        else if( nk == 1 ) then
           do inode = 1,pnode
              ipoin          = lnods(inode,ielem)
              elcod(1,inode) = coord(1,ipoin)
              elcod(2,inode) = coord(2,ipoin)
              elcod(3,inode) = coord(3,ipoin)
              PHI(inode)     = CONJG(funcion(ipoin))*funcion(ipoin)
           end do
        else
           call runend('NK VALE 2')
        end if
        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        call elmcar(pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
             elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,gphes,ielem)

        if( nderiv == 0 ) then

           do igaus = 1,pgaus
              promed = 0.0_rp
              do inode = 1,pnode
                 promed = promed + PHI(inode) * elmar(pelty)%shape(inode,igaus)
              enddo
              rta = rta + promed * gpvol(igaus)
           end do

        else

           call runend('NOT GOOD')
           do igaus = 1,pgaus
              promed = 0.0_rp
              do inode =1,pnode
                 !promed = promed + PHI(inode)*(elmar(pelty)%deriv(1,inode,igaus)+&
                 !     elmar(pelty)%deriv(2,inode,igaus)+elmar(pelty)%deriv(3,inode,igaus))
                 promed = promed + PHI(inode) &
                      &   * ( elmar(pelty)%deriv(1,inode,igaus) &
                      &   +   elmar(pelty)%deriv(2,inode,igaus) &
                      &   +   elmar(pelty)%deriv(3,inode,igaus) )
              enddo
              rta    = rta + promed * promed * gpvol(igaus)
           end do

        endif


     enddo

  end if

  call pararr('SUM',0_ip,1_ip,rta)

end subroutine qua_3dinte


subroutine qua_3dinte_R(NK,NDERIV,funcion,rta)
  !-----------------------------------------------------------------------
  !****f* 3dinte
  ! NAME 
  !    3dinte
  ! DESCRIPTION
  !    This routine make an three-dimensional integral over the domain
  ! INPUT
  !   funcion (*) ?
  !  NK  number of phionda to be integrated
  !  NK==0  we integrate Rho
  !
  !  NDERIV = 1 make the integration of the kinetic energy 
  ! 
  ! OUTPUT
  !    RTA = is the final integral value 
  ! USES
  ! USED BY
  !-----------------------------------------------------------------------
  use def_master, only : INOTMASTER
  use def_domain
  use def_quanty
  implicit none
  integer(ip) ::  NK,NDERIV 
  REAL(rp)    :: funcion(*)
  real(rp)    :: rta,promed,phi(mnode)
  real(rp)    :: gpvol(mgaus)                          
  integer(ip) :: ielem,igaus,inode,ipoin                          
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj



  rta=0.0_rp

  if( INOTMASTER ) then

     do ielem = 1,nelem 

        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        plapl = 0
        porde = lorde(pelty)
        ptopo = ltopo(pelty)

        do inode = 1,pnode
           ipoin           = lnods(inode,ielem)
           elcod(1,inode)  = coord(1,ipoin)
           elcod(2,inode)  = coord(2,ipoin)
           elcod(3,inode)  = coord(3,ipoin)
           IF( NK == 0 ) THEN
              PHI(inode) = funcion(ipoin) 
           ELSE  
              PHI(inode) = funcion(ipoin) * funcion(ipoin)
           ENDIF
        enddo

        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        call elmcar(pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
             elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,gphes,ielem)


        if( nderiv == 0 ) then

           do igaus = 1,pgaus
              promed = 0.0_rp
              do inode = 1,pnode
                 promed = promed + PHI(inode)*elmar(pelty)%shape(inode,igaus)
              enddo
              rta = rta + promed * gpvol(igaus)
           enddo

        else

           do igaus=1,pgaus
              promed = 0.0_rp
              do inode = 1,pnode
                 promed = promed + PHI(inode)*(elmar(pelty)%deriv(1,inode,igaus)+&
                      elmar(pelty)%deriv(2,inode,igaus)+elmar(pelty)%deriv(3,inode,igaus))
              enddo
              rta = rta + promed * promed * gpvol(igaus)
           enddo

        endif

     enddo

  end if

  call pararr('SUM',0_ip,1_ip,rta)

end subroutine qua_3dinte_R

subroutine integracuad(l,M,rta) 
  !-----------------------------------------------------------------------
  !****f* 3dinte
  ! NAME 
  !    3dinte
  ! DESCRIPTION
  !    This routine make an three-dimensional integral over the domain
  ! INPUT
  !   funcion (*) ?
  !  NK  number of phionda to be integrated
  !  NK==0  we integrate Rho
  !
  !  NDERIV = 1 make the integration of the kinetic energy 
  ! 
  ! OUTPUT
  !    RTA = is the final integral value 
  ! USES
  ! USED BY
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_quanty
  implicit none
  integer(ip) :: l,m
  complex(rp) :: rta,phi(mnode),armonico,promed(mgaus)
  real(rp)    :: gpvol(mgaus)                          
  integer(ip) :: ielem,igaus,inode,ipoin                          
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo,n_localn(mnode)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj

  rta=0.0_rp

  if( INOTMASTER ) then

     do ielem=1,nelem ! DO JEL=1,NE


        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        plapl = 0
        porde = lorde(pelty)
        ptopo = ltopo(pelty)

        do igaus=1,pgaus
           promed(igaus)= 0.0_rp
        enddo

        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           n_localn(inode)=ipoin
           elcod(1,inode) = coord(1,ipoin)
           elcod(2,inode) = coord(2,ipoin)
           elcod(3,inode) = coord(3,ipoin)

           if(l.eq.0) then
              PHI(inode)= RHOON(ipoin,1)*armonico(elcod(1,inode),elcod(2,inode),elcod(3,inode),l,m)
           elseif(l.eq.1) then
              PHI(inode)= RHOON(ipoin,1)*armonico(elcod(1,inode),elcod(2,inode),elcod(3,inode),l,m)* &
                   &               sqrt(elcod(1,inode)**2+elcod(2,inode)**2+elcod(3,inode)**2)
           else 
              PHI(inode)= RHOON(ipoin,1)*armonico(elcod(1,inode),elcod(2,inode),elcod(3,inode),l,m)* &
                   &               (elcod(1,inode)**2+elcod(2,inode)**2+elcod(3,inode)**2)

           endif
        enddo

        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        call elmcar(pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
             elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,gphes,ielem)

        do igaus=1,pgaus
           do inode=1,pnode
              promed(igaus)= promed(igaus) + PHI(inode)*elmar(pelty)%shape(inode,igaus)
           enddo
           rta = rta + promed(igaus)*gpvol(igaus)
        enddo

     enddo

  end if

  call pararr('SUM',0_ip,1_ip,rta)

end subroutine integracuad
