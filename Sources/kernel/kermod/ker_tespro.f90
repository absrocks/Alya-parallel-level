subroutine ker_tespro()
  !-----------------------------------------------------------------------
  !****f* Kermod/ker_tespro
  ! NAME 
  !    ker_tespro
  ! DESCRIPTION
  ! USED BY
  !    ker_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame 
  use def_master
  use def_domain
  use mod_ker_proper
  use mod_gradie
  implicit none
  integer(ip) :: ielem,igaus,pgaus,dummi,idime,kmate(nmate)
  integer(ip) :: plapl,pnode,inode,ipoin,pelty,imate
  real(rp)    :: gpden(mgaus)
  real(rp)    :: grden(ndime,mgaus)
  real(rp)    :: gpsha(mnode,mgaus)                  
  real(rp)    :: gpder(ndime,mnode,mgaus)                  
  real(rp)    :: gpcar(ndime,mnode,mgaus)            
  real(rp)    :: gphes(ntens,mnode,mgaus)            
  real(rp)    :: gpvol(mgaus)                        
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: gpcor(ndime,mgaus)
  real(rp)    :: gpcog(ndime)
  real(rp)    :: dummr(mgaus)

  if( IMASTER ) return

  !call memgen(0_ip,npoin,0_ip)
  !call memgen(0_ip,ndime,npoin)
  !do ipoin = 1,npoin
  !   gesca(ipoin) = 2.0_rp * coord(1,ipoin) - 3.0_rp * coord(2,ipoin)
  !end do
  !call grasca(gesca,gevec)
  !do ipoin = 1,npoin
  !   print*,'GRDEN=',gevec(1,ipoin), 2.0_rp
  !   print*,'GRDEN=',gevec(2,ipoin),-3.0_rp
  !end do
  !call memgen(2_ip,ndime,npoin)
  !call memgen(2_ip,npoin,0_ip)
  !stop

  kmate = 0
  do ielem = 1,nelem

     pelty = ltype(ielem)
     pgaus = ngaus(pelty)
     pnode = nnode(pelty)
     imate = lmate(ielem)
     plapl = 0

     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        do idime = 1,ndime
           elcod(idime,inode) = coord(idime,ipoin)
        end do
     end do
     gpcor = 0.0_rp
     do igaus = 1,pgaus
        do inode = 1,pnode
           do idime = 1,ndime
              gpcor(idime,igaus) = gpcor(idime,igaus) + elcod(idime,inode) * elmar(pelty) % shape(inode,igaus)
           end do
        end do
     end do
     gpcog = 0.0_rp
     do inode = 1,pnode
        do idime = 1,ndime
           gpcog(idime) = gpcog(idime) + elcod(idime,inode)
        end do
     end do
     do idime = 1,ndime
        gpcog(idime) = gpcog(idime) / real(pnode,rp)
     end do
     
     call elmca2(&
          pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
          elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
          gpder,gpcar,gphes,ielem) 
     !
     ! PGAUS
     !
     call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,gpsha,gpcar)
     call ker_proper('GRADI','PGAUS',dummi,ielem,grden,pnode,pgaus,gpsha,gpcar)
     print*,' '
     print*,'PGAUS, MATERIAL= ',imate
     print*,' '
     do igaus = 1,pgaus
        print*,'GPDEN=',gpden(igaus),2.0_rp*gpcor(1,igaus)-3.0_rp*gpcor(2,igaus)
     end do
     do igaus = 1,pgaus
        print*,'GRADI=',grden(1,igaus), 2.0_rp
        print*,'GRADI=',grden(2,igaus),-3.0_rp
     end do
     !
     ! IGAUS
     !
     igaus = 3
     call ker_proper('DENSI','IGAUS',igaus,ielem,gpden(igaus:), pnode,pgaus,gpsha,gpcar)
     call ker_proper('GRADI','IGAUS',igaus,ielem,grden(:,igaus),pnode,pgaus,gpsha,gpcar)
     print*,' '
     print*,'IGAUS, MATERIAL= ',imate
     print*,' '
     print*,'GPDEN=',gpden(igaus),2.0_rp*gpcor(1,igaus)-3.0_rp*gpcor(2,igaus)
     print*,'GRADI=',grden(1,igaus),2.0_rp
     print*,'GRADI=',grden(2,igaus),-3.0_rp
     !
     ! COG
     !
     call elmcar(&
          pnode,1_ip,plapl,elmar(pelty) % weicg,elmar(pelty) % shacg,&
          elmar(pelty) % dercg,dummr,elcod,gpvol,gpcar,&
          gphes,ielem)
     igaus = 1
     call ker_proper('DENSI','COG  ',1_ip,ielem,gpden)
     call ker_proper('GRADI','COG  ',1_ip,ielem,grden,pnode,pgaus,gpsha,gpcar)
     print*,' '
     print*,'C.O.G, MATERIAL= ',imate
     print*,' '
     print*,'GPDEN=',gpden(igaus),2.0_rp*gpcog(1)-3.0_rp*gpcog(2)
     print*,'GRADI=',grden(1,igaus), 2.0_rp
     print*,'GRADI=',grden(2,igaus),-3.0_rp     
     !
     ! PNODE
     !
     call ker_proper('DENSI','PNODE',1_ip,ielem,gpden)
     call ker_proper('GRDEN','PNODE',1_ip,ielem,grden)
     print*,' '
     print*,'PNODE, MATERIAL= ',imate
     print*,' '
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        print*,'GPDEN=',gpden(inode),2.0_rp*coord(1,ipoin)-3.0_rp*coord(2,ipoin)
     end do
     do inode = 1,pnode
        print*,'GRADI=',grden(1,inode), 2.0_rp
        print*,'GRADI=',grden(2,inode),-3.0_rp
     end do

  end do
  stop
  !
  ! NPOIN
  !
  print*,' '
  print*,'NPOIN'
  print*,' '
  call memgen(0_ip,npoin,zero)
  call ker_proper('DENSI','NPOIN',dummi,dummi,gesca)
  do ipoin = 1,npoin
     print*,'DENSI=',gesca(ipoin),2.0_rp*coord(1,ipoin)-3.0_rp*coord(2,ipoin)
  end do
  call memgen(2_ip,npoin,zero)

  print*,' '
  call memgen(0_ip,ndime,npoin)
  call ker_proper('GRADI','NPOIN',dummi,dummi,gevec)
  do ipoin = 1,npoin
     print*,'GRADI=',gevec(1,ipoin), 2.0_rp
     print*,'GRADI=',gevec(2,ipoin),-3.0_rp
  end do
  call memgen(2_ip,ndime,npoin)

end subroutine ker_tespro
