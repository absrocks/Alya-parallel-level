!subroutine hlm_source()
  !------------------------------------------------------------------------
  !****f* Helmoz/hlm_source
  ! NAME 
  !    hlm_source
  ! DESCRIPTION
  !    ORDER=1:
  !      Helmozature equation, boundary operations
  ! USES
  ! USED BY
  !    hlm_matrix 
  !***
  !------------------------------------------------------------------------
  !use def_parame
  !use def_master
  !use def_domain
  !use def_helmoz
  !implicit none
  !complex(rp) :: elrhs(mnode*nequs_hlm)
  !complex(rp) :: gpmas(mgaus)
  !complex(rp) :: gpmav(ndime,mgaus)
  !complex(rp) :: gpelv(ndime,mgaus)
  !complex(rp) :: gpels(mgaus)
  !complex(rp) :: dummx
  !real(rp)    :: elcod(ndime,mnode)
  !integer(ip) :: ielem,inode,ipoin,dummi,idime,isour,kelem
  !integer(ip) :: igaus,jnode,jelty,jelem,pnode,pelty,pgaus
  !real(rp)    :: gpvol(mgaus),xjacm(9),gpcod(3)
  !real(rp)    :: deriv(ndime*mnode),coloc(3),shapf(mnode)
  !
  ! Loop over elements  
  !
  !sources: do isour = 1,nsour_hlm

     !kelem = sourc_hlm(isour) % nelem
	
     !elements: do ielem = 1,kelem

        !pelty = sourc_hlm(isour) % ltype(ielem)
        !pnode = nnode(pelty)
        !pgaus = ngaus(pelty)
        !
        ! Inititalize
        !
        !do inode = 1,mnode*nequs_hlm
           !elrhs(inode) = CMPLX( 0.0_rp , 0.0_rp )
        !end do
        !
        ! Gather operations
        !
        !do inode = 1,pnode
           !ipoin = sourc_hlm(isour) % lnods(inode,ielem)
           !do idime = 1,ndime
              !elcod(idime,inode) = sourc_hlm(isour) % coord(idime,ipoin)
           !end do
        !end do
        !
        ! GPVOL
        !
        !if( sourc_hlm(isour) % ndime == ndime ) then
           !do igaus = 1,pgaus         
              !call jacdet(&
                   !ndime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                   !xjacm,gpvol(igaus))
           !end do
        !else
           !do igaus = 1,pgaus 
	      !call bouder(&
                   !pnode,sourc_hlm(isour) % ndime,ndimb,elmar(pelty)%deriv(1,1,igaus),&
                   !elcod,xjacm,gpvol(igaus))
           !end do
        !end if
        !
        ! Gauss point values
        !

        !if( kfl_model_hlm == 1 .or. kfl_model_hlm == 3 ) then
	      !do igaus = 1,pgaus 
		 !gpmas(igaus) = CMPLX( 0.0_rp , 0.0_rp )
              !do idime = 1,ndime
                 !gpmav(idime,igaus) = CMPLX( 0.0_rp , 0.0_rp )
              !end do
              !do inode = 1,pnode
		 !ipoin = sourc_hlm(isour) % lnods(inode,ielem)
                 !gpmas(igaus) = gpmas(igaus)                 &
                      !+ elmar(pelty) % shape(inode,igaus)    &
                      !* sourc_hlm(isour) % magsp(ipoin)
		 !do idime = 1,ndime
                    !gpmav(idime,igaus) = gpmav(idime,igaus)  &
                         !+ elmar(pelty) % shape(inode,igaus) &
                         !* sourc_hlm(isour) % magvp(idime,ipoin)
                 !end do
              !end do
           !end do
        !end if

        !if( kfl_model_hlm == 2 .or. kfl_model_hlm == 3 ) then
           !do igaus = 1,pgaus  
              !gpels(igaus) = CMPLX( 0.0_rp , 0.0_rp )
              !do idime = 1,ndime
                 !gpelv(idime,igaus) = CMPLX( 0.0_rp , 0.0_rp )
              !end do
              !do inode = 1,pnode
                 !ipoin = sourc_hlm(isour) % lnods(inode,ielem)
                 !gpels(igaus) = gpels(igaus)                 &
                      !+ elmar(pelty) % shape(inode,igaus)    &
                      !* sourc_hlm(isour) % elesp(ipoin)
                 !do idime = 1,ndime
                    !gpelv(idime,igaus) = gpelv(idime,igaus)  &
                         !+ elmar(pelty) % shape(inode,igaus) &
                         !* sourc_hlm(isour) % elevp(idime,ipoin)
                 !end do
              !end do
           !end do
        !end if
        !
        ! Assembly
        !
        !do igaus = 1,pgaus
           !gpvol(igaus) = elmar(pelty)%weigp(igaus)*gpvol(igaus)
           !do idime = 1,ndime
              !gpcod(idime) = 0.0_rp
              !do inode = 1,pnode
                 !gpcod(idime) = gpcod(idime) + elmar(pelty) % shape(inode,igaus) &
                      !* elcod(idime,inode)
              !end do
           !end do
           !
           ! Find host element
           !
           !call elsest(&
                !2_ip,1_ip,ielse,mnode,ndime,npoin,nelem,nnode,&
                !lnods,ltype,ltopo,coord,gpcod,relse,&
                !jelem,shapf,deriv,coloc,dummi)
	   !if ( ISEQUEN ) then
           	!if( jelem == 0 ) call runend('COULD NOT FIND ELEMENT FOR SOURCE GAUSS POINT')
	   !else
	   	!if( (jelem == 0) .or. (jelem == -1) ) return
 	   !end if
           !jelty = ltype(jelem)
           !jnode = nnode(jelty)
           !do inode=1,jnode
	   	!if (shapf(inode)<1.0E-10) shapf(inode)= 0.0_rp
	   !end do
           !
           ! Element source contribution
           !
           !call hlm_assour(&
                !jnode,gpmas(igaus),gpmav(1,igaus),&
                !gpels(igaus),gpelv(1,igaus),&
                !gpvol(igaus),shapf,elrhs)
           !
           ! Assembly
           !
           !call hlm_assemb(2_ip,jnode,lnods(1,jelem),dummx,elrhs)

        !end do

     !end do elements

  !end do sources

!end subroutine hlm_source
