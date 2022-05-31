
subroutine geometrical_tests()
  use def_kintyp
  use def_master
  use def_domain
  use def_elmtyp
  use mod_parall
  use def_kermod
  use mod_communications
  use mod_graphs
  use mod_elmgeo
  use mod_elsest
  use mod_postpr
  use mod_partitioning
  use mod_space_filling_curve
  use mod_unity_tests
  use mod_maths
  implicit none
  real(rp),    pointer :: xx(:)
  integer(ip)          :: pelty,ipoin,inotfound,pnode,kelem,idime,inode,np
  integer(ip)          :: ii,ifoun,ielem,ielpo,jelem,jj,kk,ll,ww,d,ipart,wt
  real(rp)             :: elcod(3,6),coglo(3),bocod(1:3,4),nn(3),yy(3),comin(3),comax(3)
  real(rp)             :: shapf_test(64),shapf(64)
  real(rp)             :: deriv(300),coloc(3),dista,toler,pvect(3)
  integer(ip)          :: n,PAR_COMM_TO_USE,boxes(3)
  integer(ip), pointer :: weight(:),lpart(:),lweight(:),lpart_cpy(:),binel(:)
  real(rp),    pointer :: xyz(:,:)
  integer(4)           :: COMM4

  return
  goto 10

  call unity_tests_polynomial_integration()
  call runend('O.K.!')


  call space_filling_curve_output(2_ip,10_ip,1.0_rp,1.0_rp,1.0_rp)

  call runend('O.K.!')
goto 10
  n = 4


  !if( inotmaster ) then
  !   allocate(lpart(npoin))
  !   allocate(lweight(npoin))
  !   allocate(xyz(ndime,npoin))
  !   lweight    = 1
  !   !lweight(1:50) = 4
  !   lpart      = 0
  !   pvect      = (/ 0.5_rp,0.5_rp,0.0_rp /)
  !   boxes(1:2) = 10
  !   call partitioning_oriented_bin(n,npoin,boxes,pvect,coord,lpart,lweight)
  !   return
  !end if

  if( inotmaster ) then
     call PAR_DEFINE_COMMUNICATOR('IN MY CODE WITHOUT MASTER',COMM4)

     allocate(lpart(nelem))
     allocate(lweight(nelem))
     allocate(xyz(ndime,nelem))
     lweight    = 1
     !lweight(1:50) = 4
     lpart      = 0
     pvect      = (/ 1.0_rp,0.0_rp,0.0_rp /)
     boxes(1:3) = 200
     do ielem = 1,nelem
        do idime = 1,ndime
           xyz(idime,ielem) = sum(coord(idime,lnods(1:lnnod(ielem),ielem))) / real(lnnod(ielem),rp)
        end do
     end do
     !call partitioning_oriented_bin(ndime,n,nelem,boxes,pvect,xyz,lpart,'LAYER',lweight,COMM4)
     call partitioning_oriented_bin(ndime,n,nelem,boxes,pvect,xyz,lpart,'LOAD BALANCE',lweight,COMM4)
  end if
  do ii = 1,n
     if( inotmaster ) then
        kelem = 0
        do ielem = 1,nelem
           if( lpart(ielem) == ii ) then
              kelem = kelem + lweight(ielem)
           end if
        end do
     end if
     call PAR_SUM(kelem)
     if(inotslave) print*,'Part n=',ii,kelem
  end do

  call postpr_right_now('LPART','SCALA','NELEM',lpart)

  call runend('O.K.!')

  if( inotmaster ) then
     call PAR_DEFINE_COMMUNICATOR('IN MY CODE WITHOUT MASTER',COMM4)
     allocate(lpart(npoin))
     allocate(lpart_cpy(npoin))
     allocate(lweight(npoin))

     lweight = 1
     !lweight(1:50) = 4
     lpart = 0
     !lpart(1:10) = -1
     call partitioning_frontal(n,npoin,r_dom,c_dom,lpart,lweight,COMM4,'NODE')
     !call partitioning_frontal(n,npoin,r_dom,c_dom,lpart,lweight)
     !call partitioning_metis(n,npoin,r_dom,c_dom,lpart,lweight,COMM4,'NODE')
     lpart_cpy = lpart
     call PAR_INTERFACE_NODE_EXCHANGE(lpart_cpy,'MAX')
     do ii = 1,npoin
        if( lpart_cpy(ii) /= lpart(ii) ) call runend('ERROR')
     end do
  end if

  ielem = 0
  do ii = 1,n
     if( inotmaster ) then
        kelem = 0
        do ipoin = 1,npoin
           if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
              if( lpart(ipoin) == ii ) then
                 kelem = kelem + lweight(ipoin)
                 ielem = ielem + lweight(ipoin)
              end if
           end if
        end do
     end if
     call PAR_SUM(kelem)
     if(imaster) print*,'Part n=',ii,kelem
  end do
  call PAR_SUM(ielem)
  if(imaster) print*,'Total n=',ielem

     !call partitioning_metis(n,npoin,r_dom,c_dom,lpart)
  call postpr_right_now('LPART','SCALA','NPOIN',lpart)
  call runend('O.K.!')

10 continue

  !----------------------------------------------------------------------
  !
  ! Mesh partitioning into np subdomains using a (n x n x n) bin
  !
  !----------------------------------------------------------------------

  n  = 16
  np = 4

  allocate(weight(n**ndime))
  allocate(binel(nelem))
  weight = 0
  wt     = 0
  boxes  = n
  do idime = 1,ndime
     comin(idime) = minval(coord(idime,:))
     comax(idime) = maxval(coord(idime,:))
  end do

  do ielem = 1,nelem

     yy = 0.0_rp
     do inode = 1,lnnod(ielem)
        ipoin = lnods(inode,ielem)
        yy(1:ndime) = yy(1:ndime) + coord(1:ndime,ipoin)
     end do
     yy = yy / real(lnnod(ielem),rp)
     call maths_mapping_coord_to_3d(ndime,boxes,comin,comax,yy,ii,jj,kk)

     ll           = maths_mapping_3d_to_1d(n,n,0_ip,ii,jj,kk)
     binel(ielem) = ll
     weight(ll)   = weight(ll) + lgaus(ielem)
     wt           = wt + lgaus(ielem)
  end do

  call space_filling_curve_output(ndime,n,1.0_rp/real(n,rp),1.0_rp/real(n,rp),1.0_rp/real(n,rp),weight,'LEXICAL')

  wt    = wt / np
  ww    = 0
  kk    = 1
  ipart = 1

  do d = 1,n**ndime
     if( ndime == 2 ) then
        call maths_sfc_1d_to2d3d_tab(n,d,ii,jj)
        ll = maths_mapping_3d_to_1d(n,n,0_ip,ii,jj,kk)
     else
        call maths_sfc_1d_to2d3d_tab(n,d,ii,jj,kk)
        ll = maths_mapping_3d_to_1d(n,n,n,ii,jj,kk)
     end if
     ww = ww + weight(ll)
     do ielem = 1,nelem
        if( binel(ielem) == ll ) then
           binel(ielem) = -ipart
        end if
     end do
     if( ww >= wt ) then
        ww    = 0
        ipart = ipart + 1
     end if
  end do
  binel = abs(binel)
  call postpr_right_now('LPART','SCALA','NELEM',binel)
  deallocate(weight)
  deallocate(binel)

  call runend('O.K.!')

  if( INOTMASTER ) then
     kfl_lele2 = 1
     kfl_lelp2 = 1
     call lelpo2()
     !call graphs_eleele(&
     !  nelem,npoin,mnode,mepoi,lnods,lnnod,&
     !  pelpo,lelpo,nedge,medge,pelel,lelel)

     do ielem = 1,nelem
        if( leinv_loc(ielem) == 60 ) then
           print*,'allo'
           do ielpo = pelel_2(ielem),pelel_2(ielem+1)-1
              jelem = lelel_2(ielpo)
              if( jelem > nelem ) then
                 print*, commd%neights(leldo(1,jelem-nelem)),leinv_loc(jelem)
              else
                 print*,kfl_paral,leinv_loc(jelem)
              end if
           end do
        end if
     end do
  end if
  call runend('O.K.!')

  elcod(1:3,1) = (/ 5.190989429460000E-002_rp,  5.559730980170000E-002_rp,  1.699681829750000E-002_rp /)
  elcod(1:3,2) = (/ 5.218503067920000E-002_rp,  5.571069138870000E-002_rp,  1.677076601660000E-002_rp /)
  elcod(1:3,3) = (/ 5.223954143100000E-002_rp,  5.566297143330001E-002_rp,  1.709512581440000E-002_rp /)
  elcod(1:3,4) = (/ 5.191602653030000E-002_rp,  5.555485050010001E-002_rp,  1.698885972610000E-002_rp /)
  elcod(1:3,5) = (/ 5.219615120440001E-002_rp,  5.566891978310001E-002_rp,  1.676463817900000E-002_rp /)
  elcod(1:3,6) = (/ 5.225132475050000E-002_rp,  5.562191528080000E-002_rp,  1.708630990100000E-002_rp /)

  coglo(1:3)   = (/ 5.200177051836418E-002_rp,  5.560847058816261E-002_rp,  1.702584393891663E-002_rp /)

  write(*,*)
  do ii = 1,13
     write(*,*) elmar(TRI03) % posgp(1:2,ii)
  end do
  stop
  print*,'a'

  pnode = 6
  pelty = 34

  elcod(1:3,1) = (/ -1.510463160000000E-002_rp,  0.259003415200000_rp,      -0.279569664200000_rp /)
  elcod(1:3,2) = (/ -1.528435380000000E-002_rp,  0.258834177700000_rp,      -0.279176316300000_rp /)
  elcod(1:3,3) = (/ -1.521732580000000E-002_rp,  0.259066797900000_rp,      -0.279486836700000_rp /)
  elcod(1:3,4) = (/ -1.518098590000000E-002_rp,  0.258992770000000_rp,      -0.279690490800000_rp /)
  elcod(1:3,5) = (/ -1.535433820000000E-002_rp,  0.258811283200000_rp,      -0.279303922500000_rp /)
  elcod(1:3,6) = (/ -1.524770640000000E-002_rp,  0.259007622500000_rp,      -0.279623231500000_rp /)
  coglo(1:3)   = (/ -1.523013400000000E-002_rp,  0.258855219600000_rp,      -0.279673926200000_rp /)

  print*,'b'
  if( elmgeo_inside_element_using_faces(pelty,elcod,coglo) ) then
     print*,'c=oui'
     call elmgeo_newrap(coglo,coloc,ndime,pnode,elcod,shapf_test,deriv)
  else
     print*,'c=non'
  end if

  print*,'d'
  !coglo(1:3) = (/ -600294.815592175_rp,       -830794.824829704_rp,  -300862.510440974_rp /)
  !bocod(1:3,1) = (/ 8.131726089999999E-002_rp ,  0.218109726900000_rp ,      -0.239409774500000_rp   /)
  !bocod(1:3,2) = (/ 8.112300190000001E-002_rp ,  0.217919886100000_rp ,      -0.239590942900000_rp   /)
  !bocod(1:3,3) = (/ 8.108361970000000E-002_rp ,  0.217949596700000_rp ,      -0.239581189900000_rp   /)
  !bocod(1:3,4) = (/ 8.127600620000000E-002_rp ,  0.218136846100000_rp ,      -0.239400959800000_rp   /)
  !call elmgeo_nearest_point_on_QUA04(bocod,coglo,yy,nn)
  !call elmgeo_natural_coordinates_on_QUA04(bocod,yy,coloc)
  !print*,'coloc=',coloc
  !stop
  !COGLO
  !TOLER
  ! 3.310787378353252E-004
  ! 1.164153218269348E-010

  stop

  inotfound = 0
  do ipoin = 1,npoin

     coglo(1:ndime) = coord(1:ndime,ipoin)
  !do ielem = 1,nelem
  !do ielem = 1,1

     !ipoin = ielem
     !pnode = lnnod(ielem)
     !do idime = 1,ndime
     !   coglo(idime) = sum(coord(idime,lnods(1:pnode,ielem)))
     !end do
     !coglo = coglo / real(pnode,rp)
     !coglo(1:3) = (/ 0.17960793999999999_rp,0.26761568000000002_rp,-0.22290830000000000_rp /)
     !coglo = (/ 20.0_rp,20.0_rp,20.0_rp /)
     !ielse(15) = 1

     call elsest_host_element(&
          ielse,relse,1_ip,meshe(ndivi),coglo,kelem,&
          shapf_test,deriv,coloc,dista)
     if(kelem==0) then
        print*,'not a=',ipoin
        print*,'not b=',coglo(1:ndime)
        inotfound = inotfound + 1
        !stop
     else
        !yy = 0.0_rp
        !do inode = 1,lnnod(kelem)
        !   yy(1:ndime) = yy(1:ndime) + shapf_test(inode) * coord(1:ndime,lnods(inode,kelem))
        !end do
        !print*,'intersection=',yy(1:ndime)
     end if
  end do
  print*,'npoin     = ',npoin
  print*,'not found = ',inotfound
  call runend('O.K.!')

  bocod(1:3,1) = (/ 0.0_rp,0.0_rp,1.0_rp /)
  bocod(1:3,2) = (/ 0.2_rp,0.0_rp,1.0_rp /)
  bocod(1:3,3) = (/ 0.2_rp,2.0_rp,1.0_rp /)
  bocod(1:3,4) = (/ 0.0_rp,0.1_rp,1.0_rp /)

  coglo(1:3)   = (/ 0.75_rp,0.5_rp,-3.0_rp /)
  call elmgeo_nearest_point_on_QUA04(bocod,coglo,yy,nn)
  print*,'projection=',yy(1:3)
  print*,'normal=',nn(1:3)

stop
  ielem = 1
  elcod(1:ndime,1:6) = coord(1:ndime,lnods(1:6,ielem))
  coglo(1) = 0.109456842300000_rp
  coglo(2) = 0.262726287900000_rp
  coglo(3) =  -0.289190207300000_rp
  pelty = 34
print*,'caca'
  print*,'result=',  elmgeo_inside_element_using_faces(pelty,elcod,coglo,1.0e-6_rp)
print*,'pipi'

  stop
  call graphs_output(npoin,r_dom,c_dom,coord,90_4,'graphs-'//trim(intost(kfl_paral)))
  return

  !print*,'a=',leinv_loc(1:nelem_2)
  print*,'a=',kfl_paral,commd % neights(leldo(1,1:2))
  return
!size(leldo,2,KIND=ip)         ! Which neighbor  this element
!call runend('O.K.!')
  !leldo(2,jelem-nelem)         ! Local numbering of JELEM in my neighbor

  !r_elm_2 => meshe(ndivi) % r_elm_2
  !c_elm_2 => meshe(ndivi) % c_elm_2
  !do ielem = 1,meshe(ndivi) % nelem
  !   write(900+kfl_paral,'(i2,a,10(1x,i2))') leinv_loc(ielem),': ',leinv_loc(c_elm_2(r_elm_2(ielem):r_elm_2(ielem+1)-1))
  !end do
  !return

  if( inotmaster ) then
     nullify(xx)
     allocate(xx(nelem_2))
     xx = 0.0_rp
     if(kfl_paral==2 ) xx(1:nelem) = kfl_paral
     !call PAR_FROM_GHOST_ELEMENT_EXCHANGE(xx,'SUM','IN MY CODE')
     call PAR_GHOST_ELEMENT_EXCHANGE(xx,'SUM','IN MY CODE')
     if( kfl_paral == 1 ) then
        print*,'a=',xx(1:nelem)
        print*,'b=',xx(nelem+1:nelem_2)
     end if
  end if

end subroutine geometrical_tests

