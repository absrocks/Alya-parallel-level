subroutine qua_exaerr(itask)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_exaerr
  ! NAME 
  !    qua_exaerr   !!! OJO name of variables
  ! DESCRIPTION
  !    This routine computes the errors with respect to the exact solution
  ! USES
  !   qua_exacso
  !   elmder
  ! USED BY
  !    qua_output
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  use mod_memchk
  use mod_postpr
  implicit none 
  integer(ip), intent(in) :: itask
  integer(ip) :: ielem,inode,ipoin,igaus,idime,pnode
  integer(ip) :: pelty,pmate
  real(rp)    :: cartd(ndime,mnode),xjaci(9),xjacm(9) 
  real(rp)    :: elcod(ndime,mnode),elqua(mnode)       
  real(rp)    :: diqua,digrt,abqua,abgrt,gpcod(3)
  real(rp)    :: gpqua,gpgrt(3),exqua,exgrt(3)
  real(rp)    :: dvolu,detjm,dummr(50)

  select case(itask)

  case(1_ip)
     !
     ! Error calculations
     !
     err01_qua = 0.0_rp 
     err02_qua = 0.0_rp

     elements_2: do ielem=1,nelem

        pelty=ltype(ielem)
        pnode=nnode(pelty)
        pmate=1
        if(nmate>1) pmate=lmate(ielem)
        !
        ! Gather operations
        !
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do idime=1,ndime
              elcod(idime,inode)=coord(idime,ipoin)
           end do
           elqua(inode)=phion(ipoin,1)
        end do

        gauss_points_2: do igaus=1,ngaus(pelty)
           !
           ! Cartesian derivatives and Jacobian
           !
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                elcod,cartd,detjm,xjacm,xjaci)
           dvolu=elmar(pelty)%weigp(igaus)*detjm  
           !
           ! Gauss point values
           !
           gpqua=0.0_rp
           do inode=1,pnode
              gpqua=gpqua+elqua(inode)*elmar(pelty)%shape(inode,igaus)
           end do
           do idime=1,ndime
              gpgrt(idime)=0.0_rp
              gpcod(idime)=0.0_rp
           end do
           do inode=1,pnode
              do idime=1,ndime
                 gpgrt(idime)=gpgrt(idime)&
                      +elqua(inode)*cartd(idime,inode)
                 gpcod(idime)=gpcod(idime)&
                      +elcod(idime,inode)*elmar(pelty)%shape(inode,igaus)
              end do
           end do

           if(kfl_naxis==1) then
              dvolu=dvolu*gpcod(1)*twopi
           end if

           call qua_exacso(&
                1_ip,gpcod,dummr,dummr,&
                dummr,dummr,dummr,dummr,exqua,exgrt,dummr)

           diqua        = dabs(gpqua-exqua)
           abqua        = dabs(exqua)
           err01_qua(1) = err01_qua(1) + diqua*dvolu
           err02_qua(1) = err02_qua(1) + diqua*diqua*dvolu
!           err0i_qua(1) = max(err0i_qua(1),diqua)
           err01_qua(2) = err01_qua(2) + abqua*dvolu
           err02_qua(2) = err02_qua(2) + exqua*exqua*dvolu
!           err0i_qua(2) = max(err0i_qua(2),abqua)

           do idime=1,ndime
              digrt        = dabs(gpgrt(idime)-exgrt(idime))
              abgrt        = dabs(exgrt(idime))
!              err11_qua(1) = err01_qua(1) + digrt*dvolu
!              err12_qua(1) = err02_qua(1) + digrt*digrt*dvolu
!              err1i_qua(1) = max(err0i_qua(1),digrt)
!              err11_qua(2) = err01_qua(2) + abgrt*dvolu
!              err12_qua(2) = err02_qua(2) + exgrt(idime)*exgrt(idime)*dvolu
!              err1i_qua(2) = max(err0i_qua(2),abgrt)
           end do

        end do gauss_points_2

     end do elements_2

     if(err01_qua(2)>zequa) err01_qua(1)=err01_qua(1)/err01_qua(2)  
     if(err02_qua(2)>zequa) err02_qua(1)=err02_qua(1)/dsqrt(err02_qua(2))
!     if(err0i_qua(2)>zequa) err0i_qua(1)=err0i_qua(1)/err0i_qua(2)  
!     if(err11_qua(2)>zequa) err11_qua(1)=err11_qua(1)/err11_qua(2)  
!     if(err12_qua(2)>zequa) err12_qua(1)=err12_qua(1)/dsqrt(err12_qua(2))
!     if(err1i_qua(2)>zequa) err1i_qua(1)=err1i_qua(1)/err1i_qua(2)       
     
	 write(lun_outpu_qua,100) ittim,cutim,err01_qua(1),err02_qua(1)


!     write(lun_outpu_qua,100) ittim,cutim,&
!          err01_qua(1),err02_qua(1),err0i_qua(1),&
!          err11_qua(1),err12_qua(1),err1i_qua(1)

  case(2_ip)
     !
     ! Compute error
     !
     do ipoin=1,npoin
        call qua_exacso(&
             1_ip,coord(1,ipoin),dummr,dummr,&
             dummr,dummr,dummr,dummr,exqua,dummr,dummr)
        gesca(ipoin)=phion(ipoin,1)-exqua
     end do

  end select

100 format(///,&
       & 5x,'FINITE ELEMENT ERRORS',/,&
       & 5x,'=====================',//,&
       & 5x,'TIME STEP= ',i5,', CURRENT TIME= ',e16.8e3,//,&
       & '     NORM       VALUE  ',/,5x,23('-'),/,&
       & '     W(0,1) ',e16.8e3,/,&
       & '     W(0,2) ',e16.8e3,/,&
       & '     W(0,i) ',e16.8e3,/,&
       & '     W(1,1) ',e16.8e3,/,&
       & '     W(1,2) ',e16.8e3,/,&
       & '     W(1,i) ',e16.8e3,/,5x,23('-'))

end subroutine qua_exaerr
