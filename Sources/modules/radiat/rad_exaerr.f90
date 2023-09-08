subroutine rad_exaerr(itask)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_exaerr
  ! NAME 
  !    rad_exaer
  ! DESCRIPTION
  !    This routine computes the errors with respect to the exact solution
  ! USES
  ! USED BY
  !    rad_output
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  use mod_memchk
  use def_elmtyp
  implicit none 
  integer(ip), intent(in) :: itask
  integer(ip) :: ielem,inode,ipoin,igaus,idime,pnode
  integer(ip) :: pelty,pmate
  real(rp)    :: cartd(ndime,mnode),xjaci(9),xjacm(9) 
  real(rp)    :: elcod(ndime,mnode),eltem(mnode)       
  real(rp)    :: ditem,digrt,abtem,abgrt,gpcod(3)
  real(rp)    :: gptem,gpgrt(3),extem,exgrt(3)
  real(rp)    :: dvolu,detjm,dummr(50)

  select case(itask)

  case(1_ip)
     !
     ! Error calculations
     !
     err01_rad = 0.0_rp 
     err02_rad = 0.0_rp
     err0i_rad = 0.0_rp
     err11_rad = 0.0_rp 
     err12_rad = 0.0_rp
     err1i_rad = 0.0_rp

     elements_2: do ielem=1,nelem

        pelty=ltype(ielem)
        if(lelch(ielem)==ELFEM .and. pelty>0)then
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
              eltem(inode)=radav_rad(ipoin,1)
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
              gptem=0.0_rp
              do inode=1,pnode
                 gptem=gptem+eltem(inode)*elmar(pelty)%shape(inode,igaus)
              end do
              do idime=1,ndime
                 gpgrt(idime)=0.0_rp
                 gpcod(idime)=0.0_rp
              end do
              do inode=1,pnode
                 do idime=1,ndime
                    gpgrt(idime)=gpgrt(idime)&
                         +eltem(inode)*cartd(idime,inode)
                    gpcod(idime)=gpcod(idime)&
                         +elcod(idime,inode)*elmar(pelty)%shape(inode,igaus)
                 end do
              end do

              if(kfl_naxis==1) then
                 dvolu=dvolu*gpcod(1)*twopi
              end if

              call rad_exacso(&
                   1_ip,gpcod,dummr,dummr,&
                   dummr,dummr,dummr,dummr,extem,exgrt,dummr)

              ditem        = abs(gptem-extem)
              abtem        = abs(extem)
              err01_rad(1) = err01_rad(1) + ditem*dvolu
              err02_rad(1) = err02_rad(1) + ditem*ditem*dvolu
              err0i_rad(1) = max(err0i_rad(1),ditem)
              err01_rad(2) = err01_rad(2) + abtem*dvolu
              err02_rad(2) = err02_rad(2) + extem*extem*dvolu
              err0i_rad(2) = max(err0i_rad(2),abtem)

              do idime=1,ndime
                 digrt        = abs(gpgrt(idime)-exgrt(idime))
                 abgrt        = abs(exgrt(idime))
                 err11_rad(1) = err11_rad(1) + digrt*dvolu
                 err12_rad(1) = err12_rad(1) + digrt*digrt*dvolu
                 err1i_rad(1) = max(err1i_rad(1),digrt)
                 err11_rad(2) = err11_rad(2) + abgrt*dvolu
                 err12_rad(2) = err12_rad(2) + exgrt(idime)*exgrt(idime)*dvolu
                 err1i_rad(2) = max(err1i_rad(2),abgrt)
              end do

           end do gauss_points_2
        end if
     end do elements_2

     err02_rad(1) = sqrt(err02_rad(1))
     err12_rad(1) = sqrt(err12_rad(1))

     if(err01_rad(2)>zetem) err01_rad(1) = err01_rad(1)/err01_rad(2)        ! L1(T)
     if(err02_rad(2)>zetem) err02_rad(1) = err02_rad(1)/sqrt(err02_rad(2))  ! L2(T)
     if(err0i_rad(2)>zetem) err0i_rad(1) = err0i_rad(1)/err0i_rad(2)        ! Li(T)
     if(err11_rad(2)>zetem) err11_rad(1) = err11_rad(1)/err11_rad(2)        ! L1(gradT)
     if(err12_rad(2)>zetem) err12_rad(1) = err12_rad(1)/sqrt(err12_rad(2))  ! L2(gradT)
     if(err1i_rad(2)>zetem) err1i_rad(1) = err1i_rad(1)/err1i_rad(2)        ! Li(gradT)
     write(momod(modul)%lun_outpu,100) ittim,cutim,&
          err01_rad(1),err02_rad(1),err0i_rad(1),&
          err11_rad(1),err12_rad(1),err1i_rad(1)

  case(2_ip)
     !
     ! Compute error
     !
     do ipoin=1,npoin
        call rad_exacso(&
             1_ip,coord(1,ipoin),dummr,dummr,&
             dummr,dummr,dummr,dummr,extem,dummr,dummr)
        gesca(ipoin)=radav_rad(ipoin,1)-extem
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

end subroutine rad_exaerr
