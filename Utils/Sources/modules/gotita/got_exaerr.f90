subroutine got_exaerr()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_exaerr
  ! NAME 
  !    got_exaerr
  ! DESCRIPTION
  !    Compute the error with respect to an exact solution
  ! USES
  ! USED BY
  !    got_output
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  use mod_memchk
  use mod_postpr
  implicit none 
  real(rp)              :: gpcar(ndime,mnode) 
  real(rp)              :: xjaci(9),xjacm(9) 
  real(rp)              :: elcod(ndime,mnode) 
  real(rp)              :: elvel(ndime,mnode) 
  real(rp)              :: elcdr(mnode),elvdr(ndime,mnode)
  real(rp)              :: gpcod(ndime) 
  real(rp)              :: gpcdr,gpvdr(ndime)
  real(rp)              :: gpvel(ndime)
  real(rp)              :: gpgcd(ndime)
  real(rp)              :: gpgvd(ndime,ndime)
  integer(ip)           :: ielem,inode,ipoin,jdime
  integer(ip)           :: igaus,idime
  integer(ip)           :: pnode,pelty,pgaus
  integer(4)            :: istat
  real(rp)              :: gpvol,gpdet,diffp,abpre,difeu,abvel
  real(rp)              :: erp01(2),erp02(2),erp0i(2)
  real(rp)              :: erp11(2),erp12(2),erp1i(2)
  real(rp)              :: eru01(2),eru02(2),eru0i(2)
  real(rp)              :: eru11(2),eru12(2),eru1i(2) 
  real(rp)              :: dummr(3)
  real(rp), allocatable :: error(:)
  !
  ! Error calculations
  !
  erp01=0.0_rp
  erp02=0.0_rp
  erp0i=0.0_rp
  erp11=0.0_rp
  erp12=0.0_rp
  erp1i=0.0_rp
  eru01=0.0_rp
  eru02=0.0_rp
  eru0i=0.0_rp
  eru11=0.0_rp
  eru12=0.0_rp
  eru1i=0.0_rp

  do ielem=1,nelem
     !
     ! Initialize
     !
     pelty=ltype(ielem)
     pnode=nnode(pelty)
     pgaus=ngaus(pelty)
     !
     ! Gather operations
     !
     do inode=1,pnode
        ipoin=lnods(inode,ielem)
        do idime=1,ndime
           elcod(idime,inode)=coord(idime,ipoin)
           elvdr(idime,inode)=vdrop(idime,ipoin,1)
        end do
        elcdr(inode)=cdrop(ipoin,1)
        call got_veloci(ipoin,elcod(1,inode),elvel(1,inode))
     end do

     do igaus=1,pgaus
        !
        ! Cartesian derivatives and Jacobian
        !
        call elmder(&
             pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
             elcod,gpcar,gpdet,xjacm,xjaci)
        gpvol=elmar(pelty)%weigp(igaus)*gpdet     
        !
        ! Gauss point values
        !
        gpcdr=0.0_rp
        do idime=1,ndime
           gpcod(idime)=0.0_rp
           gpvdr(idime)=0.0_rp
           gpvel(idime)=0.0_rp
        end do
        do inode=1,pnode
           gpcdr=gpcdr+elcdr(inode)*elmar(pelty)%shape(inode,igaus)
           do idime=1,ndime
              gpvdr(idime)=gpvdr(idime)+elvdr(idime,inode)&
                   *elmar(pelty)%shape(inode,igaus)
              gpcod(idime)=gpcod(idime)+elcod(idime,inode)&
                   *elmar(pelty)%shape(inode,igaus)
              gpvel(idime)=gpvel(idime)+elvel(idime,inode)&
                   *elmar(pelty)%shape(inode,igaus)
           end do
        end do
        do idime=1,ndime
           gpgcd(idime)=0.0_rp
           do jdime=1,ndime
              gpgvd(idime,jdime)=0.0_rp
           end do
        end do
        do inode=1,pnode
           do idime=1,ndime
              gpgcd(idime)=gpgcd(idime)&
                   +elcdr(inode)*gpcar(idime,inode)
              do jdime=1,ndime
                 gpgvd(idime,jdime)=gpgvd(idime,jdime)&
                      +elvdr(jdime,inode)*gpcar(idime,inode)                
              end do
           end do
        end do
        !
        ! Get exact solution
        !
        call got_exacso(1_ip,gpcod,gpvel,dummr)        
        !
        ! Compute errors
        !
        diffp    = abs(gpcdr-excdr_got)
        abpre    = abs(excdr_got)
        erp01(1) = erp01(1) + diffp*gpvol
        erp02(1) = erp02(1) + diffp*diffp*gpvol
        erp0i(1) = max(erp0i(1),diffp)
        erp01(2) = erp01(2) + excdr_got*gpvol
        erp02(2) = erp02(2) + excdr_got*excdr_got*gpvol
        erp0i(2) = max(erp0i(2),excdr_got)
        do idime=1,ndime
           difeu    = abs(gpvdr(idime)-exvdr_got(idime))
           abvel    = abs(exvdr_got(idime))
           eru01(1) = eru01(1) + difeu*gpvol
           eru02(1) = eru02(1) + difeu*difeu*gpvol
           eru0i(1) = max(eru0i(1),difeu)
           eru01(2) = eru01(2) + abvel*gpvol
           eru02(2) = eru02(2) + abvel*abvel*gpvol
           eru0i(2) = max(eru0i(2),abvel)
           diffp    = abs(gpgcd(idime)-exgcd_got(idime))
           abpre    = abs(exgcd_got(idime)) 
           erp11(1) = erp11(1) + diffp*gpvol
           erp12(1) = erp12(1) + diffp*diffp*gpvol
           erp1i(1) = max(erp1i(1),diffp)
           erp11(2) = erp11(2) + abpre*gpvol
           erp12(2) = erp12(2) + abpre*abpre*gpvol
           erp1i(2) = max(erp1i(2),abpre)
           do jdime=1,ndime
              difeu    = abs(gpgvd(idime,jdime)-exgvd_got(idime,jdime))
              abvel    = abs(exgvd_got(idime,jdime))
              eru11(1) = eru11(1) + difeu*gpvol
              eru12(1) = eru12(1) + difeu*difeu*gpvol
              eru1i(1) = max(eru1i(1),difeu)
              eru11(2) = eru11(2) + abvel*gpvol
              eru12(2) = eru12(2) + abvel*abvel*gpvol
              eru1i(2) = max(eru1i(2),abvel)
           end do
        end do

     end do

  end do

  erp02(1) = sqrt(erp02(1))
  erp12(1) = sqrt(erp12(1))
  eru02(1) = sqrt(eru02(1))
  eru12(1) = sqrt(eru12(1))

  if(erp01(2)>0.0_rp) erp01(1) = erp01(1)/erp01(2) 
  if(erp02(2)>0.0_rp) erp02(1) = erp02(1)/sqrt(erp02(2))
  if(erp0i(2)>0.0_rp) erp0i(1) = erp0i(1)/erp0i(2)
  if(eru01(2)>0.0_rp) eru01(1) = eru01(1)/eru01(2) 
  if(eru02(2)>0.0_rp) eru02(1) = eru02(1)/sqrt(eru02(2))
  if(eru0i(2)>0.0_rp) eru0i(1) = eru0i(1)/eru0i(2) 
  if(erp11(2)>0.0_rp) erp11(1) = erp11(1)/erp11(2) 
  if(erp12(2)>0.0_rp) erp12(1) = erp12(1)/sqrt(erp12(2))
  if(erp1i(2)>0.0_rp) erp1i(1) = erp1i(1)/erp1i(2)
  if(eru11(2)>0.0_rp) eru11(1) = eru11(1)/eru11(2) 
  if(eru12(2)>0.0_rp) eru12(1) = eru12(1)/sqrt(eru12(2))
  if(eru1i(2)>0.0_rp) eru1i(1) = eru1i(1)/eru1i(2) 

  write(momod(modul)%lun_outpu,100) ittim,itcou,&
       eru01(1),erp01(1),eru02(1),erp02(1),eru0i(1),erp0i(1),&
       eru11(1),erp11(1),eru12(1),erp12(1),eru1i(1),erp1i(1)
  !
  ! Output and Postprocess error
  !
  if(postp(1)%npp_stepi(4)/=0) then
     allocate(error(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ERROR','got_exaerr',error)
     do ipoin=1,npoin
        call got_veloci(ipoin,coord(1,ipoin),gpvel)
        call got_exacso(1_ip,gpcod,gpvel,dummr)
        error(ipoin)=abs(excdr_got-cdrop(ipoin,1))
     end do
     !call postpr(error,postp(1)%wopos(1:2,4),ittim,cutim-dtime)
     call memchk(two,istat,mem_modul(1:2,modul),'ERROR','got_exaerr',error)
     deallocate(error,stat=istat)
     if(istat/=0) call memerr(two,'ERROR','got_exaerr',0_ip)
  end if

100 format(///,10X,'FINITE ELEMENT ERRORS',                              &
       &              /,10X,'=====================',//,                           &
       &  '          TIME STEP NO.',I5,',  ITERATION NO. ',I5,/,10X,40('-'),/,    &
       &  '          NORM       VELOCITY             WATER F. ',/,10X,40('-'),/,  &
       &  '          W(0,1) ',E12.5,9X,E12.5 ,/, &
       &  '          W(0,2) ',E12.5,9X,E12.5 ,/, &
       &  '          W(0,i) ',E12.5,9X,E12.5 ,/, &
       &  '          W(1,1) ',E12.5,9X,E12.5 ,/, &
       &  '          W(1,2) ',E12.5,9X,E12.5 ,/, &
       &  '          W(1,i) ',E12.5,9X,E12.5 ,/,5X,40('-'))
101 format(///,10X,'SUBGRID SCALE NORM',&
       &     /,10X,'==================',//,                           &
       &  '          TIME STEP NO.',I5,',  ITERATION NO. ',I5,/,10X,40('-'),/,    &
       &  '          NORM       VELOCITY ',/,10X,40('-'),/,  &
       &  '          W(0,2) ',E12.5)

end subroutine got_exaerr
