!-----------------------------------------------------------------------
!> @addtogroup Immbou
!> @{
!> @file    ibm_chgaib.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Compute the host elements of the Gauss points for the IBM. 
!> @details Compute the host elements of the Gauss points for the IBM. 
!>          In parallel, if several subdomains find a host element
!>          for the same Gauss point (due to finite tolerance), only the
!>          subdomain with the lowest number keep it.
!> @} 
!-----------------------------------------------------------------------
subroutine ibm_chgaib()
  use def_master
  use def_kermod
  use def_domain
  use def_immbou
  use mod_memchk
  use mod_messages, only : livinf
  implicit none  
  integer(ip)          :: igaub,pgaub,idime,ielem,inodb,pnodb,ipoin,mdang
  integer(ip)          :: iimbo,iboun,pblty,kgaus,dummi,ielpo,idang,kgaub
  integer(4)           :: istat
  real(rp)             :: posgb(3),bocod(ndime,mnoib)
  real(rp)             :: deriv(3,64),coloc(3),shaib(mnode)
  integer(ip), save    :: ipass = 0
  integer(ip), target  :: qgaub(1),ndang(1)
  integer(ip), pointer :: danel(:)

  !
  ! KGAUB = total number of Gauss points
  !
  if( ittim /= 0 ) call livinf(165_ip,'GAUSS',0_ip) 
  if( kfl_coibm == 0 ) return

  kgaub    = 0
  qgaub(1) = 0
  posgb    = 0.0_rp
  coloc    = 0.0_rp
  
  do iimbo = 1,nimbo
     do iboun = 1,imbou(iimbo) % nboib
        pblty = imbou(iimbo) % ltyib(iboun)
        kgaub = kgaub + ngaib(pblty)
     end do
  end do

  !
  ! Host element of Gauss points. QGAUB = number of Gauss point with host element
  !
  if( INOTMASTER ) then
    
     if( ipass > 0 ) then
        call memchk(two,istat,memor_dom,'LGAIB','ibm_chgaib',lgaib)
        deallocate(lgaib,stat=istat)
        if(istat/=0) call memerr(two,'LGAIB','ibm_chgaib',0_ip)
     end if
     allocate(lgaib(kgaub),stat=istat)   
     call memchk(zero,istat,memor_dom,'LGAIB','ibm_chgaib',lgaib)
     ipass    = ipass + 1
     kgaus    = 0

     !
     ! Move nodes with mesh velocity different to zero (ALE technique)
     !
     !if( kfl_aleib_ibm > 0 ) then
     !   do ipoin = 1,npoin        
     !      do idime = 1,ndime
     !         coord(idime,ipoin) = coord(idime,ipoin) + dispm(idime,ipoin,1)
     !      end do
     !   end do
     !end if

     do iimbo = 1,nimbo
        do iboun = 1,imbou(iimbo)%nboib

           pblty    = imbou(iimbo)%ltyib(iboun)
           pnodb    = nnode(pblty)
           pgaub    = ngaib(pblty)
           do inodb = 1,pnodb
              ipoin = imbou(iimbo)%lnoib(inodb,iboun)
              do idime = 1,ndime
                 bocod(idime,inodb) = imbou(iimbo)%cooib(idime,ipoin)
              end do
           end do

           do igaub = 1,pgaub

              kgaus = kgaus + 1
              do idime = 1,ndime
                 posgb(idime) = 0.0_rp
                 do inodb = 1,pnodb
                    posgb(idime) = posgb(idime) &
                         + elmar(pblty)%shaib(inodb,igaub) * bocod(idime,inodb)
                 end do
              end do

              call runend('IBM NOT CODED')
              !call elsest(&
              !     2_ip,1_ip,ielse,mnode,ndime,npoin,nelem,nnode(1:),&
              !     lnods,ltype,ltopo,coord,posgb,relse,ielem,&
              !     shaib,deriv,coloc,dummi)

              if( ielem > 0 ) then
                 qgaub(1) = qgaub(1) + 1
                 lgaib(kgaus) = ielem
              else if( ISEQUEN ) then
                 print*,''
                 print*,'IBM_CHKAIB: NO ELEMENT FOUND',posgb(1:ndime)
                 call runend('IBM_CHKAIB: NO ELEMENT FOUND')
              end if

           end do
        end do
     end do

     !
     ! Put nodes with mesh velocity diffrente to zero in the original position (ALE technique)
     !
     !if( kfl_aleib_ibm > 0 ) then
     !   do ipoin = 1,npoin   
     !      do idime = 1,ndime
     !         coord(idime,ipoin) = coord(idime,ipoin) - dispm(idime,ipoin,1)
     !      end do
     !   end do
     !end if
  end if

  !
  ! Total number of Gauss points with host elements
  !
  call parari('SUM',0_ip,1_ip,qgaub)

  !
  ! Error message
  !
  if( qgaub(1) > kgaub ) then
     call livinf(71_ip,' ',0_ip)
  else if( qgaub(1) < kgaub ) then
     call runend('IBM_CHGAIB: HOST ELEMENTS NOT FOUND')
  end if
  !
  ! Sequential version can get out
  !
  if( ISEQUEN ) return
  !
  ! Take out Gauss point repeated over different subdomains
  !   
  if( ISLAVE .and. qgaub(1) /= kgaub ) then

     !
     ! DANEL = Identify border (dangerous) elements => DANEL = -1
     !
     allocate(danel(0:nelem),stat=istat)   
     call memchk(zero,istat,memor_dom,'DANEL','ibm_chgaib',danel)
     danel(0) = 0
     do ipoin = npoi1+1,npoin
        do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
           ielem = lelpo(ielpo)
           danel(ielem) = -1
        end do
     end do

     !
     ! NDANG = Number of dangerous Gauss points
     !
     ndang(1) = 0
     do kgaus = 1,kgaub
        ielem = lgaib(kgaus)
        if( danel(ielem) < 0 ) ndang(1) = ndang(1) + 1
     end do
     !
     ! I have NDANG damgerous Gauss point: tell it to my neighbors
     ! Subdomains with lowest number keep the Gauss point
     !


     mdang = max(1_ip,ndang(1))
     call memgen(1_ip,mdang,0_ip)
     idang = 0
     do kgaus = 1,kgaub
        if( danel(lgaib(kgaus)) < 0 ) then
           idang        =  idang + 1
           gisca(idang) =  kgaus
        end if
     end do
     parin => ndang 
     npari =  ndang(1)
     pari1 => gisca
     call par_slesca()
     !
     ! LGAIB = Put to zero Gauss points that were taken out
     !

     do idang = 1,ndang(1)
        if( gisca(idang) < 0 ) then
           kgaus        = -gisca(idang)
           lgaib(kgaus) = 0
        end if
     end do
     call memgen(1_ip,mdang,0_ip)
     call memchk(two,istat,memor_dom,'DANEL','ibm_chgaib',danel)
     deallocate(danel,stat=istat)
     if(istat/=0) call memerr(two,'DANEL','ibm_chgaib',0_ip)  
        

  end if
  !
  ! Recheck (just in case)
  !
  if( INOTMASTER ) then
     qgaub(1) = 0
     do kgaus = 1,kgaub
        if( lgaib(kgaus) > 0 ) qgaub(1) = qgaub(1) + 1
     end do
  end if
  if( IPARALL ) then
     npari =  1
     parin => qgaub
     call par_operat(3_ip)
     if( qgaub(1) /= kgaub ) then
        if( IMASTER ) print*,qgaub(1),kgaub
        call runend('IBM_CHGAIB: GAUSS POINTS MAY BE LOST OR COUNTED SEVERAL TIMES')
     end if
  end if


end subroutine ibm_chgaib
