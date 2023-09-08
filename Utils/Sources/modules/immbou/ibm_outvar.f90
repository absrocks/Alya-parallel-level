!-----------------------------------------------------------------------
!> @addtogroup Immbou
!> @{
!> @file    ibm_outvar.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Output for postprocess
!> @details Output for postprocess
!> @} 
!-----------------------------------------------------------------------
subroutine ibm_outvar(ivari)
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use def_postpr
  use def_elmtyp
  use mod_kdtree
  use mod_intpol

  implicit none
  integer(ip), intent(in) :: ivari  !> checking variable ivari 
  integer(ip)             :: iimbo,kpoin,idime,ipoib,ipoin,jpoin,ncoun,iblty
  integer(ip)             :: iwaib,pblty,inodb,iboun,kboun
  integer(ip), save       :: ipass=0
  character(5)            :: wopo2(5)
  integer(ip)             :: ibopo
  integer(ip)             :: dummi,ielem,inode
  integer(ip)             :: iinte,ninte,izdom
  real(rp)                :: deriv(3,64),coloc(3),shaib(mnode)
  real(rp)                :: dummr,dumma(ndime),propo(ndime)
  real(rp)                :: rutim,total,dista,proje(3),dmini,dmaxi

  if( ivari == 0 ) return
  !
  ! Define postprocess variable
  !
  rutim = cutim

  select case (ivari)  

  case(0_ip)

     return

  case(1_ip)

     !-------------------------------------------------------------------
     !
     ! DMESH
     !
     !-------------------------------------------------------------------

     if( INOTSLAVE ) then

        call memgen(0_ip,ndime,npoib)

        kpoin = 0
        do iimbo = 1,nimbo
           do ipoib = 1,imbou(iimbo)%npoib
              kpoin = kpoin + 1
              do idime = 1,ndime
                 gevec(idime,kpoin) = imbou(iimbo)%cooib(idime,ipoib) - imbou(iimbo)%cooin(idime,ipoib)
              end do
           end do
        end do

        !if( kfl_outfo == -1 ) then
        if( 1 == 1 ) then
           !
           ! GiD format
           !
           if( ipass == 0 ) then
              call ibm_openfi(2_ip)
              write(lun_resib,'(a)')'GiD Post Results File 1.0'
              write(lun_resib,'(a)')' '
           end if
           ipass = 1 
           write(lun_resib,2) 'DMESH','ANALYSIS',cutim,'Vector'
           write(lun_resib,3) 'DMESH_X','DMESH_Y','DMESH_Z'
           write(lun_resib,1) 'Values'
           do ipoin = 1,npoib 
              write(lun_resib,4) npoin+ipoin,(gevec(idime,ipoin),idime=1,ndime)
           end do
           write(lun_resib,1) 'End Values'

        else if( kfl_outfo == -3000 ) then
           !
           ! Alya Binary format
           !           
           call ibm_openfi(2_ip)
           write(lun_resib) ( (gevec(idime,ipoin),idime=1,ndime),ipoin=1,npoib)
           call ibm_openfi(4_ip)

        else if( kfl_outfo == -30 ) then
           !
           ! VU format
           !
           call ibm_openfi(2_ip)
           wopo2(1) = 'DMESX'
           wopo2(2) = 'DMESY'
           wopo2(3) = 'DMESZ'
           write(lun_resib,205) cutim
           iblty = imbou(1)%ltyib(1)
           ncoun = 0
           do idime = 1,ndime
              ncoun = ncoun + 4_ip
              write(lun_resib,200) trim(wopo2(idime)),trim(fil_resi2),npoib,ncoun
              ncoun = ncoun + 4_ip + rp * npoib
              write(lun_resib,210)     
              write(lun_resib,220) trim(wopo2(idime)),trim(cepos(iblty)),trim(wopo2(idime)),'ConnecIB','ZoneIB'
              write(lun_resib,230)
              write(lun_resi2) (gevec(idime,ipoin),ipoin=1,npoib)
           end do
           call ibm_openfi(4_ip)

        end if

        call memgen(2_ip,ndime,npoib)

     end if
     return

  case(2_ip)

     !-------------------------------------------------------------------
     !
     ! DISPM
     !
     !-------------------------------------------------------------------

     gevec => dispm(:,:,1)   

  case(3_ip)

     gesca => lndib_ibm

  case(4_ip)

     !-------------------------------------------------------------------
     !
     ! WALLS
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do iwaib = 1,nwaib
           do iboun = 1,twall_ibm(iwaib) % nboun
              pblty = twall_ibm(iwaib) % ltypb(iboun)
              do inodb = 1,nnode(pblty)
                 kpoin = twall_ibm(iwaib) % lnodb(inodb,iboun)
                 ipoin = twall_ibm(iwaib) % lninv(kpoin)
                 gesca(ipoin) = real(iwaib)
              end do
           end do
        end do
     end if

  case(5_ip)

     !-------------------------------------------------------------------
     !
     ! PREIB
     !
     !-------------------------------------------------------------------
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)      
        do ipoin = 1,npoin
           if (lntib(ipoin) < 0) then
              iimbo = abs(lntib(ipoin))


                 !
                 ! Find the number of free neighbor nodes
                 !
                 ninte = 0_ip
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    if ( lntib(jpoin) <= 0 ) then
                       ninte = ninte + 1_ip
                    end if
                 end do
                 !
                 ! Find the number of free neighbor nodes in other subdomains
                 !
                 if (npoin_2 > npoin .and. ipoin >= npoi1+1 .and. ipoin <= npoin) then
                    do izdom = r_dom_2(ipoin-npoi1),r_dom_2(ipoin-npoi1+1)-1
                       jpoin = c_dom_2(izdom)           
                       if ( lntib(jpoin) <= 0 ) then
                          ninte = ninte + 1
                       end if
                    end do
                 end if
                 !
                 ! Allocate the variables
                 !
                 allocate( lnint(ipoin) % lnode(ninte) )
                 allocate( lnint(ipoin) % shapl(ninte+1) )
                 lnint(ipoin) % limit = ninte
                 !
                 ! Find the free neighbor nodes
                 !
                 iinte = 0
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)                    
                    if ( lntib(jpoin) <= 0 ) then
                       iinte = iinte + 1                          
                       lnint(ipoin) % lnode(iinte) = jpoin                       
                    end if
                 end do
                 !
                 ! Find the free neighbor nodes in other subdomains
                 !
                 if (npoin_2 > npoin .and. ipoin >= npoi1+1 .and. ipoin <= npoin) then
                    do izdom = r_dom_2(ipoin-npoi1),r_dom_2(ipoin-npoi1+1)-1
                       jpoin = c_dom_2(izdom)           
                       if ( lntib(jpoin) <= 0 ) then
                          iinte = iinte + 1                          
                          lnint(ipoin) % lnode(iinte) = jpoin
                       end if
                    end do
                 end if
                 !
                 ! Find the kriging interpolation coefficients
                 !
                 iimbo    = abs(lntib(ipoin))
                 call dpopar(1_ip,coord(1:ndime,ipoin),&
                      imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                      imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                      dummr,dumma,propo,dummi,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                      imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                      imbou(iimbo) % lnele)
                 
                 call krigin(coord,propo,lnint(ipoin) % limit,lnint(ipoin) % lnode,lnint(ipoin) % shapl)
                 
                 gesca(ipoin) = 0
                 do iinte = 1,lnint(ipoin) % limit
                    jpoin = lnint(ipoin) % lnode(iinte)
                    gesca(ipoin) = gesca(ipoin) + lnint(ipoin) % shapl(iinte)*press(jpoin,1)
                 end do

              else
                 gesca(ipoin) = press(ipoin,1)
              end if
           end do
        end if

  case(6_ip)

     !-------------------------------------------------------------------
     !
     ! VELIB
     !
     !-------------------------------------------------------------------
     if( INOTMASTER ) then

        call memgen(zero,ndime,npoin)

        do ipoin = 1,npoin
           if (lntib(ipoin) < 0) then
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do                            
              do inode = 2,lnint(ipoin) % limit
                 jpoin = lnint(ipoin) % lnode(inode)                 
                 do idime = 1,ndime
                    gevec(idime,ipoin) = gevec(idime,ipoin) - &
                         veloc(idime,jpoin,1)*(lnint(ipoin) % shapl(inode)/lnint(ipoin) % shapl(1))
                 end do
              end do
              do idime = 1,ndime
                 gevec(idime,ipoin) = gevec(idime,ipoin) + &
                      veloc(idime,ipoin,1)*(1.0_rp/lnint(ipoin) % shapl(1))
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = veloc(idime,ipoin,1)
              end do
           end if
        end do       
     end if

  case(7_ip)

     !-------------------------------------------------------------------
     !
     ! DISIB
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then

        call memgen(zero,ndime,npoin)
        do ipoin = 1,npoin
           if (lntib(ipoin) < 0) then
              iimbo     = abs(lntib(ipoin))
              
              call dpopar(1_ip,coord(1:ndime,ipoin),&
                   imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                   imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                   dummr,dumma,propo,dummi,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                   imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                   imbou(iimbo) % lnele)

              do idime = 1,ndime
                 gevec(idime,ipoin) =  propo(idime) - coord(idime,ipoin) 
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end do
     end if

  case(8_ip)

     !-------------------------------------------------------------------
     !
     ! DISTANCE
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then

        dmini =  1.0e9_rp
        dmaxi = -1.0e9_rp
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = 1.0e9_rp
        end do
        do ielem = 1,nelem
           if( lelch(ielem) == ELCUT ) then
              iimbo = letib(ielem)
              if( iimbo > 0 ) then
                 do inode = 1,lnnod(ielem)
                    ipoin = lnods(inode,ielem)    
                    call dpopar(1_ip,coord(1:ndime,ipoin),&
                         imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                         imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooib,&
                         dummr,dumma,propo,dummi,imbou(iimbo) % fabox,imbou(iimbo) % sabox,&
                         imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist,&
                         imbou(iimbo) % lnele)                    
                    gesca(ipoin) = dista
                    if( dista < dmini ) dmini = dista
                    if( dista > dmaxi ) dmaxi = dista
                 end do
              end if
           end if
        end do
        do ipoin = 1,npoin        
           if( gesca(ipoin) > 1.0e8_rp ) then
              if( lntib(ipoin) > 0 ) then
                 gesca(ipoin) = dmini
              else
                 gesca(ipoin) = dmaxi
              end if
           end if
        end do
     end if

  end select

  call outvar(&
       ivari,&
       ittim,rutim,postp(1) % wopos(1,ivari))
  !
  ! GiD formats
  !
1 format(a)
2 format('Result ',a,' ',a,' ',e14.8,' ',a,' OnNodes')
3 format('ComponentNames ',a,',',a,',',a)
4 format(i9, 3(1x,e16.8E3))
  !
  ! VU format
  !
200 format('FIELD<double> ',a,'("',a,'",',i12,',',i12,');')
205 format('TEXTE Time(" ',e12.6,'");')
210 format('SOLUTION Solution( ) =',/,'{')
220 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,');')
230 format('};')

end subroutine ibm_outvar
