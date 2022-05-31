!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_reabcp.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   
!> @details 
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_reabcp(nbcod,nicod,lbcod,licod,rbcod,ricod)
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_nastal 
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip), intent(in) :: nbcod,nicod,lbcod(50,2),licod(50,2)
  real(rp),    intent(in) :: rbcod(50,10),ricod(50,5)    ! rbcod could contain also info on non-constant bc
  integer(ip)             :: ifint,krang(50,3),irang,kfixi,kfibo,ifoun,ifaul,kfoun,kfaul,ipoin,ibopo
  integer(ip)             :: kcode,nrang,icode,idime,kcods(10),nfoun,ibcod,idofn,ncodf
  integer(ip)             :: icodf,jcodf,kfini(5),kfilo(5),kpoin
  real(rp)                :: vrang(50,6,3),xfixi(10)
  character(4)            :: chfix

  ! Adapt reading to the b.c. file origin (nastal)
  ncodf = ndofn_nsa
  if (kfl_bofty_nsa == 1) ncodf = ndime

  !
  ! positional boundary conditions: reading the conditions
  !
  ifint= 0
  ifaul= 0
  kfaul= 0
  krang= 0
  vrang= 0.0_rp
  irang= 0
  
  if (nbcod == -1 .and. nicod == -1) call runend('nsa_reabcs: USER-DEFINED CODES SECTION NOT PRESENT')
  call ecoute('nsa_reabcs')
  do while(words(1)/='ENDON')
     if (words(1)=='INTER') then
        if (exists('ADDIN')) ifint= 1 
     else if (words(1)=='DEFAU') then
        ifaul = 1
        kfaul = int(param(1))
     else if (words(1)=='BCBOU') then                       
        irang = irang + 1
        do while(words(1)/='ENDBC')
           !
           ! boundary conditions inside a bounding box
           !
           if (words(1)=='SMALL') then 
              vrang(irang,1,1) = param(1)
              vrang(irang,1,2) = param(2)
              vrang(irang,1,3) = param(3)
           else if(words(1)=='LARGE') then                             
              vrang(irang,2,1) = param(1)
              vrang(irang,2,2) = param(2)
              vrang(irang,2,3) = param(3)
           else if(words(1)=='CODEN') then                             
              krang(irang,1) = 1_ip
              krang(irang,2) = int(param(1))
           end if
           call ecoute('nsa_reabcs')
        end do
        
     else if (words(1)=='RANGE') then                       
        irang = irang + 1
        do while(words(1)/='ENDRA')
           !
           ! limiting hyper-surface
           !
           if (words(1)=='CODEN') then                             
              kcode = getint('CODEN',1_ip,'#CODE NUMBER FOR THIS BODY')
              krang(irang,2) = kcode
           else if (words(1)=='COORD') then
              if (words(2)=='XCOOR') krang(irang,3) = 1
              if (words(2)=='YCOOR') krang(irang,3) = 2
              if (words(2)=='ZCOOR') krang(irang,3) = 3
           else if(words(1)=='LOWER') then                             
              krang(irang,1)   =  10
              vrang(irang,1,1) = param(1)
           else if(words(1)=='GREAT') then
              !                             krang(irang,1)   =  krang(irang,1) + 1
              krang(irang,1)   =  20
              vrang(irang,1,1) =  param(1)
           end if
           call ecoute('nsa_reabcs')
        end do
     end if
     call ecoute('nsa_reabcs')
  end do
  
  !
  ! positional boundary conditions: setting the conditions
  !
  nrang = irang

  points: do ipoin = 1,npoin
     
     ibopo = lpoty(ipoin)
     
     !                    if (coord(1,ipoin) > -0.1 .and. coord(1,ipoin) < 0.1) then
     !                       if (coord(2,ipoin) > -0.1 .and. coord(2,ipoin) < 0.1) then
     !                          write(6,*) 'hola',ipoin,ibopo,coord(1:2,ipoin)
     !                          stop
     !                       end if
     !                    end if
     
     if (ibopo > 0 .or. nicod > 0) then
        
        kfixi = -1
        kfibo =  0
        ifoun =  0
        
        conditions: do irang= 1,nrang
           
           kfoun= 0
           if (krang(irang,1) == 1) then   ! bounding box
              if (coord(1,ipoin) > vrang(irang,1,1)) then
                 if (coord(1,ipoin) < vrang(irang,2,1)) then
                    kfoun= kfoun+1
                 end if
              end if
              if ((kfoun == 1) .and. (ndime>=2)) then
                 if (coord(2,ipoin) > vrang(irang,1,2)) then
                    if (coord(2,ipoin) < vrang(irang,2,2)) then
                       kfoun= kfoun+1
                    end if
                 end if
              end if
              if ((kfoun == 2) .and. (ndime==3)) then
                 if (coord(3,ipoin) > vrang(irang,1,3)) then
                    if (coord(3,ipoin) < vrang(irang,2,3)) then
                       kfoun= kfoun+1
                    end if
                 end if
              end if
              if (kfoun == ndime) then
                 ifoun= ifoun + 1
                 kcods(ifoun) = krang(irang,2)
              end if 
             
           else if (krang(irang,1) == 10 ) then   ! lower than the coordinate value
              idime= krang(irang,3)
              if (coord(idime,ipoin) < vrang(irang,1,1)) then
                 ifoun= ifoun + 1
                 kcods(ifoun) = krang(irang,2)
              end if
           else if (krang(irang,1) == 20 ) then   ! greater than the coordinate value
              idime= krang(irang,3)
              if (coord(idime,ipoin) > vrang(irang,1,1)) then
                 ifoun= ifoun + 1
                 kcods(ifoun) = krang(irang,2)
              end if              
           end if
                      
        end do conditions

        nfoun= ifoun
        if ((ifint == 0).and.(nfoun > 1)) nfoun=1
        
        if ((nfoun == 0).and.(ifaul==1)) then     ! a default code is given
           ifoun= 1
           kcods(ifoun) = kfaul
        end if
        
        do ifoun=1,nfoun
           kcode= kcods(ifoun)
           codes_ifoun: do ibcod= 1,nbcod
              if(lbcod(ibcod,1) == kcode) then
                 kfixi          = lbcod(ibcod,2)
                 xfixi(1:ncodf) = rbcod(ibcod,1:ncodf) 
                 xfixi(1+ncodf) = rbcod(ibcod,1+ncodf)   ! needed when non-const bc
!                 write(6,*) 'ueeeeeep',xfixi(1:ncodf+1)
                 kfibo          = 1
                 exit codes_ifoun
              end if
           end do codes_ifoun
           
           if (nicod >= 0) then
              initials_ifoun: do ibcod= 1,nicod
                 if(licod(ibcod,1) == kcode) then
                    kfixi          = licod(ibcod,2)
                    xfixi(1:ncodf) = ricod(ibcod,1:ncodf) 
                    if (kfibo == 1) then
                       call runend('nsa_reabcs: USER-DEFINED INITIAL' &
                            // 'CODE ( '//adjustl(trim(chfix))        &
                            // ' ) IS USED ALREADY AS A BOUNDARY CODE')
                    else
                       kfibo          = 2
                    end if
                    exit initials_ifoun
                 end if
              end do initials_ifoun
           end if
           
           if (kfixi == -1) then
              write (chfix,'(i2)') kcode
              call runend('nsa_reabcs: USER-DEFINED FIXITY OR INITIAL' &
                   // 'CODE ( '//adjustl(trim(chfix))                  &
                   // ' ) NON DEFINED')
           end if
           
           if (kfibo == 1) then
              kfilo(1) = kfixi
              call codfix(ncodf,kfilo)
              jcodf= 0
              do icodf= 1,ncodf
                 if (kfl_fixno_nsa(icodf,ipoin)<= 0) then
                    kfl_fixno_nsa(icodf,ipoin) = kfilo(icodf)
                    bvess_nsa(icodf,ipoin,1)   = xfixi(icodf)
                 end if
              end do
              if (kfl_conbc_nsa == 0) then
                 kfl_funno_nsa(ipoin) = int(xfixi(ncodf+1))
                 if( (kfl_fixno_nsa(    1,ipoin)==1.or.&
                      kfl_fixno_nsa(    2,ipoin)==1.or.&
                      kfl_fixno_nsa(ndime,ipoin)==1).and.kfl_funno_nsa(ipoin)/=0) then
                    bvess_nsa(1:ndime,ipoin,2)=bvess_nsa(1:ndime,ipoin,1)
                 end if
              end if
              
              if(ibopo>0) then

                 kfl_fixrs_nsa(ibopo) = int(param(3))
                 if(kfl_fixno_nsa(1,ipoin)==2 .and. kfl_fixrs_nsa(ibopo)==0) &
                      kfl_fixrs_nsa(ibopo)=-1
                 if(kfl_fixno_nsa(1,ipoin)==5 .and. kfl_fixrs_nsa(ibopo)==0) &
                      kfl_fixrs_nsa(ibopo)=-1
                 if(kfl_fixno_nsa(1,ipoin)==9 .and. kfl_fixrs_nsa(ibopo)==0) &
                      kfl_fixrs_nsa(ibopo)=-1
              end if
              
              if(kfl_visco_nsa == 1 .and. kfl_fixno_nsa(1,ipoin)==2) then
                 !  when viscosity \= 0, change slip nodes to no-slip ones
                 kfl_fixno_nsa(1:ndime,ipoin  ) = 1 
                 bvess_nsa(    1:ndime,ipoin,1) = 0.0_rp                   
                 if (kfl_foreg_nsa == 0) then
                    kfl_fixno_nsa(ndofn_nsa,ipoin  ) = 1 
                    bvess_nsa(    ndofn_nsa,ipoin,1) = tstag_nsa                   
                 end if
              end if
              
           else if (kfibo == 2) then
              kfl_inifi_nsa(1)                = 1                    
              kfini(1)                    =  kfixi  
              call codfix(ncodf,kfini)
              do idofn=1,ncodf
                 if (kfini(idofn) == 1) then
                    bvess_nsa(idofn,ipoin,1)  =  xfixi(idofn)
                 end if
              end do
           end if
           
           
        end do
        
        ! correct local basis: locals are pre-eminent
        
        if (kfl_fixno_nsa(1,ipoin) == 2 .or. kfl_fixno_nsa(1,ipoin) == 9) then
           kfl_fixno_nsa(2,ipoin)     = 0
           kfl_fixno_nsa(ndime,ipoin) = 0
        end if
        
        
     end if

  end do points
  
end subroutine nsa_reabcp
