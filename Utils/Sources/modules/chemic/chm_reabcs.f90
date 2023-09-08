subroutine chm_reabcs()
  !-----------------------------------------------------------------------
  !****f* partis/chm_reabcs
  ! NAME
  !    chm_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions
  ! OUTPUT 
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_chemic
  use mod_opebcs
  use mod_ecoute, only      : ecoute
  implicit none
  integer(ip)  :: iclas,ipara,ispec
  !
  ! Allocate memory
  !

  if( kfl_icodn > 0 ) then
     call opnbcs(0_ip,nclas_chm,1_ip,0_ip,tncod_chm) 
  end if
  if( kfl_icodb > 0 ) then
     call opbbcs(0_ip,nclas_chm,1_ip,tbcod_chm)      
  end if
  call chm_membcs(3_ip)

  if( INOTSLAVE ) then

     kfl_allcl_chm = 0
     do iclas = 1,nclas_chm
        kfl_usrbc_chm(iclas) = 0      
     end do
     do iclas = 1,nspec_chm
        kfl_initi_chm(iclas) = 0   
        xinit_chm(iclas,:)     = 0.0_rp
     end do

     call ecoute('chm_reabcs')
     do while( words(1) /= 'BOUND')
        call ecoute('chm_reabcs')
     end do
     !
     ! Read data
     !
     call ecoute('chm_reabcs')
     do while( words(1) /= 'ENDBO')

        if( words(1) == 'CODES' .and. exists('NODES') ) then 
           !
           ! User-defined codes on nodes
           !
           if( exists('CLASS') ) then
              if( exists('ALL  ') ) then
                 iclas = 1
                 kfl_allcl_chm = 1
              else
                 iclas = getint('CLASS',1_ip,  '#Class number')
                 kfl_allcl_chm = 0                 
              end if
           else if( exists('SPECI') ) then
              kfl_allcl_chm = 0
              iclas=0
              do ispec=1,nspec_chm
                 if ( exists(speci(ispec)%name) ) iclas = ispec
              enddo
           else
              call runend('BOUNDARY CODES NEEDS CLASS NUMBER')
           end if

           if( iclas < 1 .or. iclas > nclas_chm ) then
              call runend('class number is wrong')
           end if
           tncod => tncod_chm(iclas:)
           call reacod(1_ip)

        else if( words(1) == 'CODES' .and. exists('BOUND') ) then
           !
           ! User-defined codes on boundaries
           !          
           if( exists('CLASS') ) then
              iclas = getint('CLASS',1_ip,  '#Class number')
              if( iclas >= 1 .and. iclas <= nclas_chm ) then
                 tbcod => tbcod_chm(iclas:)
                 call reacod(2_ip)
              end if
           else
              call runend('BOUNDARY CODES NEEDS CLASS NUMBER')
           end if

        else if( words(1) == 'INITI' ) then
           !
           ! Initial conditions
           !
           call ecoute('chm_reabcs')
           do while( words(1) /= 'ENDIN')
              if( words(1) == 'CLASS' ) then

                 if( exists('CONST') ) then

                    if( exists('CLASS') ) then
                       iclas = getint('CLASS',1_ip,  '#Class number')
                       if( exists('EQUIL') ) then
                          kfl_initi_chm(iclas) = 1
                       else
                          kfl_initi_chm(iclas) = 2
                          xinit_chm(iclas,1) = getrea('CONST',0.0_rp,'#Class constant initial value')
                       end if
                       do while( words(1) /= 'ENDCL') 
                          call ecoute('chm_reabcs')
                       end do
                    end if

                 else if( exists('USERF') ) then
                    iclas = getint('CLASS',1_ip,'#Class number')
                    kfl_usrbc_chm(iclas) = getint('USERF',1_ip,'#User function')
                    if( iclas < 1 .or. iclas > nclas_chm )&
                         call runend('CHM_REABCS: WRONG CLASS NUMBER')
                    do while( words(1) /= 'ENDCL') 
                       call ecoute('chm_reabcs')
                    end do

                 end if

              else if( words(1) == 'SPECY' ) then

                 if( exists('CONST') ) then

                    if( exists('SPECY') ) then
                       iclas = getint('SPECY',1_ip,  '#Specy number')
                       kfl_initi_chm(iclas) = 2
                       xinit_chm(iclas,1) = getrea('CONST',0.0_rp,'#Class constant initial value')
                       do while( words(1) /= 'ENDSP') 
                          call ecoute('chm_reabcs')
                       end do
                    end if
                 end if

              else if( words(1) == 'CONST' ) then

                 if( exists('ONSET') ) then
                    do iclas = 1,nspec_chm 
                       if( exists(speci(iclas)%name(1:5)) ) then
                          kfl_initi_chm(iclas) = -getint('ONSET',-1_ip,'#prescribe value on set only')
                          xinit_chm(iclas,1) = getrea(speci(iclas)%name,0.0_rp,'#Species constant initial value') 
                       end if
                    enddo
                 else  if( exists('BIPHA') .or. exists('LEVEL') ) then
                    if (exists('BIPHA')) then
                       do iclas = 1,nspec_chm 
                          kfl_initi_chm(iclas) = 4
                       enddo
                    else if (exists('LEVEL')) then
                       do iclas = 1,nspec_chm 
                          kfl_initi_chm(iclas) = 3
                       enddo
                    endif
                    do while( words(1) /= 'ENDPH') 
                       if (words(1)=='PHASE') then
                          ipara = getint('PHASE',1_ip,'Phase for initial value')
                          do iclas = 1,nspec_chm 
                             if( exists(speci(iclas)%name(1:5)) ) then
                                xinit_chm(iclas,ipara) = getrea(speci(iclas)%name,0.0_rp,'#Species constant initial value') 
                             end if
                          enddo
                       endif
                       call ecoute('chm_reabcs')
                    end do
                 else
                    do iclas = 1,nspec_chm 
                       if( exists(speci(iclas)%name(1:5)) ) then
                          kfl_initi_chm(iclas) = 2
                          xinit_chm(iclas,1) = getrea(speci(iclas)%name,0.0_rp,'#Species constant initial value') 
                       end if
                    enddo
                 end if

              end if

              call ecoute('chm_reabcs')
           end do

        else if( words(1) == 'NATUR' ) then
           !
           ! Parameters for natural boundary conditions
           !
           call ecoute('chm_reabcs')
           iclas=0
           do while( words(1) /= 'ENDNA')
              iclas=iclas+1
              if(iclas>nclas_chm)         call runend('TOO MANY CLASSES FOR NATURAL B.C. PARAMETERS')              
              if(nnpar>size(panat_chm,1)) call runend('TOO MANY PARAMETERS FOR NATURAL B.C. PARAMETERS')
              do ipara=1,nnpar
                 panat_chm(ipara,iclas)=param(ipara)
              end do
              call ecoute('chm_reabcs')
           end do

        end if
        call ecoute('chm_reabcs')
     end do

  end if

end subroutine chm_reabcs
