subroutine ibm_readim()
  !-----------------------------------------------------------------------
  !****f* Immbou/ibm_readim
  ! NAME
  !    ibm_readim
  ! DESCRIPTION
  !    Read IB
  ! OUTPUT
  ! USED BY
  !    ibm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use def_inpout
  use mod_memchk
  use mod_iofile
  use mod_ecoute, only : ecoute
  use mod_elmgeo, only : elmgeo_element_name_to_type
  implicit none
  integer(ip) :: ipara,iimbo,ielty,dummi

  if( INOTSLAVE  ) then
     !
     ! Initializations (defaults)
     !
     nimbo         = 0
     !
     ! Reach the section
     !
     call ecoute('ibm_readim')
     do while( words(1) /= 'DIMEN' )
        call ecoute('ibm_readim')
     end do
     call ecoute('ibm_readim')

     !-------------------------------------------------------------------
     !
     ! Dimensions
     !
     !-------------------------------------------------------------------

     do while (words(1) /= 'ENDDI' )

        if( words(1) == 'NUMBE'  ) then

           nimbo = getint('NUMBE',1_ip,'*NUMBER OF IBM')
           call ibm_defini(0_ip)

        else if( words(1) == 'IB   '  ) then

           iimbo = getint('NUMBE',1_ip,'*PARTICLE NUMBER')
           if( nimbo == 0 )    call runend('IBM_REAPHY: NO PARTICLE')
           if( iimbo > nimbo ) call runend('IBM_REAPHY: WRONG PARTICLE NUMBER')

           call ecoute('ibm_readim')

           do while( words(1) /= 'ENDIB' )

              if( words(1) == 'NODAL' ) then
                 !
                 ! NPOIB: Number of nodes
                 !
                 if( exists('VOLUM') .or.  imbou(iimbo) % kfl_typeb == -1 ) then
                    imbou(iimbo) % npoin = getint('NODAL',0_ip,'*NUMBER OF VOLUME IBM NODAL POINTS')    
                    imbou(iimbo) % npoib = imbou(iimbo) % npoin
                 else
                    imbou(iimbo) % npoib = getint('NODAL',0_ip,'*NUMBER OF BOUNDARY IBM NODAL POINTS')
                    npoib = npoib + imbou(iimbo) % npoib
                 end if

              else if( words(1) == 'BOUND' ) then
                 !
                 ! NPOIB: Number of boundaries
                 !
                 imbou(iimbo) % nboib = getint('BOUND',0_ip,'*NUMBER OF IBM BOUNDARY ELEMENTS')
                 nboib = nboib + imbou(iimbo) % nboib

              else if( words(1) == 'ELEME' ) then
                 !
                 ! NELEM: Number of volume elements
                 !
                 imbou(iimbo) % nelem = getint('ELEME',0_ip,'*NUMBER OF IBM BOUNDARY ELEMENTS')

              else if( words(1) == 'TYPES' ) then
                 !
                 ! LTYIB: types of boundaries
                 !
                 if( nnpar == 0 ) then
                    ipara =  2
                    do while( trim(words(ipara)) /= "" )
                       call elmgeo_element_name_to_type(words(ipara),ielty)
                       if( ielty >= 2 .and. ielty <= nelty ) then
                          lexib(ielty) = 1
                       end if
                       ipara = ipara+1
                    end do
                 else
                    do ipara = 1,nnpar
                       ielty = int(param(ipara))
                       if( ielty >= 2 .and. ielty <= nelty ) then
                          lexib(ielty) = 1
                       end if 
                    end do
                 end if

              else if( words(1) == 'TYPE ' ) then
                 !
                 ! Type of IB: boundary, element or body fitted
                 !
                 if( words(2) == 'EMBED' ) then           ! Embedded with boundary mesh
                    imbou(iimbo) % kfl_typeb =  0
                 else if( words(2) == 'BOUND' ) then      ! Embedded with boundary mesh
                    imbou(iimbo) % kfl_typeb =  0
                 else if( words(2) == 'FORCE' ) then      ! Embedded but does not create holes
                    imbou(iimbo) % kfl_typeb = -2
                 else if( words(2) == 'ELEME' ) then      ! IB with volume mesh
                    imbou(iimbo) % kfl_typeb = -1
                 else if( words(2) == 'VOLUM' ) then      ! IB with volume mesh
                    imbou(iimbo) % kfl_typeb = -1
                 else if( words(2) == 'BODYF' ) then      ! IB is a boudnary set
                    imbou(iimbo) % kfl_typeb = getint('SET  ',1_ip,'*BODY FITTED PARTICLE SET')
                 else 
                    imbou(iimbo) % kfl_typeb =  0
                 end if

              else if( words(1) == 'MODEL' ) then
                 !
                 ! Model of IB: Wind turbine, etc.
                 !
                 if( words(2) == 'WINDT' ) then           ! Wind turbine
                    imbou(iimbo) % kfl_model = 1
                    call ibm_model_wind_turbine(1_ip,iimbo)
                 end if

              else if( words(1) == 'COUPL' ) then
                 !
                 ! Type of coupling: Dirichlet or force term
                 !
                 if( words(2) == 'DIRIC' ) then
                    imbou(iimbo) % kfl_coupl = 0
                 else if( words(2) == 'FORCE' ) then
                    imbou(iimbo) % kfl_coupl = 1
                 end if

              end if

              call ecoute('ibm_readim')
           end do
        end if

        call ecoute('ibm_readim')

     end do
     
     if( nimbo == 0 ) call runend('IBM_READIM: NO IB HAS BEEN DEFINED')
     !
     ! NDIME: Guess dimension (as it is not yet known!)     
     !
#ifdef NDIMEPAR

#else
     ndime = 2 
#endif
     iimbo = 1 
     if( nimbo > 0 ) then
        if( imbou(iimbo) % kfl_typeb == 0 ) then
           do ielty = 10,29
#ifdef NDIMEPAR

#else
              if( lexib(ielty) == 1 ) ndime = 3
#endif           
           end do
        end if
     end if
     ndimb = ndime - 1
     !
     ! IESTA_DOM, etc: Where element and boundary element types start and stop
     !
     if( ndime == 1 ) then
        iesta_dom =  2
        iesto_dom =  9
        ibsta_dom =  1
        ibsto_dom =  1
     else if( ndime == 2 ) then
        iesta_dom = 10
        iesto_dom = 29
        ibsta_dom =  2
        ibsto_dom =  9
     else
        iesta_dom = 30
        iesto_dom = 50
        ibsta_dom = 10
        ibsto_dom = 29
     end if
     !
     ! Treat elements of dimension ndime-1: NGAUS, LQUAD and LEXIS
     !
     do ielty = iesta_dom,iesto_dom
        if( lexib(ielty) == 1 ) then
           call bouele(nelty,dummi,dummi,-ielty,ngaus,dummi,lexib)
        end if
     end do
     !
     ! MNODI and MNOIB
     !
     mnodi = -1
     mnoib = -1
     do ielty = iesta_dom,iesto_dom
        if( lexib(ielty) == 1 ) then
           mnodi = max(mnodi,nnode(ielty))
        end if
     end do
     do ielty = ibsta_dom,ibsto_dom
        if( lexib(ielty) == 1 ) then
           mnoib = max(mnoib,nnode(ielty))
        end if
     end do
     mnodi = max(mnoib,mnodi)
     !if( mnodi < 0 ) mnodi = 9

  end if

end subroutine ibm_readim
