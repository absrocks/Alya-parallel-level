subroutine lev_parall(itask)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_parall
  ! NAME
  !    lev_parall
  ! DESCRIPTION
  !    This routine is a bridge to Parall service  
  ! USED BY
  !    lev_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_levels
  use      def_parall
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ii

  if( IPARALL ) then

     select case(itask)

     case(1)
        !
        ! Exchange data read in lev_reaphy, lev_reanut and lev_reaous
        ! always using MPI, even if this is a partition restart
        !
        call lev_sendat(1_ip)

     case(2)
        !
        ! Exchange data read in lev_reabcs
        !
        call lev_sendat(2_ip)

     case(3) 
        !
        ! Sum up residual contribution of slave neighbors
        !
        if( ISLAVE ) then
           call vocabu(NPOIN_REAL_1DIM,0_ip,0_ip)
           parr1 => unkno
           call par_slexch()
        end if

     case(4)
        !
        ! Sum up residual contribution of slave neighbors
        !
        if(kfl_paral>0) then
           call vocabu(NPOIN_REAL_1DIM,0_ip,0_ip)
           parr1 => rhsid
           call par_slexch()
        end if

        !    case(5)
        !      !
        !      ! Modified mass matrix for absorbing b.c.'s
        !      !
        !      if(kfl_paral>0) then
        !         party =  3                     ! nodes
        !         pardi =  1                     ! 1-dimensions
        !         parki =  5                     ! Real number (2 dimensions in one)
        !        pard1 =  1
        !        parr1 => vmass_lev(1:npoin)
        !         call par_slexch()
        !         parr1 => null()
        !      end if

     case(6)

        ! 
        ! nredm_lev writing
        ! master will know the interface data (vector size)
        ! of each slave 
        if( IMASTER ) then
           
           pard1=npart_par                                     ! Bridge
           npart=npart_par

           ! master allocates vector 
           ! nredm_lev to receive interface data
           ! from processor npart_par 
           !
           ! Master: NREDM_LEV(2,PARD1)
           !
           call lev_memall(5_ip)
           nredi_lev = 0

           ! loop to receive from each slave  
           do pard2 = 1,pard1
              npari =  2
              kfl_desti_par=pard2
              call parari('RCV',0_ip,2_ip,nredm_lev(1,pard2))
           end do

        else if(kfl_paral>0) then


           !  each slave sends its data  
           npari =  2
           parin => nredi_lev
           kfl_desti_par=0                                     ! Send to master
           call par_sendin()
            
        end if


     case(7)
        !
        ! Discrete interface vector size reduction
        ! (sizes for the discrete interface on the whole domain
        ! have to be computed)
        npari =  2
        parin => nredt_lev

        call par_operat(3_ip)

     case(8)

        ! master gathers total interface points coordinates vector 
        if( IMASTER ) then
           ! master prepares to receive
           ! each interface points coordinates vector part
           ! from each slave
           ii=1
           pard1=npart_par                                     ! Bridge
           npart=npart_par
           
           do pard2=1,pard1 

              nparr =  ndime*nredm_lev(2,pard2)
              kfl_desti_par=pard2
              if( nredm_lev(2,pard2) /= 0 ) then
                 call pararr('RCV',0_ip,ndime*nredm_lev(2,pard2),coord_lev(1,ii))
              end if
              ii=ii+nredm_lev(2,pard2)
           enddo


        else if( ISLAVE ) then

           ! each slave  prepares to send
           ! its interface points coordinates vector part
           if( nredi_lev(2) /= 0 ) then
              call pararr('S2M',0_ip,ndime*nredi_lev(2),coorp_lev)
           end if

        end if

     case(9)

        ! master sends total interface points coordinates vector 
        ! to each slave
        call pararr('BCT',0_ip,ndime*npoin_lev,coord_lev)
 
     case(10)

        ! master gathers total information about which element (whole mesh numeration) each surface belongs to 
        if( IMASTER ) then
           ! master prepares to receive
           ! each  part
           ! from each slave
           ii=1
           pard1=npart_par                                     ! Bridge
           npart=npart_par

           do pard2=1,pard1 

              nparr =  nredm_lev(1,pard2)    ! number of surfaces of each subdomain
              kfl_desti_par=pard2
              if(  nredm_lev(1,pard2) /= 0 ) then          ! ask guillaume why in case(8) they do no use nparr directly    
                 call parari('RCV',0_ip,nredm_lev(1,pard2),lebsu_lev(ii))     ! here also
              end if
              ii=ii+nredm_lev(1,pard2)    ! and here, check if what I did is right
           enddo


        else if( ISLAVE ) then

           ! each slave  prepares to send
           ! its vector part
           if( nredi_lev(1) /= 0 ) then
              call parari('S2M',0_ip,nredi_lev(1),lebsp_lev)
           end if

        end if

     case(11)

        ! master sends total interface points coordinates vector 
        ! to each slave
        call parari('BCT',0_ip,nboun_lev,lebsu_lev)

     end select

  end if


end subroutine lev_parall
