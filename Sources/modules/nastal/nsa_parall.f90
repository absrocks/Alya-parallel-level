subroutine nsa_parall(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsa_parall
  ! NAME
  !    nsa_parall
  ! DESCRIPTION
  !    This routine is a bridge to Parall service  
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  use      def_parall
  use mod_communications, only       :  PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only       :  PAR_BARRIER
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipart,istar,istop

  call nsa_chkpar  

  if ( weparal) then

     select case(itask)

     case(1_ip)
        !
        ! Exchange data read in nsa_reaphy, nsa_reanut and nsa_reaous
        ! always using MPI, even if this is a partition restart
        !
        call nsa_sendat(one)

     case(2_ip)    !! DMM This has to be removed, not used anymore!
        !
        ! Exchange data read in nsa_reabcs
        !
        pard1=nproc_par-1                                   ! Bridge
        pard2=nsire_par
        istar = pard1
        istop = pard2
        do ipart = 0,istar,istop
           if( kfl_paral >= ipart .and. kfl_paral < ipart+istop ) then
              call nsa_sendat(two)
           end if
           call PAR_BARRIER()
        end do

     case(3_ip)
        !
        ! Sum up residual contribution of slave neighbors
        !
        if ( islave ) then
           call vocabu(NPOIN_REAL_12DI,ndime,0_ip)
           parr1 => rhsid(1:ndime*npoin)
           call par_slexch()
        end if

     case(4_ip)
        !
        ! Sum up residual contribution of slave neighbors
        !
        if ( islave ) then
           call vocabu(NPOIN_REAL_12DI,ndofn_nsa,0_ip)
           parr1 => rhsid(1:ndofn_nsa*npoin)
           call par_slexch()
        end if

     case(5_ip)
        !
        ! Sum up vorticity and gradient invariants contribution of slave neighbors
        !
        if ( islave ) then
           call vocabu(NPOIN_REAL_2DIM,ndime+1,0_ip)
           parr2 => vorti(1:ndime+1,1:npoin)
           call par_slexch()
        end if

     case(6_ip)
        !
        ! Sum up aerodynamic body forces
        !
        nparr =  20         ! total area + frontal areas(3) + press forces(3) + visco forces(3) + free space
        parre => ccoef_nsa
        call par_operat(3_ip)                     
        
     case(7_ip)
        !
        ! Sum up subscale contribution of slave neighbors
        !
        if ( islave ) then
           call rhsmod(ndime,umoss_nsa(1,1,2))
           call rhsmod(1_ip,denss_nsa(1,2))
           call rhsmod(1_ip,eness_nsa(1,2))
           call rhsmod(1_ip,frequ_nsa)
        end if

     case(8_ip)
        !
        ! Sum up nodal time steps
        !
        if ( islave ) then
           call vocabu(NPOIN_REAL_2DIM,3,0_ip)
           parr2 => dtieq_nsa(1:3,1:npoin,1)
           call par_slexch()
        end if

     case(9_ip)
        !
        ! Interchange corrected fixnos
        !
        if ( islave ) then
!           call PAR_INTERFACE_NODE_EXCHANGE(ndofn_nsa,kfl_fixno_nsa,'MAX','IN MY CODE','SYNCHRONOUS')
           call vocabu(NPOIN_INTE_2DIM,ndofn_nsa,0_ip)
           pari2 => kfl_fixno_nsa
           call par_slexch()
        end if

     end select

  end if

end subroutine nsa_parall
