subroutine tem_updbcs(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_updbcs
  ! NAME 
  !    tem_updbcs
  ! DESCRIPTION
  !    This routine updates the temperature boundary conditions:
  !    1. Before a time step begins
  !    2. Before a global iteration begins
  !    3. Before an inner iteration begins
  ! USED BY
  !    tem_begste
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame, only       :  pi
  use def_domain
  use def_kermod
  use def_temper
  use mod_ker_space_time_function
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: jtask
  integer(ip)             :: ipoin,ibopo,iboun,ipnat,ifunc,ivalu
  integer(ip)             :: kk,idime,icomp
  real(rp)                :: xx,rad
  real(rp)                :: tenew,venor,teold,cploc(6,2),dummr
  real(rp), external      :: funcre
  real(rp)                :: dummy

!  if( itask == -1 ) then
!        !  
!        ! Before a time step
!        !    
!        if( 1==1 ) then
!           icomp = min(3,ncomp_tem)
!           do ipoin = 1,npoin
!              rad = 0.0_rp
!              do idime = 1,ndime
!                 rad = rad + (coord(idime,ipoin)-0.005_rp)**2
!              end do
!              rad = sqrt(rad)
!              if( rad<0.0005_rp) then
!                 therm(ipoin,ncomp_tem) = 1.0_rp
!              end if
!           end do
!        end if
!        return
!     end if


  jtask = abs(itask)

  if( INOTMASTER ) then

     select case(jtask)

     case(1)
        !  
        ! Before a time step
        !    
        if( kfl_intbc_tem /= 0 ) then
           call tem_intbcs()
        end if

        if( kfl_conbc_tem == 0 ) then
           !
           ! Non-constant bc
           !
           if( kfl_discr_tem == NODAL_SCHEME ) then
              do ipoin = 1,npoin 
                 if(kfl_fixno_tem(1,ipoin)==1) then
                    if( kfl_funno_tem(ipoin) > 0 .and. kfl_funno_tem(ipoin) <= interval_funno ) then
                       tenew = bvess_tem(1,ipoin,2)&
                            *funcre(funpa_tem(1,kfl_funno_tem(ipoin)),6,&
                            kfl_funty_tem(kfl_funno_tem(ipoin)),cutim)
                       if( kfl_timei_tem /= 0 .and. kfl_tiacc_tem == 2 .and. kfl_tisch_tem == 1 ) then
                          bvess_tem(1,ipoin,1) = 0.50_rp*tenew+0.50_rp*therm(ipoin,ncomp_tem)
                          therm(ipoin,1)       = bvess_tem(1,ipoin,1)
                          therm(ipoin,2)       = bvess_tem(1,ipoin,1)
                       else
                          bvess_tem(1,ipoin,1) = tenew                 
                          therm(ipoin,1)       = bvess_tem(1,ipoin,1)
                          therm(ipoin,2)       = bvess_tem(1,ipoin,1)
                       end if
                    end if
                 end if
              end do
           end if

           do iboun = 1,nboun
              if( kfl_fixbo_tem(iboun) /= 0 ) then
                 if( kfl_funbo_tem(iboun) > 0 ) then
                    do ipnat = 1,npnat_tem
                       tenew = bvnat_tem(ipnat,iboun,2)&
                            *funcre(funpa_tem(1,kfl_funbo_tem(iboun)),6,&
                            kfl_funty_tem(kfl_funbo_tem(iboun)),cutim)
                       if( kfl_timei_tem /= 0 .and. kfl_tiacc_tem == 2 .and. kfl_tisch_tem == 1 ) then
                          bvnat_tem(ipnat,iboun,1) = 0.50_rp*tenew+0.50_rp*bvnat_tem(ipnat,iboun,1)
                       else
                          bvnat_tem(ipnat,iboun,1) = tenew                 
                       end if
                    end do
                 end if
              end if
           end do

           if( number_space_time_function > 0 ) then
              do ipoin = 1,npoin
                 if( kfl_funno_tem(ipoin) < 0 ) then 
                    ifunc = -kfl_funno_tem(ipoin)            
                    call ker_space_time_function(&
                         ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,tenew)
                    !                    tenew = tenew * bvess_tem(1,ipoin,2)
                    if( kfl_timei_tem /= 0 .and. kfl_tiacc_tem == 2 .and. kfl_tisch_tem == 1 ) then
                       teold = therm(ipoin,ncomp_tem)
                       bvess_tem(1,ipoin,1) = 0.50_rp*(tenew+teold)
                    else                               
                       bvess_tem(1,ipoin,1) = tenew
                    end if
                 end if
              end do
           end if

           if( kexist_tran_fiel > 0  ) then
              do ipoin=1,npoin
                 if( kfl_funno_tem(ipoin) > interval_funno ) then ! interval_funno =1000 
                    ifunc = kfl_funno_tem(ipoin) - interval_funno
                    kk = k_tran_fiel(ifunc) !indicates to which interval the current time belongs.
                    xx = x_tran_fiel(ifunc) !indicates the position between the begining and end of the interval. 
                    tenew = xfiel(ifunc) % a(1,ipoin,kk) * xx + xfiel(ifunc) % a(1,ipoin,kk+1) * (1.0_rp-xx)
                    !
                    ! These lines are identical to the ones from space time function perhaps it would be nicer to create a small subroutine  
                    !
                    if( kfl_timei_tem /= 0 .and. kfl_tiacc_tem == 2 .and. kfl_tisch_tem == 1 ) then
                       teold = therm(ipoin,ncomp_tem)
                       bvess_tem(1,ipoin,1) = 0.50_rp * ( tenew + teold )
                    else
                       bvess_tem(1,ipoin,1) = tenew
                    end if
                 end if
              end do
           end if

           if( kfl_exist_fixi7_tem == 1 .and. itask==1 ) then
              !
              ! Open or closes point depending on angle between normal and VELOC
              !
              call memgen(1_ip,npoin,0_ip)   !allocate gisca
              call open_close(2,kfl_fixno_tem,dummy,1_ip) ! calculates gisca(1:npoin)  ! The last parameter (1_ip) corresponds to the first dimension of kfl_fixno_tem
              do ipoin=1,npoin
                 if(abs(kfl_fixno_tem(1,ipoin)) == 7 ) then
                    if( gisca(ipoin) > 0 ) then  
                       kfl_fixno_tem(1,ipoin) = 7
                    else
                       kfl_fixno_tem(1,ipoin) = -7
                    end if
                 end if
              end do
              call memgen(3_ip,npoin,0_ip)   !deallocate gisca
           end if

        end if

        if (kfl_regim_tem==4 .and. kfl_plepp_tem /= 4) call tem_calcEnthalpyBC()
        !
        ! Impose Dirichlet bc
        !
        do ipoin = 1,npoin
           if( kfl_fixno_tem(1,ipoin) > 0 ) then
              tempe(ipoin,2) = bvess_tem(1,ipoin,1)
           end if
        end do

     case(2)
        !
        ! Before a global iteration
        !  
        do ipoin=1,npoin
           if(kfl_fixno_tem(1,ipoin)==4.or.kfl_fixno_tem(1,ipoin)==-4) then
              ibopo=lpoty(ipoin)
              if(ibopo>0) then
                 venor=dot_product(veloc(1:ndime,ipoin,1),exnor(1:ndime,1,ibopo))
                 if(venor<=0.0_rp) then
                    kfl_fixno_tem(1,ipoin)= 4
                 else
                    kfl_fixno_tem(1,ipoin)=-4
                 end if
              else
                 kfl_fixno_tem(1,ipoin)=-4
              end if
           end if
        end do

        if (kfl_regim_tem==4 .and. kfl_plepp_tem /= 4) then
           do ipoin=1,npoin
              if(kfl_fixno_tem(1,ipoin)==1) then
                 do ivalu = 1,6
                    cploc(ivalu,1) = sphec(ipoin,ivalu,1)
                    cploc(ivalu,2) = sphec(ipoin,ivalu,2)
                 end do
                 dummr = 0.0_rp
                 call tem_comput(2_ip,bvess_tem(1,ipoin,2),dummr,cploc,tenew)
                 bvess_tem(1,ipoin,1) = tenew
                 therm(ipoin,1)       = bvess_tem(1,ipoin,1)
                 therm(ipoin,2)       = bvess_tem(1,ipoin,1)
              end if
           end do
        end if


     case (3)
        !
        ! Before an inner iteration
        ! 

     end select


  end if

end subroutine tem_updbcs
