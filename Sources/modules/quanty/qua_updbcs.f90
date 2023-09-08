subroutine qua_updbcs(itask)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_updbcs
  ! NAME 
  !    qua_updbcs
  ! DESCRIPTION
  !    This routine updates the Schrodinger boundary conditions:
  !    1. Before a time step begins
  !    2. Before a global iteration begins
  !    3. Before an inner iteration begins
  ! USES
  !  qua_exabcs
  !  qua_intbcs
  !  
  ! USED BY
  !    qua_begste
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_quanty
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ibopo,iboun,ipnat
  real(rp)                :: tenew,venor
  real(rp), external      :: funcre

  if(kfl_paral/=0) then

     select case(itask)

     case(1)

        if(kfl_exacs_qua/=0) then
           !
           ! Exact solution
           !
           !call qua_exabcs(2_ip)

        else
           !  
           ! Before a time step
           !     
           if(kfl_intbc_qua/=0) then
          !    call qua_intbcs()
           end if

           if(kfl_conbc_qua==0) then
              !
              ! Non-constant bc
              !
              do ipoin=1,npoin 
                 if(kfl_fixno_qua(1,ipoin)==1) then
                    if(kfl_funno_qua(ipoin)/=0) then
                       tenew=bvess_qua(ipoin,2)&
                            *funcre(funpa_qua(1,kfl_funno_qua(ipoin)),6,&
                            kfl_funty_qua(kfl_funno_qua(ipoin)),cutim)
                       if(kfl_timei_qua/=0.and.kfl_tiacc_qua==2.and.kfl_tisch_qua==1) then
                          bvess_qua(ipoin,1)=0.50_rp*tenew+0.50_rp*phion(ipoin,ncomp_qua)
                       else
                          bvess_qua(ipoin,1)=tenew                 
                       end if
                    end if
                 end if
              end do


           end if

        end if

     case(2)
        !
        ! Before a global iteration
        !  
        do ipoin=1,npoin
           if(kfl_fixno_qua(1,ipoin)==4.or.kfl_fixno_qua(1,ipoin)==-4) then
              ibopo=lpoty(ipoin)
              if(ibopo>0) then
                 venor=dot_product(veloc(1:ndime,ipoin,1),exnor(1:ndime,1,ibopo))
                 if(venor<=0.0_rp) then
                    kfl_fixno_qua(1,ipoin)= 4
                 else
                    kfl_fixno_qua(1,ipoin)=-4
                 end if
              else
                 kfl_fixno_qua(1,ipoin)=-4
              end if
           end if
        end do

     case (3)
        !
        ! Before an inner iteration
        ! 

     end select

  end if

end subroutine qua_updbcs
