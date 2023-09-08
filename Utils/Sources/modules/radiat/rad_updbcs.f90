subroutine rad_updbcs(itask)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_updbcs
  ! NAME 
  !    rad_updbcs
  ! DESCRIPTION
  !    This routine updates the radiation boundary conditions:
  !    1. Before a time step begins
  !    2. Before a global iteration begins
  !    3. Before an inner iteration begins
  ! USED BY
  !    rad_begste
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_radiat
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ibopo,iboun,ipnat
  real(rp)                :: tenew,venor
  real(rp), external      :: funcre

  if(kfl_paral/=0) then

     select case(itask)

     case(1)

        if(kfl_exacs_rad/=0) then
           !
           ! Exact solution
           !
           call rad_exabcs(2_ip)

!!$        else
!!$           !  
!!$           ! Before a time step
!!$           !     
!!$           if(kfl_intbc_rad/=0) then
!!$              call rad_intbcs()
!!$           end if
!!$
!!$           if(kfl_conbc_rad==0) then
!!$              !
!!$              ! Non-constant bc
!!$              !
!!$              do ipoin=1,npoin 
!!$                 if(kfl_fixno_rad(1,ipoin)==1) then
!!$                    if(kfl_funno_rad(ipoin)/=0) then
!!$                       genew=bvess_rad(ipoin,2)&
!!$                            *funcre(funpa_rad(1,kfl_funno_rad(ipoin)),6,&
!!$                            kfl_funty_rad(kfl_funno_rad(ipoin)),cutim)
!!$                       
!!$                       bvess_rad(ipoin,1) = genew                 
!!$                       radav_rad(ipoin,1)     = bvess_rad(ipoin,1)
!!$                       radav_rad(ipoin,2)     = bvess_rad(ipoin,1)
!!$                    end if
!!$                 end if
!!$              end do
!!$
!!$              do iboun=1,nboun
!!$                 if(kfl_fixbo_rad(iboun)/=0) then
!!$                    if(kfl_funbo_rad(iboun)/=0) then
!!$                       do ipnat=1,npnat_rad
!!$                          genew=bvnat_rad(ipnat,iboun,2)&
!!$                               *funcre(funpa_rad(1,kfl_funbo_rad(iboun)),6,&
!!$                               kfl_funty_rad(kfl_funbo_rad(iboun)),cutim)
!!$                          bvnat_rad(ipnat,iboun,1)=genew                 
!!$                       end do
!!$                    end if
!!$                 end if
!!$              end do
!!$
!!$           end if

        end if

     case(2)
        !
        ! Before a global iteration
        !  
!!$        do ipoin=1,npoin
!!$           if(kfl_fixno_rad(1,ipoin)==4.or.kfl_fixno_rad(1,ipoin)==-4) then
!!$              ibopo=lpoty(ipoin)
!!$              if(ibopo>0) then
!!$                 venor=dot_product(veloc(1:ndime,ipoin,1),exnor(1:ndime,1,ibopo))
!!$                 if(venor<=0.0_rp) then
!!$                    kfl_fixno_rad(1,ipoin)= 4
!!$                 else
!!$                    kfl_fixno_rad(1,ipoin)=-4
!!$                 end if
!!$              else
!!$                 kfl_fixno_rad(1,ipoin)=-4
!!$              end if
!!$           end if
!!$        end do

     case (3)
        !
        ! Before an inner iteration
        ! 

     end select

  end if

end subroutine rad_updbcs
