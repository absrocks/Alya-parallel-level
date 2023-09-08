subroutine wav_explic
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_explic
  ! NAME 
  !    wav_explic
  ! DESCRIPTION
  !    Explicit advance
  !    g0*Mu^{n+1} + (A+g1*M)u^{n} + g2*Mu^{n-1} = b^n
  !    <=>
  !    u^{n+1}=1/g0*[-g1*u^{n}-g2*u^{n-1}+M^{-1}*[b^n-Au^{n)]]
  !
  ! USES
  ! USED BY
  !    wav_solite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_wavequ
  implicit none
  integer(ip) :: ipoin
  real(rp)    :: gamm0,gamm1,gamm2,dtn,dtn1,fact1,fact2
  !
  ! Au^{n}
  !
  call bcsrax(1_ip,npoin,1_ip,amatr,c_sol,r_sol,wavam(1,3),unkno)
  !
  ! Service Parall
  !
  call wav_parall(3_ip) ! Au^{n}
  call wav_parall(4_ip) ! RHS
  !
  ! Advance in time
  !
  if(kfl_tisch_wav==2) then
     dtn   =  1.0_rp/dtinv_wav 
     dtn1  =  dtold(1)
     gamm0 =  2.0_rp/(dtn*(dtn+dtn1))
     gamm1 = -2.0_rp/(dtn*dtn1)
     gamm2 =  2.0_rp/(dtn1*(dtn+dtn1))
     fact1 =  1.0_rp/gamm0
     if(kfl_absor_wav==0) then
        do ipoin=1,npoin
           unkno(ipoin)=fact1&
                &       *(-gamm1*wavam(ipoin,3)&
                &         -gamm2*wavam(ipoin,4)&
                &         +(rhsid(ipoin)-unkno(ipoin))/vmass(ipoin))
        end do
     else
        do ipoin=1,npoin
           unkno(ipoin)=fact1/vmass_wav(ipoin)&
                &       *(-gamm1*vmass(ipoin)*wavam(ipoin,3)&
                &         -gamm2*vmass(ipoin)*wavam(ipoin,4)&
                &         +rhsid(ipoin)-unkno(ipoin))
        end do
     end if
  end if
  !
  ! Impose boundary conditions
  !
  if(kfl_onnod_wav==1) then
     do ipoin=1,npoin
        if(kfl_fixno_wav(ipoin)==1) unkno(ipoin)=bvess_wav(ipoin)
     end do
  end if
  !
  ! Update subgrid scale
  !
  if(kfl_subgs_wav==1) call wav_elmope(3_ip)

end subroutine wav_explic
