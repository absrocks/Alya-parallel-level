subroutine wav_matrix
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_matrix
  ! NAME 
  !    wav_matrix
  ! DESCRIPTION
  !    Compute matrix and RHS
  ! USES
  !    wav_elmope
  ! USED BY
  !    wav_solite
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_wavequ
  use      def_domain
  use      def_solver
  implicit none
  integer(ip), save :: ipass=0
  integer(ip)       :: izmat,ipoin
  real(rp)          :: dt2
  !
  ! LHS assembly
  !
  if(ipass==0) then

     ipass=1
     do izmat=1,nzmat
        amatr(izmat)=0.0_rp
     end do
     call wav_elmope(1_ip)
     call wav_bouope(1_ip)

     if(kfl_tisch_wav==3.and.kfl_massm_wav==0) then
        !
        ! Newmark: beta*A+M/dt^2
        !
        dt2=dtinv_wav*dtinv_wav
        call wav_addmas(npoin,amatr,vmass,dt2,c_sol,r_sol)        
     end if

     if(kfl_absor_wav==1) then
        !
        ! Exchange boundary terms of boundary mass matrix
        !
        call wav_parall(5_ip)
        !
        ! Add original mass matrix
        !
        do ipoin=1,npoin
           vmass_wav(ipoin)=vmass_wav(ipoin)+vmass(ipoin)
        end do
     end if

  end if
  !
  ! RHS assembly
  !
  do ipoin=1,npoin
     rhsid(ipoin)=0.0_rp
  end do
  call wav_elmope(2_ip)
  call wav_bouope(2_ip)

  if(kfl_tisch_wav==3.and.kfl_massm_wav==0) then
     !
     ! Newmark: RHS <= RHS + M/dt^2*u^n + M/dt*v^n + (0.5-beta)*M*a^n
     !
     dt2=dtinv_wav*dtinv_wav
     do ipoin=1,npoin
        rhsid(ipoin)=rhsid(ipoin)+vmass(ipoin)&
             *(dt2*wavam(ipoin,3)+dtinv_wav*wavve_wav(ipoin,2)+(0.5_rp-nebet_wav)*wavac_wav(ipoin,2))
     end do
  end if
  
end subroutine wav_matrix
