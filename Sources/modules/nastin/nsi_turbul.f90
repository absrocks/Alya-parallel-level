subroutine nsi_turbul(&
     itask,jtask,pnode,pgaus,igaui,igauf,&
     gpsha,gpcar,elvel,gpden,gpvis,gpmut, &
     gpgvi,grvis,gpgve,ielem,kfl_kxmod)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_turbul
  ! NAME 
  !    nsi_turbul
  ! DESCRIPTION
  !    Compute viscosity and its gradient due to turbulence
  !    GBMUT ....... mut
  !    JTASK = 1 ... Compute gradient GRVIS = grad(mut)
  !            0 ... Do not compute gradient GRVIS = grad(mut)
  !    ITASK = 1 ... GPVIS <= GPVIS + GPMUT
  !            0 ... Do not change GPVIS
  !    GPGVI = Viscosity gradient grad(mu+mut)
  ! USES
  !
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,nvart  
  use def_domain, only       :  mnode,ndime, lnods, walld
  use def_master, only       :  zeror, untur, postp
  use def_nastin, only       :  turbu_nsi,zensi,turmu_nsi,kfl_grvir_nsi
  use def_kermod, only       :  turmu_ker, kfl_logva
  use mod_ker_regularization, only : regul_k, regul_e, kfl_regularization
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  implicit none 
  integer(ip),  intent(in)   :: itask,jtask,pnode,pgaus,igaui,igauf
  integer(ip),  intent(in)   :: kfl_kxmod, ielem
  real(rp),     intent(in)   :: gpsha(pnode,pgaus)
  real(rp),     intent(in)   :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)   :: elvel(ndime,pnode)
  real(rp),     intent(in)   :: gpden(pgaus)
  real(rp),     intent(out)  :: gpvis(pgaus)
  real(rp),     intent(inout):: gpmut(pgaus)
  real(rp),     intent(inout):: gpgvi(ndime,pgaus)
  real(rp),     intent(out)  :: grvis(ndime,pgaus)
  real(rp),     intent(out)  :: gpgve(ndime,ndime,pgaus)
  integer(ip)                :: igaus,inode,idime,jdime,kdime, ipoin
  real(rp)                   :: uaste,eltur(2, mnode),seci4
  real(rp)                   :: A0, As, W, phi, cmu, divve, simgr(3,3), kinen 
  real(rp)                   :: epsil, gpwal, omega, gpvor, a1, F2, Cr, f0, f02, Wsqr6
  real(rp)                   :: regue, reguk, sigmr


  if(turmu_ker % kfl_exist /= 0_ip) then ! matias suggested I use this - seems good idea 
     
     !
     ! Compute mu_t = rho * nu_t for LES
     !
     gpmut(igaui:igauf) = gpden(igaui:igauf) * gpmut(igaui:igauf)    ! DEJARLO CONSISTENTE   en element operations esta fuera!!!!!!
     !
     ! Compute effective viscosity & gradient of viscosity: 
     !    mu_eff        = mu + mu_t 
     !    grad (mu_eff) = grad(mu) + grad(mu_t)
     !
     if (itask == 1) then
        do igaus = igaui,igauf
           gpvis(igaus) = gpvis(igaus) + gpmut(igaus)
           do idime = 1,ndime
              gpgvi(idime,igaus) = gpgvi(idime,igaus) + grvis(idime,igaus)
           end do
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! Needed for postprocessing turbulent viscosity
     !
     !----------------------------------------------------------------------
     if ( jtask == 1 ) then
        if(         output_postprocess_check_variable_postprocess(74_ip)  &
             & .or. output_postprocess_check_variable_postprocess(75_ip)  &
             & .or. output_postprocess_check_variable_postprocess(76_ip)  &
             & .or. output_postprocess_check_variable_postprocess(77_ip)  &
             & .or. output_postprocess_check_variable_postprocess(78_ip)  &
             & .or. output_postprocess_check_variable_postprocess(92_ip)  &
             & .or. output_postprocess_check_variable_postprocess(93_ip)  &
             & .or. output_postprocess_check_variable_postprocess(94_ip)  &
             & .or. output_postprocess_check_variable_postprocess(95_ip)  ) then
           do igaus = igaui,igauf
              turmu_nsi(ielem)%a(1,igaus,1) = gpmut(igaus)
           end do
        end if
     end if

  end if

end subroutine nsi_turbul
