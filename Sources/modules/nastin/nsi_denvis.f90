subroutine nsi_denvis()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_denvis
  ! NAME 
  !    nsi_denvis
  ! DESCRIPTION
  !    Smoothing of density and viscosity
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    nsi_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_proper 
  implicit none
  integer(ip)       :: ipoin,ndofn,dummi
  real(rp)          :: xfact
  real(rp), pointer :: prope_tmp(:) 

  if( INOTMASTER ) then 

     nullify ( prope_tmp )
     allocate( prope_tmp(npoin) )
     call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)
     do ipoin = 1,npoin
        prope_nsi(1,ipoin) = prope_tmp(ipoin)
     end do
     call ker_proper('VISCO','NPOIN',dummi,dummi,prope_tmp)
     do ipoin = 1,npoin
        prope_nsi(2,ipoin) = prope_tmp(ipoin)
     end do
     deallocate( prope_tmp )
     
  end if

end subroutine nsi_denvis
