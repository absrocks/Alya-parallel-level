subroutine hlm_fivecs()

  !-----------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_fivecs.f90
  ! NAME 
  !    hlm_fivecs
  ! DESCRIPTION
  !    This routine calculates values of field vectors at given sites  
  !    from FE-computed potentials using numerical differentiation 
  !    based on the Moving Least-Squares Interpolation (MLSI) sheme.
  !    Also, it outputs values of field vectors at given sites in 
  !    an output file.
  ! USES
  !    hlm_mlsint
  !    hlm_fivout
  ! USED BY
  !    hlm_output
  !-----------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain
  use def_helmoz
  use mod_postpr
  use mod_iofile

  implicit none

  integer(ip) :: nn,ii
  real(rp)    :: invmu
  real(rp)    :: tpoin(3)
  complex(rp) :: omg,Ax,Ay,Az,P
  complex(rp) :: Asx(4),Asy(4),Asz(4),Psis(4)

  nn = nmlsi_hlm                      !Number of closest mesh nodes to a site needed for the MLSI
  invmu = 1.0_rp / perma_hlm(1)       !1/mu0
  omg = cmplx(0.0_rp, anguf_hlm, kind=rp)      ! JELENA: -i*w

  do ii = 1,nsite_hlm

        !Test point
        tpoin(1) = site_hlm(1,ii)
        tpoin(2) = site_hlm(2,ii)
        tpoin(3) = site_hlm(3,ii)

  	call hlm_mlsint(tpoin,nn,ii,Asx,Asy,Asz,Psis)       !The Moving Least-Squares Interpolation

        !Calculate secondary magnetic field intensity vector, Hs = 1/mu0 * rot(As)
        magfi_hlm(1,ii) = invmu * (Asz(2) - Asy(3))
        magfi_hlm(2,ii) = invmu * (Asx(3) - Asz(1))
        magfi_hlm(3,ii) = invmu * (Asy(1) - Asx(2))

        !Calculate secondary electric field intensity vector, Es = i * w * (As + grad(Psis))
        elefi_hlm(1,ii) = omg * (Asx(1) * tpoin(1) + Asx(2) * tpoin(2) + Asx(3) * tpoin(3) + Asx(4) + Psis(1))
        elefi_hlm(2,ii) = omg * (Asy(1) * tpoin(1) + Asy(2) * tpoin(2) + Asy(3) * tpoin(3) + Asy(4) + Psis(2))
        elefi_hlm(3,ii) = omg * (Asz(1) * tpoin(1) + Asz(2) * tpoin(2) + Asz(3) * tpoin(3) + Asz(4) + Psis(3))

  end do

  !Output values of field vectors at given sites in an output file 
  call hlm_fivout()
  
end subroutine hlm_fivecs
