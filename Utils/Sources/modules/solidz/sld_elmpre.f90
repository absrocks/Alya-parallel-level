subroutine sld_elmpre(&
     ielem,pnode,pgaus,pmate,eldis,eldip,elvel,elacc,eldix,elvex,elcod,&
     elmof,elrst,gpsha,gpcar,gphea,gpgi0,gpdet,gpvol,&
     elepo,elepo_new_inverse,elepo_new_determinant,&
     gpvel,gpacc,&
     gprat,gpdis,gpcau,gpigd,gpcod,gpgdi,&
     gpigd_eps,gpgdi_eps,gpdet_eps,&
     gpmof,gprestre,gpdet_min,&
     ielem_min,trpn_sld_temp)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_elmpre
  ! NAME
  !    sld_elmpre
  ! DESCRIPTION
  !    Compute some Gauss values
  ! INPUT
  !    GPGI0 ... Previous deformation gradient tensor ............ F(n)
  ! OUTPUT
  !    GPDIS ... Deformation ..................................... phi (phidot is elvel)
  !    GPVEL ... Deformation velocity............................. phidot
  !    GPGDI ... Updated deformation gradient .....................F(n+1) = grad(phi)
  !    GPIGD ... Inverse of updated deformation gradient tensor .. F^{-1}(n+1)
  !    GPRAT ... Rate of Deformation   ........................... Fdot = grad(phidot)
  !    GPCAU ... Right Cauchy-Green deformation tensor ........... C = F^t x F
  !    GPCOD ... Physical coordinates of the gauss points
  !    GPRESTRE ... Residual stresses
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp,lg
  use def_domain, only       :  mnode,ndime,mgaus,lnods, vmass
  use def_solidz, only       :  volum_sld
  use def_solidz, only       :  kfl_xfeme_sld, kfl_ninex_sld
  use def_solidz, only       :  kfl_gdepo,covar_sld,kfl_coupt_sld
  use def_solidz, only       :  kfl_moduf_sld,nvgij_inv_sld,kfl_restr_sld,kfl_prest_sld
  use def_master, only       :  kfl_paral,coupling, ittim, kfl_eccty, ITER_K, TIME_N, stateland
  implicit none
  integer(ip), intent(in)    :: ielem
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pgaus
  integer(ip), intent(in)    :: pmate
  real(rp),    intent(in)    :: eldis(ndime,pnode,*)
  real(rp),    intent(in)    :: eldip(ndime,pnode  )
  real(rp),    intent(in)    :: elvel(ndime,pnode)
  real(rp),    intent(in)    :: elacc(ndime,pnode,3)
  real(rp),    intent(in)    :: elrst(    6,pnode)
  real(rp),    intent(in)    :: eldix(ndime,pnode,*)
  real(rp),    intent(in)    :: elvex(ndime,pnode)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(in)    :: elmof(pnode)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gphea(pnode,pgaus)
  real(rp),    intent(in)    :: gpgi0(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpdet(pgaus)
  real(rp),    intent(in)    :: gpdet_eps(pgaus)
  real(rp),    intent(inout) :: gpvol(pgaus)
  real(rp),    intent(out)   :: elepo(ndime,ndime,pnode)
  real(rp),    intent(inout) :: trpn_sld_temp(pnode)
  real(rp),    intent(out)   :: elepo_new_inverse(ndime,ndime,pnode)
  real(rp),    intent(out)   :: elepo_new_determinant(pnode)
  real(rp),    intent(out)   :: gpvel(ndime,pgaus,*)
  real(rp),    intent(out)   :: gpacc(ndime,pgaus,3)
  real(rp),    intent(out)   :: gprat(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpdis(ndime,pgaus,*)       !< Deformation: phi
  real(rp),    intent(out)   :: gpcau(ndime,ndime,pgaus)   !< Right Cauchy-Green deformation tensor: C = F^t x F
  real(rp),    intent(out)   :: gpcod(ndime,mgaus)
  real(rp),    intent(out)   :: gpigd(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpigd_eps(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpgdi(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpgdi_eps(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpmof(pgaus)
  real(rp),    intent(out)   :: gprestre(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpdet_min
  integer(ip), intent(out)   :: ielem_min
  real(rp)                   :: detme,volel,bidon,xmean
  integer(ip)                :: inode,igaus,idime,jdime
  integer(ip)                :: kdime,ipoin,kauxi,ivoig
  real(rp)                   :: gpig0(ndime,ndime), gpgd1(ndime,ndime)
  real(rp)                   :: gpde0
  !
  ! GPDIS and GPVEL: displacement and velocity
  !
  kauxi= 0
  gpacc = 0.0_rp
  elepo = 0.0_rp
  elepo_new_determinant = 0.0_rp
  elepo_new_inverse = 0.0_rp

  do idime = 1,ndime
     do igaus = 1,pgaus
        gpcod(idime,igaus)   = 0.0_rp
        gpdis(idime,igaus,ITER_K) = 0.0_rp
        gpdis(idime,igaus,TIME_N) = 0.0_rp
        gpvel(idime,igaus,ITER_K) = 0.0_rp
        do inode = 1,pnode
           gpdis(idime,igaus,ITER_K) = gpdis(idime,igaus,ITER_K) + eldis(idime,inode,ITER_K)*gpsha(inode,igaus)
           gpdis(idime,igaus,TIME_N) = gpdis(idime,igaus,TIME_N) + eldis(idime,inode,TIME_N)*gpsha(inode,igaus)
           gpvel(idime,igaus,ITER_K) = gpvel(idime,igaus,ITER_K) + elvel(idime,inode)*gpsha(inode,igaus)
           gpacc(idime,igaus,ITER_K) = gpacc(idime,igaus,ITER_K) + elacc(idime,inode,ITER_K)*gpsha(inode,igaus)
           gpacc(idime,igaus,TIME_N) = gpacc(idime,igaus,TIME_N) + elacc(idime,inode,TIME_N)*gpsha(inode,igaus)
           gpcod(idime,igaus)   = gpcod(idime,igaus)   + elcod(idime,inode)*gpsha(inode,igaus)
        end do
        do jdime = 1,ndime
           gprestre(idime,jdime,igaus) = 0.0_rp
           gpgdi(idime,jdime,igaus) = 0.0_rp
           gpgdi_eps(idime,jdime,igaus) = 0.0_rp
           gprat(idime,jdime,igaus) = 0.0_rp
        end do
     end do
  end do
  !
  ! GPDIS and GPVEL: displacement and velocity
  ! to include the contribution of the enrichment dofs
  !
  if (kfl_xfeme_sld == 1) then
     do idime=1,ndime
        do igaus=1,pgaus
           do inode=1,pnode
              gpdis(idime,igaus,ITER_K)=gpdis(idime,igaus,ITER_K)&
                   +eldix(idime,inode,ITER_K)*gpsha(inode,igaus)*gphea(inode,igaus)
              gpdis(idime,igaus,TIME_N)=gpdis(idime,igaus,TIME_N)&
                   +eldix(idime,inode,TIME_N)*gpsha(inode,igaus)*gphea(inode,igaus)
              gpvel(idime,igaus,ITER_K)=gpvel(idime,igaus,ITER_K)&
                   +elvex(idime,inode)*gpsha(inode,igaus)*gphea(inode,igaus)
           end do
        end do
     end do
  end if

  !
  ! GPGDI: Deformation gradients
  !
  if (kfl_prest_sld == 0) then
     do igaus=1,pgaus
        do idime=1,ndime
           do jdime=1,ndime
              do inode=1,pnode
                 gpgdi(idime,jdime,igaus)=gpgdi(idime,jdime,igaus)&
                      +eldis(idime,inode,ITER_K)*gpcar(jdime,inode,igaus)
                 gprat(idime,jdime,igaus)=gprat(idime,jdime,igaus)&
                      +elvel(idime,inode)*gpcar(jdime,inode,igaus)
                 gpgdi_eps(idime,jdime,igaus)=gpgdi_eps(idime,jdime,igaus)&
                      +(eldis(idime,inode,ITER_K)+eldip(idime,inode))*gpcar(jdime,inode,igaus)
              end do
           end do
        end do
     end do
  end if

  !
  !PRESTRESS
  !
  if (kfl_prest_sld == 1) then
     do igaus=1,pgaus
        if (ittim < 3) then
           do idime=1,ndime
              do jdime=1,ndime
                 do inode=1,pnode
                    gpgdi(idime,jdime,igaus)=gpgdi(idime,jdime,igaus)&
                         +eldis(idime,inode,ITER_K)*gpcar(jdime,inode,igaus)
                 end do
              end do
           end do
           gpgdi(1,1,igaus)= gpgdi(1,1,igaus) + 1.0_rp
           gpgdi(2,2,igaus)= gpgdi(2,2,igaus) + 1.0_rp
           if (ndime==3) gpgdi(3,3,igaus)= gpgdi(3,3,igaus) + 1.0_rp
        else
           call invmtx(gpgi0(1,1,igaus),gpig0(1,1),gpde0,ndime)
           do idime=1,ndime
              do jdime=1,ndime
                 gpgd1(idime,jdime) =0.0_rp
                 do inode=1,pnode
                    gpgd1(idime,jdime)=gpgd1(idime,jdime)&
                         +eldis(idime,inode,ITER_K)*gpcar(jdime,inode,igaus)/gpde0
                 end do
              end do
           end do
           gpgd1(1,1)= gpgd1(1,1) + 1.0_rp
           gpgd1(2,2)= gpgd1(2,2) + 1.0_rp
           if (ndime==3) gpgd1(3,3)= gpgd1(3,3) + 1.0_rp
           do idime=1,ndime
              do jdime=1,ndime
                  do kdime=1,ndime
                     gpgdi(idime,jdime,igaus) = gpgdi(idime,jdime,igaus) + gpgd1(idime,kdime)*gpgi0(kdime,jdime,igaus)
                  end do
              end do
           end do
        end if
     end do
  endif

  !
  ! GPRESTRE: Residual stresses
  !
  if (kfl_restr_sld < 0) then
     do igaus=1,pgaus
        do idime=1,ndime
           do jdime=1,ndime
              do inode=1,pnode
                 ivoig= nvgij_inv_sld(idime,jdime)
                 gprestre(idime,jdime,igaus)=gprestre(idime,jdime,igaus)&
                      +elrst(ivoig,inode)*gpsha(inode,igaus)
              end do
           end do
        end do
     end do
  end if
  !
  ! GPGDI: Deformation gradients for xfem case
  ! to include the contribution of the enrichment dofs  !
  ! MARIANO: REVISAR ESTO PARA EL NINEX
  if (kfl_xfeme_sld == 1) then
     do igaus=1,pgaus
        do idime=1,ndime
           do jdime=1,ndime
              do inode=1,pnode
                 gpgdi(idime,jdime,igaus)=gpgdi(idime,jdime,igaus)&
                      +eldix(idime,inode,1)*gpcar(jdime,inode,igaus)*gphea(inode,igaus)
                 gprat(idime,jdime,igaus)=gprat(idime,jdime,igaus)&
                      +elvex(idime,inode)*gpcar(jdime,inode,igaus)*gphea(inode,igaus)
              end do
           end do
        end do
     end do
  end if
  !
  ! GPMOF: Modulating fields, elmof is set to 1 when there are no fields (then, gpmof will be 1)
  !
  if( kfl_moduf_sld(1) > 0 ) then
     do igaus = 1,pgaus
        gpmof(igaus) = 0.0_rp
        do inode = 1,pnode
           gpmof(igaus) = gpmof(igaus) + elmof(inode) * gpsha(inode,igaus)
        end do
     end do
  end if

  detme= 0.0_rp
  volel= 0.0_rp
  gpigd= 0.0_rp
  gpigd_eps= 0.0_rp

  do igaus=1,pgaus
     if(kfl_prest_sld == 0) then   !without prestress
        gpgdi(1,1,igaus)= gpgdi(1,1,igaus) + 1.0_rp
        gpgdi(2,2,igaus)= gpgdi(2,2,igaus) + 1.0_rp
        if (ndime==3) gpgdi(3,3,igaus)= gpgdi(3,3,igaus) + 1.0_rp
        gpgdi_eps(1,1,igaus)= gpgdi_eps(1,1,igaus) + 1.0_rp
        gpgdi_eps(2,2,igaus)= gpgdi_eps(2,2,igaus) + 1.0_rp
        if (ndime==3) gpgdi_eps(3,3,igaus)= gpgdi_eps(3,3,igaus) + 1.0_rp
     end if
     !Compute J
     call invmtx(gpgdi(1,1,igaus),gpigd(1,1,igaus),gpdet(igaus),ndime)
     if (kfl_ninex_sld == 1) call invmtx(gpgdi_eps(1,1,igaus),gpigd_eps(1,1,igaus),gpdet_eps(igaus),ndime)

     !if( kfl_exacs_sld == 0 ) then

     if( gpdet(igaus) < 1.0e-12_rp ) then
        !write(*,*) 'SLD_ELMPRE:  **WARNING** GPDET < 1.0E-12. Check element ', ielem,' subdomain ', kfl_paral
        !write(*,*) '      value: ',gpdet(igaus)
        !if( kfl_xfeme_sld == 1 ) then
        !   write(*,*) 'The element is enriched'
        !end if
     else if( gpdet(igaus) < 1.0e-4_rp ) then
        kauxi = 1
     end if

     if( gpdet(igaus) <= 0.0_rp .and. gpdet(igaus) < gpdet_min ) then
        !$OMP CRITICAL (compute_gpdet_min)
        if( gpdet(igaus) < gpdet_min ) then
           !
           ! GPDET < GPDET_MIN must be checked just in case another thread has just modified it!
           !
           gpdet_min = gpdet(igaus)
           ielem_min = ielem
        end if
        !$OMP END CRITICAL (compute_gpdet_min)

        !write(*,*) &
        !     '----> SLD_ELMPRE:  **ERROR** NEGATIVE GPDET. Check element ', ielem, &
        !     ' subdomain ', kfl_paral
        !write(*,*) &
        !     '      value: ',gpdet(igaus)
        !write(*,*) &
        !     '      Nodal coordinates: '
        !if (kfl_xfeme_sld==1) then
        !   write(*,*) &
        !        'The element is enriched'
        !end if
        !do inode=1,pnode
        !   write(*,*) &
        !        '      node=',inode,'  position=',elcod(1:ndime,inode)
        !end do
        !call runend ('SLD_ELMPRE: NEGATIVE JACOBIAN DETERMINANT')
     end if
     !end if
     !
     ! Witness element (debugging) connected with sld_inivar
     !
     !     if (ielem==kfl_ielem) then
!!!     if (ielem == 1167) then        ! ALE Tezduyar rigidization parameters
!!!        bidon= gpvol(igaus)
!!!        gpvol(igaus) = bidon * (1.0_rp/gpdet(igaus))**2.0_rp
!!!     end if
!!!          if (ielem==1167) then
!!!             write(4112,*) cutim,gpdet(igaus),gpvol(igaus)
!!!          end if


     volel = volel + gpvol(igaus)
     detme = detme + gpvol(igaus) * gpdet(igaus)

     !
     ! Land ECC model (coupling variable with ORd)
     !
     if ((( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 )) &
      .and. (kfl_eccty(pmate) ==  3 .or. kfl_eccty(pmate) == 4)) then
       trpn_sld_temp = 0.0_rp
       do inode = 1,pnode
         xmean = gpsha(inode,igaus) * gpvol(igaus)
         trpn_sld_temp(inode) = trpn_sld_temp(inode)+xmean*stateland(3,ielem,igaus,2)
       end do
     end if

     
     !
     ! GDEPO: Push forward
     !
     if( kfl_gdepo /= 0 ) then
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           xmean = gpsha(inode,igaus) * gpvol(igaus)
           elepo_new_determinant(inode) = elepo_new_determinant(inode) + xmean * gpdet(igaus)
           do idime= 1,ndime
              do jdime= 1,ndime
                 elepo(idime,jdime,inode) = elepo(idime,jdime,inode) + xmean * gpgi0(idime,jdime,igaus)
                 elepo_new_inverse(idime,jdime,inode) = elepo_new_inverse(idime,jdime,inode) + xmean * gpigd(idime,jdime,igaus)
              end do
           end do
        end do
     end if

  end do
  volum_sld(1) = volum_sld(1) + volel         ! reference
  volum_sld(2) = volum_sld(2) + detme         ! deformed

  !
  ! GPCAU: Cauchy tensor
  !

  do igaus=1,pgaus
     do kdime=1,ndime
        do jdime=1,ndime
           gpcau(jdime,kdime,igaus)=0.0_rp
           do idime=1,ndime
              gpcau(jdime,kdime,igaus)= gpcau(jdime,kdime,igaus) + gpgdi(idime,jdime,igaus)*gpgdi(idime,kdime,igaus)
           end do
        end do
     end do
  end do

  if (kfl_coupt_sld(pmate) == 11) then
     covar_sld(ielem)= 1.0_rp
     if (kauxi==1) covar_sld(ielem)= 0.0_rp
  end if

100 format (2(F16.8,','))
end subroutine sld_elmpre

