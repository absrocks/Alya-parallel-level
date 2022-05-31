!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_ker_ILSA_sgsmodel.f90
!> @author  Oriol Lehmkuhl
!> @date    04/01/2018
!> @brief   ILSA sgs model
!> @details ILSA sgs model
!> see Dynamic subfilter-scale stress model for large-eddy simulations,
!>     A. Rouhi, U. Piomelli and B.J. Geurts, PHYSICAL REVIEW FLUIDS 1, 044401 (2016) 
!> @} 
!-----------------------------------------------------------------------
module mod_ker_ILSA_sgsmodel
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo

  implicit none

  real(rp),   pointer, save :: ck_square(:,:)   ! ck^2
  real(rp),   pointer, save :: ave_kres(:,:)    ! time average resolved kinetic energy 
  real(rp),   pointer, save :: ave_etot(:,:)    ! time average total dissipation 

  real(rp),   pointer, save :: ave_x1(:,:)      ! time average x1 eq(11) 
  real(rp),   pointer, save :: ave_x2(:,:)      ! time average x2 eq(11) 
  real(rp),   pointer, save :: ave_x3(:,:)      ! time average x3 eq(11) 
  real(rp),   pointer, save :: ave_x4(:,:)      ! time average x4 eq(11) 

  real(rp),   pointer, save :: ave_sij(:,:,:,:) ! time average  sij
  real(rp),   pointer, save :: ave_vel(:,:,:)   ! time average  vel

  real(rp),   pointer, save :: gpvsfs(:,:)      ! sfs viscosity
  real(rp),   pointer, save :: gplest(:,:)      ! Lest
  real(rp),   pointer, save :: gpstrain(:,:)    ! |strain|
  real(rp),   pointer, save :: eliti(:)         ! local averaging time

  real(rp),   pointer, save :: delta_f(:,:)     ! Estimated delta of the filter
  real(rp),   pointer, save :: ratio_l_df(:,:)  ! ratio between L/Delta_f
  real(rp),   pointer, save :: ratio_df_h(:,:)   ! ratio between Delta_f/h
  real(rp),   pointer, save :: ave_grad(:,:)    ! average ratio between (Grad·Grad)
  real(rp),   pointer, save :: ave_rij(:,:)     ! average ratio between (Rij·Rij)
  real(rp),   pointer, save :: d_stau(:,:)       ! dynamic_stau

  private

  public :: ker_ILSA_sgs_viscosity
  public :: ker_ILSA_temporal_average
  public :: ker_ILSA_solve_ck
  public :: ker_ILSA_get_ratio_ldf
  public :: ker_ILSA_get_ratio_dfh
  public :: ker_ILSA_get_lest
  public :: ker_ILSA_get_deltaf
  public :: ker_ILSA_get_ck

  real(rp), save :: s_tau
  real(rp), save :: dfilter
  real(rp), save :: T_ave ! a de venir de fora
  real(rp), save :: gzero

contains 


  subroutine ker_ILSA_allocate()

    integer(ip), save :: ipass = 0
    integer(ip)   :: igaus,idime,jdime,ielem

    !s_tau = 0.022_rp
    gzero = 1.0e-6_rp ! zeror
    !gzero = 1.0e-10_rp ! zeror
  

    if( ipass == 0 ) then

       ipass = 1

       nullify(ck_square)
       nullify(ave_kres)
       nullify(ave_etot)

       nullify(ave_x1)
       nullify(ave_x2)
       nullify(ave_x3)
       nullify(ave_x4)

       nullify(ave_sij)
       nullify(ave_vel)

       nullify(gpvsfs)
       nullify(gplest)
       nullify(gpstrain)

       nullify(delta_f)
       nullify(ratio_l_df)
       nullify(ratio_df_h)
       nullify(ave_grad)
       nullify(ave_rij)
       nullify(d_stau)

       call memory_alloca(mem_modul(1:2,modul),'CKSQUARE','nsi_memall',ck_square,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEKRES','nsi_memall',ave_kres,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEETOT','nsi_memall',ave_etot,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEX1','nsi_memall',ave_x1,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEX2','nsi_memall',ave_x2,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEX3','nsi_memall',ave_x3,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEX4','nsi_memall',ave_x4,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVGRAD','nsi_memall',ave_grad,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVRIJ','nsi_memall',ave_rij,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'GPVSFS','nsi_memall',gpvsfs,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'GPLEST','nsi_memall',gplest,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'GPSTRAIN','nsi_memall',gpstrain,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'AVEVEL','nsi_memall',ave_vel,nelem,mgaus,ndime)
       call memory_alloca(mem_modul(1:2,modul),'AVESIJ','nsi_memall',ave_sij,nelem,mgaus,ndime,ndime)
       call memory_alloca(mem_modul(1:2,modul),'ELITI','nsi_memall',eliti,nelem)
       call memory_alloca(mem_modul(1:2,modul),'DELTAF','nsi_memall',delta_f,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'RATIOLDF','nsi_memall',ratio_l_df,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'RATIODFH','nsi_memall',ratio_df_h,nelem,mgaus)
       call memory_alloca(mem_modul(1:2,modul),'DSTAU','nsi_memall',d_stau,nelem,mgaus)
       
       do ielem = 1,nelem
          eliti    (ielem) = 0.0_rp
          do igaus = 1,mgaus
             ck_square (ielem,igaus) = 0.0_rp
             ave_kres  (ielem,igaus) = 0.0_rp
             ave_etot  (ielem,igaus) = 0.0_rp  
             ave_x1    (ielem,igaus) = 0.0_rp  
             ave_x2    (ielem,igaus) = 0.0_rp 
             ave_x3    (ielem,igaus) = 0.0_rp   
             ave_x4    (ielem,igaus) = 0.0_rp    
             gpvsfs    (ielem,igaus) = 0.0_rp
             gplest    (ielem,igaus) = 0.0_rp
             gpstrain  (ielem,igaus) = 0.0_rp
             delta_f   (ielem,igaus) = 0.0_rp
             ratio_l_df(ielem,igaus) = 0.0_rp
             ratio_df_h(ielem,igaus) = 0.0_rp
             ave_grad  (ielem,igaus) = 0.0_rp    
             ave_rij   (ielem,igaus) = 0.0_rp    
             d_stau    (ielem,igaus) = 0.0_rp
             do idime = 1,ndime
                ave_vel(ielem,igaus,idime) = 0.0_rp
                do jdime = 1,ndime
                   ave_sij(ielem,igaus,idime,jdime) = 0.0_rp
                end do
             end do
          end do
       end do

    end if

 end subroutine ker_ILSA_allocate

  subroutine ker_ILSA_sgs_viscosity(pgaus,igaui,igauf,gpvis,gpmut,gpvel,gpgve,ielem,const1,const2,hleng)

     integer(ip),  intent(in)    :: pgaus,igaui,igauf
     integer(ip),  intent(in)    :: ielem
     real(rp),     intent(in)    :: gpvel(ndime,pgaus)
     real(rp),     intent(in)    :: gpgve(ndime,ndime,pgaus)
     real(rp),     intent(in)    :: gpvis(pgaus)
     real(rp),     intent(inout) :: gpmut(pgaus)
     real(rp),     intent(in)  :: const1,const2
     real(rp),     intent(in)  :: hleng(3)
     integer(ip)   :: igaus

     T_ave = const1 
     s_tau = const2 
     !dfilter = const2

     call ker_ILSA_allocate()
     call ker_ILSA_temporal_average(pgaus,igaui,igauf,gpvis,gpmut,gpvel,gpgve,ielem,hleng)
     call ker_ILSA_solve_ck(pgaus,igaui,igauf,gpvis,gpmut,gpvel,gpgve,ielem)

     ! loop para calcular la mu_sgs

     do igaus = igaui,igauf
        !gpmut(igaus)        = max(ck_square(ielem,igaus)*gpstrain(ielem,igaus),0.0_rp)
        gpmut(igaus)        = max(ck_square(ielem,igaus)*((ave_kres(ielem,igaus)**3_ip)/(ave_etot(ielem,igaus)**2_ip+gzero))*gpstrain(ielem,igaus),0.0_rp)
        !gpmut(igaus)        = min(gpvis(igaus)*1e6_rp, gpmut(igaus))
        gpvsfs(ielem,igaus) = gpmut(igaus)
     end do

  end subroutine ker_ILSA_sgs_viscosity   

  subroutine ker_ILSA_temporal_average(pgaus,igaui,igauf,gpvis,gpmut,gpvel,gpgve,ielem,hleng)

     integer(ip),  intent(in)  :: pgaus,igaui,igauf
     integer(ip),  intent(in)  :: ielem
     real(rp),     intent(in) :: gpvel(ndime,pgaus)
     real(rp),     intent(in) :: gpgve(ndime,ndime,pgaus)
     real(rp),     intent(in) :: gpvis(pgaus)
     real(rp),     intent(inout) :: gpmut(pgaus)
     real(rp),     intent(in)  :: hleng(3)
     integer(ip)   :: igaus,inode,idime,jdime,kdime
     real(rp)      ::  gpkres(pgaus),gpvelf(ndime,pgaus)
     real(rp)      :: gpsij(ndime,ndime,pgaus),gpepst(pgaus),gprij(ndime,ndime,pgaus)
     real(rp)      :: ave_epsilon, dt, gpsijf(ndime,ndime,pgaus),gpx1(pgaus),gpx2(pgaus),gpx3(pgaus)
     real(rp)      :: gpgrad(ndime,ndime,pgaus)
     real(rp)                  :: xmile,seci4,g2_ij(ndime,ndime), &
       g2_kk,seci5,sd_ij,denom,G__ij(3,3), aux_strain, &
       G_val(3),vdumm(3,3),sigma(3),ynorm,fvaDr,Smagc,coeff,xnutu
     real(rp)      :: op,alpha,Bbeta
     integer(ip)   :: idumi

     dt = 1.0_rp/(dtinv+gzero)
     ave_epsilon = dt/(T_ave+dt+gzero)
     if(eliti(ielem) > T_ave) then
        eliti(ielem) = 0.0_rp
     else
        eliti(ielem) = eliti(ielem) + dt
     end if

     do igaus = igaui,igauf
     ! get u and sij at gauss points
     ! gpgvel evaluated on turbu no se si duplicar
        do idime = 1,ndime
           do jdime = 1,ndime
              gpsij(idime,jdime,igaus) = 0.5_rp*(gpgve(jdime,idime,igaus)+gpgve(idime,jdime,igaus))
           end do
        end do

        gpstrain(ielem,igaus) = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
            gpstrain(ielem,igaus) = gpstrain(ielem,igaus) + (gpsij(idime,jdime,igaus)*gpsij(idime,jdime,igaus))
           end do
        end do
        gpstrain(ielem,igaus) = sqrt(2.0_rp*gpstrain(ielem,igaus))

     ! u mean and Sij mean
        do idime = 1,ndime
           ave_vel(ielem,igaus,idime) = (dt*gpvel(idime,igaus) + eliti(ielem)*ave_vel(ielem,igaus,idime))/(dt+eliti(ielem))
           do jdime = 1,ndime
              ave_sij(ielem,igaus,idime,jdime) = (dt*gpsij(idime,jdime,igaus) + eliti(ielem)*ave_sij(ielem,igaus,idime,jdime))/(dt+eliti(ielem))
           end do
        end do

     ! k res, epsilon total and R^a_{ij}

        ! fluctuating fields
        gpkres(igaus) = 0.0_rp
        do idime = 1,ndime
           gpvelf(idime,igaus) = gpvel(idime,igaus)-ave_vel(ielem,igaus,idime)
           gpkres(igaus)       = gpkres(igaus) + 0.5_rp*gpvelf(idime,igaus)*gpvelf(idime,igaus)
           do jdime = 1,ndime
               gpsijf(idime,jdime,igaus) = gpsij(idime,jdime,igaus)-ave_sij(ielem,igaus,idime,jdime) 
           end do
        end do

        ! average k and epsilon
        ave_kres(ielem,igaus) = (dt*gpkres(igaus) + eliti(ielem)*ave_kres(ielem,igaus))/(dt+eliti(ielem))

        gpepst(igaus) = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              gpepst(igaus) = gpepst(igaus) + 2.0_rp*(gpvis(igaus)+gpvsfs(ielem,igaus))*gpsijf(idime,jdime,igaus)*gpsijf(idime,jdime,igaus) 
           end do
        end do
        ave_etot(ielem,igaus) = (dt*gpepst(igaus) + eliti(ielem)*ave_etot(ielem,igaus))/(dt+eliti(ielem))

        ! estimated integral length scale
        !gplest(ielem,igaus) = ((ave_kres(ielem,igaus)**3_ip)/(ave_etot(ielem,igaus)**2_ip+gzero))
        gplest(ielem,igaus) = sqrt(ave_kres(ielem,igaus)**3_ip)/(gzero+ave_etot(ielem,igaus))

        ! R^a_{ij}
        gprij(1:ndime,1:ndime,igaus) =  0.0_rp

        do idime = 1,ndime
           gprij(idime,idime,igaus) =  - (1.0_rp/3.0_rp)*gpvelf(idime,igaus)*gpvelf(idime,igaus)
           do jdime = 1,ndime
              gprij(idime,jdime,igaus) = gprij(idime,jdime,igaus) + gpvelf(idime,igaus)*gpvelf(jdime,igaus)
           end do
        end do

        ! X1,X2 and X3 

        ! test wale
        ! g2_ij = 0.0_rp
        ! g2_kk = 0.0_rp

        ! do idime = 1,ndime
        !    do jdime = 1,ndime
        !       do kdime = 1,ndime
        !          g2_ij(idime,jdime) = g2_ij(idime,jdime) + gpgve(idime,kdime,igaus)*gpgve(kdime,jdime,igaus)
        !       end do
        !    end do
        ! end do

        ! g2_kk = g2_ij(1,1) + g2_ij(2,2) + g2_ij(3,3)

        ! seci4 = 0.0_rp
        ! seci5 = 0.0_rp
        ! sd_ij = 0.0_rp
        ! denom = 0.0_rp

        ! do idime = 1,ndime
        !    do jdime = 1,ndime
        !       seci4 = seci4 + 0.5_rp*gpgve(idime,jdime,igaus) * &          ! S_ij : S_ij
        !            (gpgve(idime,jdime,igaus) + gpgve(jdime,idime,igaus)) 
        !       sd_ij = 0.5_rp * ( g2_ij(idime,jdime) + g2_ij(jdime,idime) ) 
        !       seci5 = seci5 + sd_ij*sd_ij     ! Sd_ij:Sd_ij
        !    end do

        ! end do
        ! seci5 = seci5 - g2_kk*g2_kk/3.0_rp

        ! denom = max ( zeror , seci4**2.5_rp + seci5**1.25_rp )

        ! if(denom > gzero) then
        !    op = (seci5**1.5_rp)/denom
        ! else
        !    op = 0.0_rp
        ! end if

          G__ij = 0.0_rp ! G = g^T*g
          alpha = 0.0
          do idime = 1_ip,3_ip
             do jdime = 1_ip,3_ip
                G__ij(idime,jdime) = G__ij(idime,jdime) + gpgve(idime,jdime,igaus)*gpgve(jdime,idime,igaus)
                alpha = alpha + gpgve(jdime,idime,igaus)*gpgve(idime,jdime,igaus)
             end do
          end do

          !alpha = (G__ij(1_ip,1_ip) + G__ij(2_ip,2_ip) + G__ij(3_ip,3_ip))
          Bbeta =  G__ij(1_ip,1_ip)*G__ij(2_ip,2_ip) + G__ij(2_ip,2_ip)*G__ij(3_ip,3_ip) + G__ij(3_ip,3_ip)*G__ij(1_ip,1_ip) &
               - G__ij(1_ip,2_ip)*G__ij(1_ip,2_ip) - G__ij(2_ip,3_ip)*G__ij(2_ip,3_ip) - G__ij(1_ip,3_ip)*G__ij(1_ip,3_ip)

          if ( alpha > gzero )then                         ! Avoid divide by zero
             op = sqrt (max(( Bbeta ) / ( alpha ), 0.0_rp))
          else
             op = 0.0_rp
          end if

        ! end test wale

        gpx1(igaus) = 2.0_rp*(gplest(ielem,igaus)**4_ip)*(gpstrain(ielem,igaus)**2_ip)*(op**2_ip)
        gpx2(igaus) = 0.0_rp
        gpx3(igaus) = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              gpx2(igaus) = gpx2(igaus) + 4.0_rp*(gplest(ielem,igaus)**2_ip)*op*(gpsij(idime,jdime,igaus)*gprij(idime,jdime,igaus))
              gpx3(igaus) = gpx3(igaus) + (gprij(idime,jdime,igaus)*gprij(idime,jdime,igaus))
           end do
        end do

        ave_x1(ielem,igaus) = (dt*gpx1(igaus) + eliti(ielem)*ave_x1(ielem,igaus))/(dt+eliti(ielem))
        ave_x2(ielem,igaus) = (dt*gpx2(igaus) + eliti(ielem)*ave_x2(ielem,igaus))/(dt+eliti(ielem))
        ave_x3(ielem,igaus) = (dt*gpx3(igaus) + eliti(ielem)*ave_x3(ielem,igaus))/(dt+eliti(ielem))

        gpstrain(ielem,igaus) = op

        ! Estimated Delta_f and ratio L/Delta

        !gpgrad(1:ndime,1:ndime,igaus) =  0.0_rp
        !gprij(1:ndime,1:ndime,igaus) =  0.0_rp
        !do idime = 1,ndime
        !   do jdime = 1,ndime
        !         gpgrad(idime,jdime,igaus) = gpgrad(idime,jdime,igaus) + gpgve(jdime,idime,igaus)*gpgve(idime,jdime,igaus)
        !      gprij(idime,jdime,igaus) = gprij(idime,jdime,igaus) + gpvelf(idime,igaus)*gpvelf(jdime,igaus) - 2.0_rp*gpvsfs(ielem,igaus)*gpsij(idime,jdime,igaus)
        !   end do
        !end do

        !gpx2(igaus) = 0.0_rp
        !gpx3(igaus) = 0.0_rp
        !do idime = 1,ndime
        !   do jdime = 1,ndime
        !      gpx2(igaus) = gpx2(igaus) + (gpgrad(idime,jdime,igaus)*gpgrad(idime,jdime,igaus))
        !      gpx3(igaus) = gpx3(igaus) + (gprij(idime,jdime,igaus)*gprij(idime,jdime,igaus))
        !   end do
        !end do
        !ave_grad(ielem,igaus) = (dt*(gpx2(igaus)) + eliti(ielem)*ave_grad(ielem,igaus))/(dt+eliti(ielem))
        !ave_rij(ielem,igaus)  = (dt*(gpx3(igaus)) + eliti(ielem)*ave_rij(ielem,igaus))/(dt+eliti(ielem))

        !if(ave_rij(ielem,igaus)>gzero) then
        !    d_stau(ielem,igaus) = min(2.0_rp*dfilter**2*sqrt(ave_grad(ielem,igaus)/ave_rij(ielem,igaus)), 0.16)
        !else
        !    d_stau(ielem,igaus) = 0.02_rp
        !end if
        d_stau(ielem,igaus) = s_tau

     end do
  end subroutine ker_ILSA_temporal_average  

  subroutine ker_ILSA_solve_ck(pgaus,igaui,igauf,gpvis,gpmut,gpvel,gpgve,ielem)

     integer(ip),  intent(in)  :: pgaus,igaui,igauf
     integer(ip),  intent(in)  :: ielem
     real(rp),     intent(in) :: gpvel(ndime,pgaus)
     real(rp),     intent(in) :: gpvis(pgaus)
     real(rp),     intent(in) :: gpgve(ndime,ndime,pgaus)
     real(rp),     intent(inout) :: gpmut(pgaus)
     real(rp) :: a,b,c,d,x1,x2
     integer(ip)   :: igaus,inode

     a = 0.0_rp
     b = 0.0_rp
     c = 0.0_rp

     ! ax^2 + bx + c = 0

     do igaus = igaui,igauf

        ! d = ave_x2(ielem,igaus)**2_ip + 4.0_rp*((1.0/s_tau)**2_ip - 1.0_rp)*ave_x1(ielem,igaus)*ave_x3(ielem,igaus) 

        ! if(abs(ave_x1(ielem,igaus)) > gzero) then
        !    ck_square(ielem,igaus) = (sqrt(d)-ave_x2(ielem,igaus))/(2.0_rp*ave_x1(ielem,igaus)*((1.0_rp/s_tau)**2_ip - 1.0_rp))
        ! else
        !    if(abs(ave_x2(ielem,igaus))>gzero) then
        !       ck_square(ielem,igaus) = ave_x3(ielem,igaus)/ave_x2(ielem,igaus) 
        !    else
        !       ck_square(ielem,igaus) = 0.0_rp
        !    end if
        ! end if
        ! ck_square(ielem,igaus) = max(ck_square(ielem,igaus), 0.0_rp)

        a = ave_x1(ielem,igaus)*(1.0_rp - (1.0_rp/d_stau(ielem,igaus))**2_ip)
        b = -ave_x2(ielem,igaus)
        c = ave_x3(ielem,igaus)

        !if (abs(a)<gzero) then
        !   if(abs(c)<gzero) then
        !      x1 = 0.0_rp
        !      x2 = 0.0_rp
        !   else
        !      x1 = 2.0_rp*c/(-b + sqrt(b*b - 4.0_rp*a*c))
        !      x2 = 2.0_rp*c/(-b - sqrt(b*b - 4.0_rp*a*c))
        !   end if
        !else
        !   x1 = (-b + sqrt((b*b)-(4.0_rp*a*c))) / ((2.0_rp*a))
        !   x2 = (-b - sqrt((b*b)-(4.0_rp*a*c))) / ((2.0_rp*a))
        !end if
        !d = (b*b)/(b*b - 4.0_rp*a*c)
        !if(abs(d)<10.0_rp) then
        !   if(abs(x1)>abs(x2)) then
        !      x2 = (c/a)/x1
        !   else 
        !      x1 = (c/a)/x2
        !   end if
        !end if

        if(b> 0.0_rp) then
           d = -0.5*(b + sqrt(b*b - 4.0_rp*a*c))
        else 
           d = -0.5*(b - sqrt(b*b - 4.0_rp*a*c))
        end if

        if(abs(d)> gzero) then
           if(abs(a)>gzero) then 
              x1 = d / a
              x2 = c / d   
           else
              x1 = 0.0_rp
              x2 = c / d
           end if
        else
           x1 = 0.0_rp
           x2 = 0.0_rp 
        end if

        ck_square(ielem,igaus) = max(max(x1,x2), 0.0_rp)
     end do
  end subroutine ker_ILSA_solve_ck 

  subroutine ker_ILSA_get_ratio_ldf(pgaus,igaui,igauf,ielem,gpratio)

    integer(ip),  intent(in)  :: pgaus,igaui,igauf
    integer(ip),  intent(in)  :: ielem
    real(rp),     intent(inout) :: gpratio(pgaus)
    integer(ip)   :: igaus

    do igaus = igaui,igauf
       gpratio(igaus) = ave_etot(ielem,igaus)
    end do
 end subroutine ker_ILSA_get_ratio_ldf 

 subroutine ker_ILSA_get_ratio_dfh(pgaus,igaui,igauf,ielem,gpratio)

    integer(ip),  intent(in)  :: pgaus,igaui,igauf
    integer(ip),  intent(in)  :: ielem
    real(rp),     intent(inout) :: gpratio(pgaus)
    integer(ip)   :: igaus

   do igaus = igaui,igauf
       gpratio(igaus) = ave_etot(ielem,igaus)
    end do
 end subroutine ker_ILSA_get_ratio_dfh 

  subroutine ker_ILSA_get_lest(pgaus,igaui,igauf,ielem,gpratio)

    integer(ip),  intent(in)  :: pgaus,igaui,igauf
    integer(ip),  intent(in)  :: ielem
    real(rp),     intent(inout) :: gpratio(pgaus)
    integer(ip)   :: igaus

    do igaus = igaui,igauf
       gpratio(igaus) = gplest(ielem,igaus)
    end do
 end subroutine ker_ILSA_get_lest 

 subroutine ker_ILSA_get_deltaf(pgaus,igaui,igauf,ielem,gpratio)

    integer(ip),  intent(in)  :: pgaus,igaui,igauf
    integer(ip),  intent(in)  :: ielem
    real(rp),     intent(inout) :: gpratio(pgaus)
    integer(ip)   :: igaus

   do igaus = igaui,igauf
       gpratio(igaus) = ave_kres(ielem,igaus)
    end do
 end subroutine ker_ILSA_get_deltaf 

 subroutine ker_ILSA_get_ck(pgaus,igaui,igauf,ielem,gpratio)

    integer(ip),  intent(in)  :: pgaus,igaui,igauf
    integer(ip),  intent(in)  :: ielem
    real(rp),     intent(inout) :: gpratio(pgaus)
    integer(ip)   :: igaus

   do igaus = igaui,igauf
       gpratio(igaus) = sqrt(ck_square(ielem,igaus))
    end do
 end subroutine ker_ILSA_get_ck 

end module mod_ker_ILSA_sgsmodel
