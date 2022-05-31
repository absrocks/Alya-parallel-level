!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_ker_nsw_viscf90
!> @author  Herbert Owen
!> @brief   Obtains viscosity for no slip wall law 
!> @details - 
!>          - 
!> @} 
!-----------------------------------------------------------------------
module mod_ker_nsw_visc
  use def_kintyp,                   only : ip,rp
#ifndef VECTOR_SIZE
  use def_master,                   only :  VECTOR_SIZE
#endif
  Use def_master,                   only : ittim,kfl_paral    ! solucion trucha
  use def_kermod,                   only : kfl_nswel_ker,avwei_ker,avta1_nsw_ker,fact_nsw_ker,lnsw_exch,kfl_rough,rough_dom
  use def_domain,                   only : lnods,ltype,nnode,rough,lpoty,mnode

  implicit none
  private
  public :: ker_nsw_visc

  contains   

subroutine ker_nsw_visc(ndime,pnode,pgaus,list_elements,elavv,gpcar,elnnsw,gpvis,gpden,gpmut,elibopo,elywal,elcod,gpvis_nsw)   ! elcod just for debugging REMOVE

  integer(ip), intent(in)    :: ndime,pnode,pgaus
  integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
  real(rp),    intent(in)    :: elavv(VECTOR_SIZE,ndime,pnode)
  real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
  real(rp),    intent(in)    :: elnnsw(VECTOR_SIZE,ndime)
  real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
  real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
  real(rp),    intent(in)    :: gpmut(VECTOR_SIZE,pgaus)
  real(rp),    intent(in)    :: elibopo(VECTOR_SIZE,pnode)
  real(rp),    intent(in)    :: elywal(VECTOR_SIZE)
  real(rp),    intent(in)    :: elcod(VECTOR_SIZE,ndime,pnode)
  real(rp),    intent(out)   :: gpvis_nsw(VECTOR_SIZE,pgaus)


  integer(ip)      :: kount(VECTOR_SIZE)
  real(rp)         :: gpgnavv(VECTOR_SIZE,ndime,pgaus)  ! Gauss Point Gradient in Normal dir of AVerage Velocity 
  real(rp)         :: elgnavv(VECTOR_SIZE,ndime)        ! ELement Gradient in Normal dir of AVerage Velocity 
  real(rp)         :: elgnavvt(VECTOR_SIZE)             ! ELement Gradient in Normal dir of AVerage Tangent Velocity 

  real(rp)         :: auxvi(VECTOR_SIZE)
  real(rp)         :: tveno_aux(VECTOR_SIZE)
  real(rp)         :: auxde(VECTOR_SIZE)
  real(rp)         :: auxmut(VECTOR_SIZE)
  real(rp)         :: avelavv(VECTOR_SIZE,ndime)        ! AVerage ELement AVerage Velocity
  real(rp)         :: avta1_aux(VECTOR_SIZE,ndime)
  real(rp)         :: avtan_aux(VECTOR_SIZE,ndime)
  real(rp)         :: auxi(VECTOR_SIZE)
  real(rp)         :: velfr(VECTOR_SIZE)
  real(rp)         :: fact_aux(VECTOR_SIZE)

  real(rp)         :: avtan_fric_grad_based(VECTOR_SIZE)    ! just the magnitude
  real(rp)         :: av_mu_mut(VECTOR_SIZE)


  
  real(rp)         :: rough_aux,vikin,tveno,auxi2,kount1
  integer(ip)      :: inode,ivect,idime,jdime,igaus,ielem,ipoin,kk,ibopo
  integer(ip),parameter      :: imethod = 4   !now the default will be method 4 that is teh one taht is working best


#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

#ifdef OPENACC
    do ivect = 1,VECTOR_SIZE
#endif 
       !----------------------------------------------------------------------
       !
       ! Gauss point values
       !
       !----------------------------------------------------------------------
       !
       ! GPGNAVV = dj uav_i nj =  dj N_I_i nj  Uav_I  ! gauss point gradient in normal direction of the average velocity 
       !
       gpgnavv(DEF_VECT,:,:)   = 0.0_rp

       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgnavv(DEF_VECT,idime,igaus) = gpgnavv(DEF_VECT,idime,igaus) + elavv(DEF_VECT,idime,inode) * &
                        gpcar(DEF_VECT,jdime,inode,igaus) * elnnsw(DEF_VECT,jdime)
                end do
             end do
          end do
       end do
       !
       ! Obtain the average over gauss points of the normal gradient of the average velocity for all gauss points
       ! I could merge this loop with the upper one and obtain elgnavv directly
       ! Also average density and viscosity
       !
       elgnavv  = 0.0_rp
       auxvi = 0.0_rp
       auxde = 0.0_rp
       auxmut = 0.0_rp
       do igaus = 1,pgaus
          do idime = 1,ndime
!             elgnavvn(DEF_VECT) = elgnavvn(DEF_VECT) + gpgnavv(DEF_VECT,idime,igaus) * elnnsw(DEF_VECT,idime)
             elgnavv(DEF_VECT,idime)  = elgnavv(DEF_VECT,idime)  + gpgnavv(DEF_VECT,idime,igaus)
          end do
          auxvi(DEF_VECT)  = auxvi(DEF_VECT)  + gpvis(DEF_VECT,igaus)
          auxde(DEF_VECT)  = auxde(DEF_VECT)  + gpden(DEF_VECT,igaus)
          auxmut(DEF_VECT) = auxmut(DEF_VECT) + gpmut(DEF_VECT,igaus)
       end do
       do idime = 1,ndime
          elgnavv(DEF_VECT,idime) = elgnavv(DEF_VECT,idime) / real(pgaus,rp)
       end do
       auxvi(DEF_VECT)  = auxvi(DEF_VECT)  / real(pgaus,rp)
       auxde(DEF_VECT)  = auxde(DEF_VECT)  / real(pgaus,rp)
       auxmut(DEF_VECT) = auxmut(DEF_VECT) / real(pgaus,rp)
       !
       ! Substract normal component to keep only tangential one
       !
       !
       auxi(DEF_VECT) = 0.0_rp
       do idime=1,ndime
          auxi(DEF_VECT) = auxi(DEF_VECT) + elgnavv(DEF_VECT,idime) * elnnsw(DEF_VECT,idime)
       end do
       do idime=1,ndime
          elgnavv(DEF_VECT,idime) = elgnavv(DEF_VECT,idime) - auxi(DEF_VECT) * elnnsw(DEF_VECT,idime)
       end do
       !
       ! obtain modulus of tangential component
       !
       auxi(DEF_VECT) = 0.0_rp
       do idime=1,ndime
          auxi(DEF_VECT) = auxi(DEF_VECT) + elgnavv(DEF_VECT,idime) * elgnavv(DEF_VECT,idime)
       end do
       elgnavvt(DEF_VECT) = sqrt (auxi(DEF_VECT))




       if (.false.) then   ! this part will no longer be used
          !
          ! Obtain the average velocity for non boundary nodes
          !
          avelavv(DEF_VECT,:) = 0.0_rp
          kount = 0_ip
          do inode = 1,pnode
             !
             ! I do not need an if ibopo > 0 that would not be suitable in the vector case
             ! Instead I usar elibopo elibopo(DEF_VECT,inode)  that is =1.0 if ibopo>0 and 0.0 else
             !
             do idime=1,ndime
                avelavv(DEF_VECT,idime) = avelavv(DEF_VECT,idime) + elavv(DEF_VECT,idime,inode) * ( 1.0_rp - elibopo(DEF_VECT,inode) ) 
             end do
             kount(DEF_VECT) = kount(DEF_VECT) + nint(elibopo(DEF_VECT,inode))
          end do

          do idime=1,ndime
             avelavv(DEF_VECT,idime) = avelavv(DEF_VECT,idime) / (real(kount,rp)+1.0e-30_rp)   ! OJO trucheada para safar de que estoy haciendo cuentas en elementso interiores al pedo
          end do                                                                               ! esta relacionado con que no se como poner un if paar que haga algunos elem y otrso no
          !
          ! Substract normal component to keep only tangential one
          !
          !
          auxi(DEF_VECT) = 0.0_rp
          do idime=1,ndime
             auxi(DEF_VECT) = auxi(DEF_VECT) + avelavv(DEF_VECT,idime) * elnnsw(DEF_VECT,idime)
          end do
          do idime=1,ndime
             avelavv(DEF_VECT,idime) = avelavv(DEF_VECT,idime) - auxi(DEF_VECT) * elnnsw(DEF_VECT,idime)
          end do
       else  ! NEW OPTION
          !This comes directly in lnsw_exch(ielem)%velav(1:ndime)  - obtained in nsi_wallav -- ojo se hace al final elpaso de tiempo en restrat habra que volver a calcularlo - ya lo corregi
          ! for the moment it is not vectorized - not sure how easy/practical it would be 
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if (ielem/=0_ip) then
                avelavv(ivect,1:ndime) = lnsw_exch(ielem)%velav(1:ndime)
             else
                avelavv(ivect,1:ndime) = 0.0_rp
             end if
          end do
          !
          ! Substract normal component to keep only tangential one
          !
          !
          auxi(DEF_VECT) = 0.0_rp
          do idime=1,ndime
             auxi(DEF_VECT) = auxi(DEF_VECT) + avelavv(DEF_VECT,idime) * elnnsw(DEF_VECT,idime)
          end do
          do idime=1,ndime
             avelavv(DEF_VECT,idime) = avelavv(DEF_VECT,idime) - auxi(DEF_VECT) * elnnsw(DEF_VECT,idime)
          end do
       end if

#ifdef OPENACC
    end do
#endif
    !
    ! Obtain tange from wall_law.
    ! Compute U*: VELFR   this part is not vectorized for the moment  I would need a vector frivel , not difficult
    ! also Time average of (mu+mut) d u_t / dn  - to be used later
    !
    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then
          if(kfl_nswel_ker(ielem) >0) then
             vikin = auxvi(ivect) / auxde(ivect)                           ! nu
             call vecnor(avelavv(ivect,1:ndime),ndime,tveno,2_ip)          ! |u_tan-u_fix_tan|
             if(tveno > 1.0e-20_rp) then
                tveno_aux(ivect) = tveno    ! I save it for later use
             else
                tveno_aux(ivect) = 1.0_rp
             end if

             if( kfl_rough == 0 ) then 
                !
                ! Constant roughness
                !
                rough_aux = rough_dom

             else if( kfl_rough > 0 ) then
                rough_aux = 0.0_rp
                kount1 = 0_ip
                do inode = 1,pnode
                   ipoin =  lnods(inode,ielem)
                   ibopo = lpoty(ipoin)
                   if (ibopo /= 0) then
                      kount1 = kount1 + 1_ip
                      rough_aux = rough_aux + rough(ipoin)
                   end if
                end do
                rough_aux = rough_aux / kount1
!             else
!                rough_aux = 1.0_rp   ! this value will not be used    - but debug might complain if no value is given
             end if
             
             call frivel(elywal(ivect),rough_aux,tveno,vikin,velfr(ivect))      ! U*
             !
             ! Time average of (mu+mut) d u_t / dn  - to be used later
             !
             avta1_aux(ivect,:) = avta1_nsw_ker(:,kfl_nswel_ker(ielem))
             !
             ! fact_aux to be used later
             !
             kk = 0_ip
             auxi2 = 0.0_rp
             do inode=1,pnode
                ipoin = lnods(inode,ielem)
                if (fact_nsw_ker(ipoin) > 0.0_rp) then    ! boundary node
                   kk = kk + 1_ip 
                   auxi2 = auxi2 + fact_nsw_ker(ipoin)
                end if
                if (kk > 0_ip ) then
                   fact_aux(ivect) = auxi2 / real(kk,rp)
                else
                   fact_aux(ivect) = 1.0_rp
                end if
             end do   

          else   ! Element is not boundary element
             velfr(ivect) = 0.0_rp
             tveno_aux(ivect) = 1.0_rp   ! just some value so that it does not give problems
             avta1_aux(ivect,:) = 0.0_rp 
          end if
       else   ! Element number is null
          velfr(ivect) = 0.0_rp
          tveno_aux(ivect) = 1.0_rp      ! just some value so that it does not give problems
          avta1_aux(ivect,:) = 0.0_rp
       end if
    end do

#ifdef OPENACC
    do ivect = 1,VECTOR_SIZE
#endif
       if (imethod==1) then    ! this is the way I had it originally
          !
          ! mu_nsw = (du/dy)**(-1) ( (fact*Tau ) - (mu+mut) (du/dy) )   
          !
          !
          ! fact * Tau - (mu+mut) * du_t/dn   - all average
          !   
          auxi(DEF_VECT) = fact_aux(DEF_VECT) * velfr(DEF_VECT) * velfr(DEF_VECT) * auxde(DEF_VECT) / tveno_aux(DEF_VECT)
          !
          ! Estos nombres enstan muy chotos -- mejorarlos
          ! avtan_aux -> avtan_muaddit
          ! avta1_aux -> avtan_mu_mut
          ! -auxi(DEF_VECT) * avelavv(DEF_VECT,idime)  -- crear un vector avtan_frivel
          !
          do idime=1,ndime
             avtan_aux(DEF_VECT,idime) = - auxi(DEF_VECT) * avelavv(DEF_VECT,idime) - avta1_aux(DEF_VECT,idime)    ! tau - (mu+mut) * du_t/dn   - all average
          end do
          !
          ! magnitude ( fact * tau - (mu+mut) * du_t/dn )
          !
          auxi(DEF_VECT) = 0.0_rp
          do idime=1,ndime
             auxi(DEF_VECT) = auxi(DEF_VECT) + avtan_aux(DEF_VECT,idime) * avtan_aux(DEF_VECT,idime)
          end do
          auxi(DEF_VECT) = sqrt(auxi(DEF_VECT))
          !
          ! divide by (du/dy)
          !
          do igaus=1,pgaus   ! constant all the same value
             gpvis_nsw(DEF_VECT,igaus) =  auxi(DEF_VECT) / abs(elgnavvt(DEF_VECT)+1.0e-30_rp)   ! to avoid divide by zero perhaps there is a better option
          end do

       else if (imethod==2) then

          !
          ! mu_nsw = (du/dy)**(-1) ( (fact*Tau ) - (mu+mut) (du/dy) )
          ! I simplify it to
          ! mu_nsw = (du/dy)**(-1) (fact*Tau) - (mu+mut)
          !
          ! Now thing get simpler I can simply get mod(tau)  and avoid working with vectors I directly obtain the magnitude of tauin auxi
          !rethink if in some case I may be commiting errors for not working vectorly
          !
          ! I have realized that the two options are not totally equivalent -- when I substarct  avta1_aux(DEF_VECT,idime)   on top this is teh average value of (mu+mut)du/dy
          ! While here I am just substarctig the instantaneous value of mu+mut   -- I need to rethink which is teh most correct option
          ! moreover whne I do a restart their values are significantly diferent I need to recheck everything.
          ! creo que lo correcto es lo anterior . fact = AVGTR/AVNTR   - TAU_FRIC es tambien el valor average
          ! con lo cual lo que obtengo en AV_TAU_FRIC = AVNTR = AVGTR/fact = (1/fact)  * (AVGTR_(mu+mut) + AVGTR_(muadic))
          ! de donde puedo sacar AVGTR_(muadic) = fact * AV_TAU_FRIC  - AVGTR_(mu+mut)
          ! de donde, usando  AVGTR_(muadic) = av_muaddic * AV_du/dy ,  se obtiene finalmente
          ! av_muaddic = ( 1/AV_du/dy ) * ( fact * AV_TAU_FRIC  - AVGTR_(mu+mut) )
          !

          
          !
          ! Tangential force that comes from the friction velocity - gradient based because it has already been transformed using fact
          !
          avtan_fric_grad_based(DEF_VECT) = fact_aux(DEF_VECT) * velfr(DEF_VECT) * velfr(DEF_VECT) * auxde(DEF_VECT)
 !         if (velfr(1)> 1.0e-12) write(kfl_paral+2000,*) ' avtan_fric_grad_based(1) , fact_aux(1) , velfr(1), auxde(1)', avtan_fric_grad_based(1) , fact_aux(1) , velfr(1), auxde(1)
          !
          ! divide by (du/dy) and substract mu+mut  actually to work in an average sense I can not substract mu+mut that are instantaneous  
          !
          ! Actaully it might be better to obtain the value for gp=1 and at the end set all gp to the same value
!          do igaus=1,pgaus   ! constant all the same value
           gpvis_nsw(DEF_VECT,1) = ( avtan_fric_grad_based(DEF_VECT) / abs(elgnavvt(DEF_VECT)+1.0e-30_rp) ) !!!- (auxmut(DEF_VECT) + auxvi(DEF_VECT))
!          end do
          !
          ! Here I am going to sustract average mu+mut --- av_mu_mut
          ! to obtain it I first obtain modulus of avta1_aux(DEF_VECT,idime)
          !
          av_mu_mut(DEF_VECT) = avta1_aux(DEF_VECT,1) * avta1_aux(DEF_VECT,1) +  avta1_aux(DEF_VECT,2) * avta1_aux(DEF_VECT,2)
          if (ndime==3) av_mu_mut(DEF_VECT) = av_mu_mut(DEF_VECT) +  avta1_aux(DEF_VECT,3) * avta1_aux(DEF_VECT,3)
          av_mu_mut(DEF_VECT) = sqrt(av_mu_mut(DEF_VECT)) / abs(elgnavvt(DEF_VECT)+1.0e-30_rp)
          !
          ! Now I substract -- moreover I take into account that the total traction can not be smaller than the one coming from mu+mut
          ! ESTE ES EL PASO QUE STAB FALTANDO (el max) - aca se supone que ambas tienne la misma direccion pero creo que es algo lógico
          ! pero se podría repensar
          !
          gpvis_nsw(DEF_VECT,1) = max( 0.0_rp , gpvis_nsw(DEF_VECT,1) - av_mu_mut(DEF_VECT) )
          !
          ! Put the same value in all gauss points
          !
          do igaus=2,pgaus   
             gpvis_nsw(DEF_VECT,igaus) = gpvis_nsw(DEF_VECT,1)
          end do

       else if (imethod==3  .or. imethod==4) then

          !
          ! Idem 2 but without  max( 0.0_rp ,    -- that is, I allow a negative aditional viscosity this can be dangerous but I think it is teh only way to recover the same as the normal wal law
          ! To make it less dangerous I will change max( 0.0_rp ,  to max( -gpmut ,
          !
          
          !
          ! Tangential force that comes from the friction velocity - gradient based because it has already been transformed using fact
          !
          if (imethod==3) then   ! here is the only diference between method 3 and 4 
             avtan_fric_grad_based(DEF_VECT) = fact_aux(DEF_VECT) * velfr(DEF_VECT) * velfr(DEF_VECT) * auxde(DEF_VECT)
          else if (imethod==4) then
             avtan_fric_grad_based(DEF_VECT) =                       velfr(DEF_VECT) * velfr(DEF_VECT) * auxde(DEF_VECT)
          end if

 !         if (velfr(1)> 1.0e-12) write(kfl_paral+2000,*) ' avtan_fric_grad_based(1) , fact_aux(1) , velfr(1), auxde(1)', avtan_fric_grad_based(1) , fact_aux(1) , velfr(1), auxde(1)
          !
          ! divide by (du/dy) and substract mu+mut  actually to work in an average sense I can not substract mu+mut that are instantaneous  
          !
          ! Actaully it might be better to obtain the value for gp=1 and at the end set all gp to the same value
!          do igaus=1,pgaus   ! constant all the same value
          gpvis_nsw(DEF_VECT,1) = ( avtan_fric_grad_based(DEF_VECT) / (abs(elgnavvt(DEF_VECT))+1.0e-30_rp) ) !!!- (auxmut(DEF_VECT) + auxvi(DEF_VECT))
           !          end do
!           write(kfl_paral+2000,*) 'gpvis_nsw(DEF_VECT,1),avtan_fric_grad_based(DEF_VECT),abs(elgnavvt(DEF_VECT))',&
!                gpvis_nsw(DEF_VECT,1),avtan_fric_grad_based(DEF_VECT),abs(elgnavvt(DEF_VECT))
           
          !
          ! Here I am going to sustract average mu+mut --- av_mu_mut
          ! to obtain it I first obtain modulus of avta1_aux(DEF_VECT,idime)
          !
          av_mu_mut(DEF_VECT) = avta1_aux(DEF_VECT,1) * avta1_aux(DEF_VECT,1) +  avta1_aux(DEF_VECT,2) * avta1_aux(DEF_VECT,2)
          if (ndime==3) av_mu_mut(DEF_VECT) = av_mu_mut(DEF_VECT) +  avta1_aux(DEF_VECT,3) * avta1_aux(DEF_VECT,3)
          av_mu_mut(DEF_VECT) = sqrt(av_mu_mut(DEF_VECT)) / (abs(elgnavvt(DEF_VECT))+1.0e-30_rp)
          !
          ! Now I substract -- moreover I take into account that the total traction can not be smaller than the one coming from mu+mut
          ! ESTE ES EL PASO QUE STAB FALTANDO (el max) - aca se supone que ambas tienne la misma direccion pero creo que es algo lógico
          ! pero se podría repensar
          !
          gpvis_nsw(DEF_VECT,1) = gpvis_nsw(DEF_VECT,1) - av_mu_mut(DEF_VECT)
!          write(kfl_paral+3000,*) 'gpvis_nsw(DEF_VECT,1),av_mu_mut(DEF_VECT),(gpmut(DEF_VECT,igaus) + gpvis(DEF_VECT,igaus))',&
!               gpvis_nsw(DEF_VECT,1),av_mu_mut(DEF_VECT),(gpmut(DEF_VECT,igaus) + gpvis(DEF_VECT,igaus))
          !
          ! max between  - gpmut(DEF_VECT,igaus) and the calculated values that for the moment is in igaus =1  
          !
          do igaus=2,pgaus   
             gpvis_nsw(DEF_VECT,igaus) = max ( - gpmut(DEF_VECT,igaus)  , gpvis_nsw(DEF_VECT,1) )
          end do
          gpvis_nsw(DEF_VECT,1) = max ( - gpmut(DEF_VECT,1)  , gpvis_nsw(DEF_VECT,1) )
!          if (ittim == 3_ip) then
!             write(kfl_paral+5000,*)  gpvis_nsw(DEF_VECT,1)
!          end if

       else if (imethod==5) then

          !
          ! Similar  to method 3 but  now teh idea is that in teh first element nsw_vsic takes care of everything
          ! thus I calculate it directly. then I substract mu+mu_les (because they will be added later)
          ! finally I just correct so that it is not smaller than the laminar viscosity
          ! This has teh added effect taht teh toal viscosity is homogeneous in the element
          !

          
          !
          ! Tangential force that comes from the friction velocity - gradient based because it has already been transformed using fact
          !
          avtan_fric_grad_based(DEF_VECT) = fact_aux(DEF_VECT) * velfr(DEF_VECT) * velfr(DEF_VECT) * auxde(DEF_VECT)
 !         if (velfr(1)> 1.0e-12) write(kfl_paral+2000,*) ' avtan_fric_grad_based(1) , fact_aux(1) , velfr(1), auxde(1)', avtan_fric_grad_based(1) , fact_aux(1) , velfr(1), auxde(1)
          !
          ! divide by (du/dy) and substract mu+mut  actually to work in an average sense I can not substract mu+mut that are instantaneous  
          !
          ! Actaully it might be better to obtain the value for gp=1 and at the end set all gp to the same value
          !          do igaus=1,pgaus   ! constant all the same value
          gpvis_nsw(DEF_VECT,1) = ( avtan_fric_grad_based(DEF_VECT) / abs(elgnavvt(DEF_VECT)+1.0e-30_rp) ) !!!- (auxmut(DEF_VECT) + auxvi(DEF_VECT))
          !          end do
          !
          ! Up to here it is the same as method 3 but here I will subtract mu+mut    while there I subtracted   average(mu+mut)
          ! so tthe total vsicosity close to teh wall will be much more cosntant in time in this case than in 3
          ! in 3 it had a more variable behaviour due to the addition of mut later
          ! what will be teh effect is not totally clear to me but I guess it will be nice to see

          !
          ! Now I substract -- moreover I take into account that the total traction can not be smaller than the one coming from mu+mut
          ! ESTE ES EL PASO QUE STAB FALTANDO (el max) - aca se supone que ambas tienne la misma direccion pero creo que es algo lógico
          ! pero se podría repensar
          !
  
          !
          ! substarct mu+mut 
          !
          do igaus=1,pgaus   
             gpvis_nsw(DEF_VECT,igaus) = gpvis_nsw(DEF_VECT,1) - (gpmut(DEF_VECT,igaus) + gpvis(DEF_VECT,igaus))
             !
             ! Limit so that it is not smaller than the laminar viscosity
             !
             gpvis_nsw(DEF_VECT,igaus) = max ( gpvis(DEF_VECT,igaus)  , gpvis_nsw(DEF_VECT,igaus) )
          end do

 
       end if


       
#ifdef OPENACC
    end do
#endif

  end subroutine ker_nsw_visc

end module mod_ker_nsw_visc
