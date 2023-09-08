!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_updtss.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Compute time step
!> @details Compute time step
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_updtss(itype_dtime)
  use  def_master
  use  def_domain
  use  def_parame
  use  def_nastal
  use  def_kermod
  use  mod_ker_proper 


  implicit none 

  integer(ip)      :: itype_dtime  !> either DT_PHYSICAL or DT_PSEUDO 
  integer(ip)      :: ielem,inode,ipoin,kpoin,pelty,pgaus,pnode,igaus,idime
  integer(ip)      :: icomp,i,iauxi
  real(rp)         :: velmo,dtmi1,dtmi2,timom(2),tiene(2),ticon(2),sound,rdumy(ndime,ndime)
  real(rp)         :: auxvi(ndofn_nsa) 

  real(rp)         :: hconv,hsoun,xvisc,xdens,xener,xenth,xtemp,xpres,cflcompute,cfl
     
  real(rp)         :: hleti,hmini,hmaxi,qufac,adgam
  real(rp)         :: elcod(ndime,mnode  )
  real(rp)         :: elvel(ndime,mnode),elpre(mnode),eltem(mnode),elden(mnode),elene(mnode), &
       elent(mnode),elumo(ndime,mnode)

  real(rp)         :: dvolu(mgaus)
  real(rp)         :: detjm,xshap(mgaus),dvelo,reate
  real(rp)         :: xjaci(ndime,ndime),xjacm(ndime,ndime),cartd(ndime,mnode,mgaus)
  
  real(rp)         :: tragl(ndime,ndime)                      ! Stabilization
  real(rp)         :: xvelo(ndime),xumom(ndime),chave(ndime,2),hleng(ndime),tauxi(2,5),xdtau, &
                      xsube(ndofn_nsa)

  real(rp)    :: veige,seige,prope_tmp(1),elwme(mnode),safety_factor(2), cfl_veloc
  real(rp)    :: aomp1, aomp2, aomp3,xmowe,rgasc,xhecp,xhecv,heats,dtbou(3),dtinv_phys_local

#ifdef _OPENMP
  INTEGER     :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
  INTEGER     :: thrid, numthr
  real(rp)    :: dtcri_nsaOMP(16,4),dtmax_nsaOMP(16,4)
#endif

  envol_nsa = 0.0_rp
  devol_nsa = 0.0_rp
  auxvi     = 0.0_rp
  elent     = 0.0_rp

  ! OOOOOJOOOOOOOO
!!  if (ittim == 20) safet_nsa = 5.0*safet_nsa
!!  if (ittim == 30) then
!!     safet_nsa = 2.0*safet_nsa
!!!     dtlim_nsa = 0.1_rp
!!  end if

  if( ittim > nisaf_nsa ) safet_nsa = min(safet_nsa*safex_nsa, safma_nsa)

  if(kfl_safet_table_nsa > 0) then
     safet_nsa = safet_table_nsa(1,1)
     if (ittim > 1) then
        do iauxi= 1,kfl_safet_table_nsa
           if (resid_nsa(4) < safet_table_nsa(2,iauxi)) then
              safet_nsa = safet_table_nsa(1,iauxi)           
              cycle
           end if
        end do
     end if
  end if

  if( INOTMASTER ) then

     safety_factor(DT_PHYSICAL)= safet_nsa
     safety_factor(DT_PSEUDO)  = safet_pseud_nsa

     dtmi1 = 0.0_rp
     dtmi2 = 0.0_rp
     
     if (kfl_tiext_nsa == 1) then  ! Externally fixed time step
        dtinv = dtext_nsa
        dtinv = 1.0_rp / dtinv
     else
        do ipoin = 1,npoin
!!!!!           dtieq_nsa(1:ndtdf_nsa,ipoin,2) = dtieq_nsa(1:ndtdf_nsa,ipoin,itype_dtime) 
           dtieq_nsa(1:ndtdf_nsa,ipoin,DT_PHYSICAL) = 0.0_rp
           dtieq_nsa(1:ndtdf_nsa,ipoin,DT_PSEUDO)   = 0.0_rp
        end do
     end if
     
     dtcri_nsa = 1.0_rp / zensa
     dtmax_nsa =          zensa
     
     icomp = TIME_N
     
!!!     if (itype_dtime == DT_PSEUDO) then ! Pseudo-time step is computed (called from: nsa_begite and nsa_endite(two))
!!!        safety_factor = safet_pseud_nsa
!!!        icomp = ITER_K
!!!     end if

#ifdef _OPENMP
     dtcri_nsaOMP = dtcri_nsa(1)
     dtmax_nsaOMP = dtmax_nsa(1)
#endif
!
! Loop over elements
!
!$*OMP   PARALLEL DO SCHEDULE (GUIDED)   &
!$*OMP   DEFAULT (NONE)                  &
!$*OMP   PRIVATE (alpre, bepre, cartd, chale, chave, detjm, dvelo, dvolu, &
!$*OMP            elcod, elden, elpre, eltem, elvel, elvis, elent, eplop, gamso, &
!$*OMP            hconv, hsoun, hleng, hleti, hmini, hmaxi, ielem, igaus, inode, &
!$*OMP            ipoin, kexvi, pelty, pgaus, pnode, qufac, reate, rdumy, &
!$*OMP            seige, sound, tauxi, ticon, tiene, timom, tragl, veaux, &
!$*OMP            veige, velmo, velop, velvi, xdens, xenth,xdenp, xdent, xdtau, &
!$*OMP            xjaci, xjacm, xpres, xshap, xtemp, xvelo, xvisc, aomp1, &
!$*OMP            aomp2, aomp3, thrid, i, itype_dtime, dtinv_phys_local) &
!$*OMP   SHARED  (adgam_nsa, coord, cpcoe_nsa, densi, dtcri_nsaOMP, dtmax_nsaOMP,           &
!$*OMP            dtieq_nsa, elmar, elwgt, hnatu, icomp, kfl_foreg_nsa,       &
!$*OMP            kfl_hconv_nsa, kfl_lopre_nsa, kfl_taufa_nsa, kfl_turbu_nsa, kfl_locti_nsa, & 
!$*OMP            kfl_visco_nsa, lnods, ltype, ndime, ndtdf_nsa, nelem, ngaus,&
!$*OMP            nnode, press, safet_nsa, tempe, veloc, visco,dtinv)

     elementary_loop: do ielem=1,nelem

#ifdef _OPENMP
        thrid = OMP_GET_THREAD_NUM() + 1
#endif
 

        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        
        qufac = 1.0_rp
        if((ndime==2).and.(pnode>=4)) qufac = 0.5_rp
        if((ndime==3).and.(pnode>=5)) qufac = 0.5_rp
        
        xenth = 0.0_rp
        xdens = 0.0_rp
        xener = 0.0_rp
        xtemp = 0.0_rp
        xpres = 0.0_rp
        chave = 0.0_rp
        xvelo = 0.0_rp
        dvelo = 0.0_rp
        
        if (kfl_relat_nsa == 1) then
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              elent(inode)  = entha(ipoin,icomp)
           end do
        end if
        
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do idime=1,ndime
              elvel(idime,inode)  = veloc(idime,ipoin,icomp)
              elcod(idime,inode)  = coord(idime,ipoin)
              elumo(idime,inode)  = umome(idime,ipoin,1)
           end do
           elden(inode)  = densi(ipoin,icomp)
           !elvis(inode)  = densi(ipoin,icomp)
           elene(inode)  = energ(ipoin,icomp)
           elpre(inode)  = press(ipoin,icomp)
           eltem(inode)  = tempe(ipoin,icomp)
           
           if (kfl_coupl(ID_NASTAL,ID_CHEMIC) >= 1 ) then
              elwme(inode)  = wmean(ipoin,1)
           else
              elwme(inode)  = mowei_nsa
           endif
        end do
        
        ! hleng and tragl at center of gravity
        call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)
        
        hmini = hleng(ndime)  ! hleng(ndime) is the smallest
        hmaxi = hleng(1)      ! hleng(1) is the largest
        hconv = hmini
        hsoun = hmini
        hleti = hmini*hmaxi           
        if (kfl_higha_nsa == 1) hleti= hmini*hmini

        do igaus=1,pgaus
           call elmder(pnode,ndime,elmar(pelty)%deriv(1,1,igaus),elcod,cartd(1,1,igaus),detjm,xjacm,xjaci)
           dvolu(igaus)=elmar(pelty)%weigp(igaus)*detjm        
           xdens = 0.0_rp
           xener = 0.0_rp 
           xenth = 0.0_rp
           xtemp = 0.0_rp
           xpres = 0.0_rp
           xvisc = 0.0_rp
           xmowe = 0.0_rp
           xhecp = 0.0_rp
           xhecv = 0.0_rp
           do idime=1,ndime
              xvelo(idime) = 0.0_rp
              xumom(idime) = 0.0_rp
              xsube(idime) = umosg_nsa(idime,ielem,igaus,1)
              if (kfl_track_nsa == 0) xsube(idime) = 0.0_rp
           end do
           xsube(ndime+1) = densg_nsa(ielem,igaus,1)
           xsube(ndime+2) = enesg_nsa(ielem,igaus,1)
           if (kfl_track_nsa == 0) then
              xsube(ndime+1) = 0.0_rp
              xsube(ndime+2) = 0.0_rp
           end if
           dvelo = 0.0_rp  
           
           do inode=1,pnode
              xshap(igaus) = elmar(pelty)%shape(inode,igaus)
              xdens = xdens + elden(inode) * xshap(igaus)
              xener = xener + elene(inode) * xshap(igaus)
              xenth = xenth + elent(inode) * xshap(igaus)
              xtemp = xtemp + eltem(inode) * xshap(igaus)
              xpres = xpres + elpre(inode) * xshap(igaus)
              do idime=1,ndime
                 xvelo(idime) = xvelo(idime) + elvel(idime,inode)*xshap(igaus)
                 xumom(idime) = xumom(idime) + elumo(idime,inode)*xshap(igaus)
                 dvelo        = dvelo + cartd(idime,inode,igaus) * elvel(idime,inode)
              end do
              xmowe = xmowe + elwme(inode) * xshap(igaus)  ! Molecular weight at gauss points
           end do
           
           velmo = 0.0_rp 
           do idime=1,ndime
              ! Velocity module
              velmo = velmo + xvelo(idime)*xvelo(idime)
           end do
           
           velmo = sqrt(velmo)
           
           !
           ! Compute the laminar viscosity
           !
           if (kfl_prope /= 0 ) then
              call ker_proper('VISCO','IGAUS',1_ip,ielem,prope_tmp,pnode,pgaus,elmar(pelty)%shape(:,igaus),cartd(:,:,igaus)) 
              xvisc = prope_tmp(1)
              call ker_proper('SPHEA','IGAUS',1_ip,ielem,prope_tmp,pnode,pgaus,elmar(pelty)%shape(:,igaus),cartd(:,:,igaus))
              xhecp = prope_tmp(1)
              rgasc = runiv_nsa / xmowe  
           else 
              call nsa_lawvis(-1,icomp,xvisc,xtemp,rdumy(ndime,ndime))
              xhecp = cpcoe_nsa
              rgasc = runiv_nsa / xmowe
           endif
           
           xhecv = xhecp - rgasc
           
           ! Compute the sound speed
           
           sound = 0.0_rp
           
           adgam = xhecp / (xhecp - rgasc)
           
           sound= sqrt(adgam * xpres / (xdens+xsube(ndime+1)))
           
           ! Compute the different time increments (MUY INCOMPLETO TODAVIA...!!!!!)
           
           seige = sound
           
           if (velmo < zensa) then           
              if (ittim <= 0 .and. kfl_zero_initial_velocity_nsa == 1 .and. kfl_lopre_nsa > 1) then
                 velmo = seige
              else
                 velmo = zensa
              end if
           end if
           
           veige = velmo
           
           reate = (1.0_rp + rgasc / xhecv)*dvelo
           
           if (kfl_coupl(ID_NASTAL,ID_CHEMIC) >= 1 ) then
              heats = chemical_heat(ielem)%a(igaus)
              reate = reate + heats / xener
           endif
           
           if (dvelo <= 0.0_rp) then
              reate= -dvelo
              ! reate= 0.0_rp
           end if

           if (kfl_reate_nsa == 0) reate=0.0_rp

           tauxi = 0.0_rp           
           xdtau = xdens
           if (kfl_taufa_nsa(1,1) == 1) tauxi(DT_PHYSICAL,1) = veige/hconv
           if (kfl_taufa_nsa(2,1) == 1) tauxi(DT_PHYSICAL,2) = seige/hsoun
           if (kfl_taufa_nsa(3,1) == 1) tauxi(DT_PHYSICAL,3) = reate 
           if (kfl_taufa_nsa(4,1) == 1) tauxi(DT_PHYSICAL,4) = 4.0_rp*xvisc/hleti/xdtau

           tauxi(DT_PSEUDO,1) = veige/hconv
           tauxi(DT_PSEUDO,2) = tauxi(DT_PHYSICAL,2)
           tauxi(DT_PSEUDO,3) = tauxi(DT_PHYSICAL,3) 
                      
           if (kfl_zevel_nsa(1) == 0) then
              !
              ! velocity is non-zero, so eliminate sound from the pseudo-dt computation
              !
              if (veige .gt. 1.0e-10) then
                 tauxi(DT_PSEUDO,2) = 0.0_rp 
                 tauxi(DT_PSEUDO,3) = 0.0_rp
              end if
              if (kfl_lopre_nsa == 3) then ! CHOI & MERKLE PRECONDITIONER IS APPLIED
                 if (kfl_pseud_nsa == 1 .and. kfl_taufa_nsa(1,1) == 1) then
                    cfl_veloc = velmo / max(dtinv,1.0e-12_rp) / hconv
!                    dtinv_phys_local=  (tauxi(DT_PHYSICAL,1) + tauxi(DT_PHYSICAL,2) + tauxi(DT_PHYSICAL,4))/qufac
!                    cfl_veloc = velmo / dtinv_phys_local / hconv * safety_factor(DT_PHYSICAL)

                    if (cfl_veloc < zensa) cfl_veloc = zensa
                    cfl_veloc = sqrt(1.0_rp + 1.0_rp / cfl_veloc / cfl_veloc)
                    
                    tauxi(DT_PSEUDO,1) = tauxi(DT_PSEUDO,1) * cfl_veloc / 2.0_rp
                 end if
              end if              
              if (kfl_lopre_nsa > 1 .and. kfl_pseud_nsa == 0) then 
                 tauxi(DT_PHYSICAL,2) = 0.0_rp
                 tauxi(DT_PHYSICAL,3) = 0.0_rp
              end if
           end if


!!           if (kfl_lopre_nsa > 1 .and. kfl_pseud_nsa == 0) then 
!!              if (kfl_zevel_nsa(1) == 0) then
!!                 tauxi(DT_PHYSICAL,2) = 0.0_rp
!!                 tauxi(DT_PHYSICAL,3) = 0.0_rp
!!              end if
!!            !!tauxi(DT_PHYSICAL,4) = 0.0_rp 
!!           end if
!!           !           if (itype_dtime == DT_PSEUDO) then ! Pseudo-time step is computed
!!           if (kfl_lopre_nsa > 1) then 
!!              if (kfl_zevel_nsa(1) == 0) then ! velocity is different than zero somewhere
!!                 tauxi(DT_PHYSICAL,2) = 0.0_rp
!!                 tauxi(DT_PHYSICAL,3) = 0.0_rp
!!                 if (kfl_lopre_nsa == 3) then ! CHOI & MERKLE PRECONDITIONER IS APPLIED
!!                    if (kfl_pseud_nsa == 1 .and. kfl_taufa_nsa(1,1) == 1) then
!!                       cfl_veloc = velmo / dtinv / hconv
!!                       tauxi(DT_PHYSICAL,1) = tauxi(DT_PHYSICAL,1) * sqrt(1.0_rp + 1.0_rp / cfl_veloc / cfl_veloc) / 2.0_rp
!!                    end if
!!                 end if
!!              end if
!!           end if
!!           !           end if

           
           do i=1,2
              !! timom(i) = tauxi(i,1) + tauxi(i,2) + tauxi(i,3) + tauxi(i,4)
              timom(i) = tauxi(i,1) + tauxi(i,2) + tauxi(i,4)
              timom(i) = qufac / timom(i)          
              tiene(i)= timom(i)
              ticon(i)= timom(i)
           end do
           ! compute volume integral of extensive quantities
           envol_nsa = envol_nsa + xener * dvolu(igaus)
           devol_nsa = devol_nsa + xdens * dvolu(igaus)
           
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              
              do i=1,2
                 aomp1 = timom(i) * dvolu(igaus) * xshap(igaus) / vmass(ipoin) 
                 aomp2 = ticon(i) * dvolu(igaus) * xshap(igaus) / vmass(ipoin) 
                 aomp3 = tiene(i) * dvolu(igaus) * xshap(igaus) / vmass(ipoin) 
                 
                 if (kfl_cfllo_nsa == 1) then
                    cfl = cflcompute(veige/seige)
                    aomp1 = aomp1 * cfl
                    aomp2 = aomp2 * cfl
                    aomp3 = aomp3 * cfl
                 end if
                 
                 !$*OMP      CRITICAL (lock_dtieq_nsa)
                 dtieq_nsa(1,ipoin,i) = dtieq_nsa(1,ipoin,i) + aomp1
                 dtieq_nsa(2,ipoin,i) = dtieq_nsa(2,ipoin,i) + aomp2
                 dtieq_nsa(3,ipoin,i) = dtieq_nsa(3,ipoin,i) + aomp3
                 !$*OMP      END CRITICAL (lock_dtieq_nsa)
              end do

           end do
           
        end do
        
     end do elementary_loop
    
    !
    ! Sum up nodal time steps
    !
    call nsa_parall(eight)
    
    !
    ! Loop over the nodes to compute dtcri_nsa
    !
    nodal_loop: do ipoin = 1,npoin   
       
       dtieq_nsa(1,ipoin,DT_PHYSICAL) = dtieq_nsa(1,ipoin,DT_PHYSICAL) * safety_factor(DT_PHYSICAL)
       dtieq_nsa(2,ipoin,DT_PHYSICAL) = dtieq_nsa(2,ipoin,DT_PHYSICAL) * safety_factor(DT_PHYSICAL)
       dtieq_nsa(3,ipoin,DT_PHYSICAL) = dtieq_nsa(3,ipoin,DT_PHYSICAL) * safety_factor(DT_PHYSICAL)

       if (dtieq_nsa(1,ipoin,DT_PHYSICAL) < dtcri_nsa(1)) dtcri_nsa(1) = dtieq_nsa(1,ipoin,DT_PHYSICAL)
       if (dtieq_nsa(2,ipoin,DT_PHYSICAL) < dtcri_nsa(2)) dtcri_nsa(2) = dtieq_nsa(2,ipoin,DT_PHYSICAL)
       if ((ndtdf_nsa == 3) .and. (dtieq_nsa(3,ipoin,DT_PHYSICAL) < dtcri_nsa(3))) &
            dtcri_nsa(3) = dtieq_nsa(3,ipoin,DT_PHYSICAL)

       if (dtieq_nsa(1,ipoin,DT_PHYSICAL) > dtmax_nsa(1)) dtmax_nsa(1) = dtieq_nsa(1,ipoin,DT_PHYSICAL)
       if (dtieq_nsa(2,ipoin,DT_PHYSICAL) > dtmax_nsa(2)) dtmax_nsa(2) = dtieq_nsa(2,ipoin,DT_PHYSICAL)
       if ((ndtdf_nsa == 3) .and. (dtieq_nsa(3,ipoin,DT_PHYSICAL) > dtmax_nsa(3))) &
            dtmax_nsa(3) = dtieq_nsa(3,ipoin,DT_PHYSICAL)

       dtieq_nsa(1,ipoin,DT_PSEUDO) = dtieq_nsa(1,ipoin,DT_PSEUDO) * safety_factor(DT_PSEUDO)
       dtieq_nsa(2,ipoin,DT_PSEUDO) = dtieq_nsa(2,ipoin,DT_PSEUDO) * safety_factor(DT_PSEUDO)
       dtieq_nsa(3,ipoin,DT_PSEUDO) = dtieq_nsa(3,ipoin,DT_PSEUDO) * safety_factor(DT_PSEUDO)
       
    end do nodal_loop

 end if

 
 !
 ! For the parallel case, compute the all-mesh minimum
 !
 call pararr('MIN',0_ip,ndtdf_nsa,dtcri_nsa)
 
 if (INOTMASTER) then
    ! Global time step corresponds to dtlim_nsa = 0.0_rp
    ! Fully Local time step corresponds to dtlim_nsa = 1.0_rp
    dtbou(1) = dtcri_nsa(1) + dtlim_nsa * (dtmax_nsa(1) - dtcri_nsa(1))
    dtbou(2) = dtcri_nsa(2) + dtlim_nsa * (dtmax_nsa(2) - dtcri_nsa(2))
    dtbou(3) = dtcri_nsa(3) + dtlim_nsa * (dtmax_nsa(3) - dtcri_nsa(3))
    nodal_loop_2: do ipoin = 1,npoin   
       if (dtieq_nsa(1,ipoin,DT_PHYSICAL) > dtbou(1)) dtieq_nsa(1,ipoin,DT_PHYSICAL) = dtbou(1)
       if (dtieq_nsa(2,ipoin,DT_PHYSICAL) > dtbou(2)) dtieq_nsa(2,ipoin,DT_PHYSICAL) = dtbou(2)
       if (dtieq_nsa(3,ipoin,DT_PHYSICAL) > dtbou(3)) dtieq_nsa(3,ipoin,DT_PHYSICAL) = dtbou(3)
    end do nodal_loop_2
 end if
 
 dtinv_nsa = dtcri_nsa(1) 
 if (dtcri_nsa(2) < dtinv_nsa) dtinv_nsa = dtcri_nsa(2)
 if (dtcri_nsa(3) < dtinv_nsa) dtinv_nsa = dtcri_nsa(3)
 
 dtinv_nsa = 1.0_rp / dtinv_nsa
 
 if (kfl_timco == 1) then
    dtinv=max(dtinv,dtinv_nsa)
 end if
 
 if( INOTMASTER ) then     
    if (kfl_tiext_nsa == 1) then  ! Externally fixed time step
       dtinv_nsa = dtinv
    end if
    ! For the first iter
!    if (itinn(modul)==0) then
!       do ipoin = 1,npoin
!          dtieq_nsa(1:ndtdf_nsa,ipoin,DT_FIRST_ITER) = 1.0_rp/dtinv_nsa
!       end do
!    end if

 end if
 
end subroutine nsa_updtss


function cflcompute(xmach)
    use  def_master
    use  def_nastal 

    real(rp)  :: xmach,cflcompute,cmax,cmin
    integer(ip) :: icfl

    cmin= fcfll_nsa(        1,1)
    cmax= fcfll_nsa(mcfll_nsa,1)


    if (xmach > cmax) then
       cflcompute= fcfll_nsa(mcfll_nsa,2)
       
    else if (xmach < cmin) then
       cflcompute= fcfll_nsa(        1,2)

    else
      ! ! add one because the first value corresponds to zero
      ! icfl = 1_ip + int(real(mcfll_nsa) * (xmach - cmin)/(cmax - cmin))
       
      ! cmin= fcfll_nsa(icfl-1,1)
      ! cmax= fcfll_nsa(icfl  ,1)
       
      ! cflcompute= fcfll_nsa(icfl-1,2) + (fcfll_nsa(icfl,2)-fcfll_nsa(icfl-1,2))*(xmach - cmin)/(cmax - cmin)
      
       icfl = 1_ip
       do 
          if (xmach < fcfll_nsa(icfl,1)) then
             cflcompute = (fcfll_nsa(icfl,2)) 
             exit
          else 
             icfl=icfl+1_ip
          endif
       end do
    end if
    
  end function cflcompute
