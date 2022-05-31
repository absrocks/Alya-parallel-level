!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_updtss_cdr.f90
!> @author  Margarida Moragues
!> @date    16/11/1966
!> @brief   Compute time step size for the cdr-in-nastal
!> @details Compute time step size for the cdr-in-nastal
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_updtss_cdr
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_updtss
  ! NAME 
  !    nsa_updtss
  ! DESCRIPTION
  !    This routine computes the time step size. It is computed at the
  !    element center of gravity. 
  !    The time step is computed looping over the elements, but the local
  !    time step is assigned on a per-node basis. Each node takes the 
  !    smallest time step of all those corresponding to the elements 
  !    containing this node.
  !
  ! USED BY
  !    nsa_begste
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_parame

  use      def_nastal
  use      mod_communications, only : PAR_MIN

  implicit none 

  integer(ip)      :: ielem,inode,ipoin,kpoin,pelty,pgaus,pnode,igaus,ittot,idofn,idime
  integer(ip)      :: iequa,icomp
  real(rp)         :: velmo,dtmi1,dtmi2,timom,tiene,ticon,sound,rdumy(ndime,ndime),auxvi(ndofn_nsa)
  real(rp)         :: hconv,htran,xvisc,xdens,xener,xenth,xtemp,xpres
     
  real(rp)         :: hleti,hmini,qufac
  real(rp)         :: elcod(ndime,mnode  )
  real(rp)         :: elvel(ndime,mnode),elpre(mnode),eltem(mnode),elden(mnode),elene(mnode), &
       elvis(mnode),elent(mnode),elumo(ndime,mnode)

  real(rp)         :: dvolu(mgaus)
  real(rp)              :: detjm,xshap(mgaus),dvelo,reate
  real(rp)              :: xjaci(ndime,ndime),xjacm(ndime,ndime),cartd(ndime,mnode,mgaus)
  
  real(rp)         :: tragl(ndime,ndime)                      ! Stabilization
  real(rp)         :: xvelo(ndime),xumom(ndime),chave(ndime,2),hleng(ndime),tauxi(2,5),xdtau, &
       xsube(ndofn_nsa)

  real(rp), target :: dttar(10)

  real(rp)    :: veige,seige
  real(rp)    :: aomp1, aomp2, aomp3

#ifdef _OPENMP
  INTEGER     :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
  INTEGER     :: thrid, numthr
!  real(rp)    :: dtcri_nsaOMP(3,4)
  real(rp)    :: dtcri_nsaOMP(16,4)
#endif

  if( INOTMASTER ) then

     dtmi1 = 0.0_rp
     dtmi2 = 0.0_rp
     auxvi = 0.0_rp

     if (kfl_tiext_nsa == 1 ) then                             ! Externally fixed time step
        do ipoin = 1,npoin
           dtieq_nsa(1:ndtdf_nsa,ipoin,2) = dtieq_nsa(1:ndtdf_nsa,ipoin,1) 
           dtieq_nsa(1:ndtdf_nsa,ipoin,1) = dtext_nsa
        end do
        dtcri_nsa = dtext_nsa
     else
        do ipoin = 1,npoin
           dtieq_nsa(1:ndtdf_nsa,ipoin,2) = dtieq_nsa(1:ndtdf_nsa,ipoin,1) 
           dtieq_nsa(1:ndtdf_nsa,ipoin,1) = 0.0_rp
        end do
        dtcri_nsa = 1000000000_rp
     end if

     !
     ! Loop over elements
     !
     icomp=min(3_ip,ncomp_nsa)

#ifdef _OPENMP
     dtcri_nsaOMP = dtcri_nsa(1)
#endif


!
! Loop over elements
!

!$*OMP   PARALLEL DO SCHEDULE (GUIDED)   &
!$*OMP   DEFAULT (NONE)                  &
!$*OMP   PRIVATE (alpre, bepre, cartd, chale, chave, detjm, dvelo, dvolu, &
!$*OMP            elcod, elden, elpre, eltem, elvel, elvis, elent, eplop, gamso, &
!$*OMP            hconv, hleng, hleti, hmini, htran, ielem, igaus, inode, &
!$*OMP            ipoin, kexvi, pelty, pgaus, pnode, qufac, reate, rdumy, &
!$*OMP            seige, sound, tauxi, ticon, tiene, timom, tragl, veaux, &
!$*OMP            veige, velmo, velop, velvi, xdens, xenth,xdenp, xdent, xdtau, &
!$*OMP            xjaci, xjacm, xpres, xshap, xtemp, xvelo, xvisc, aomp1, &
!$*OMP            aomp2, aomp3, thrid) &
!$*OMP   SHARED  (adgam_nsa, coord, cpcoe_nsa, densi, dtcri_nsaOMP,           &
!$*OMP            dtieq_nsa, elmar, elwgt, hnatu, icomp, kfl_foreg_nsa,       &
!$*OMP            kfl_hconv_nsa, kfl_lopre_nsa, kfl_taufa_nsa, kfl_turbu_nsa, kfl_locti_nsa, & 
!$*OMP            kfl_visco_nsa, lnods, ltype, ndime, ndtdf_nsa, nelem, ngaus,&
!$*OMP            nnode, press, safet_nsa, tempe, veloc, visco)


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
        xvisc = 0.0_rp
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
              elumo(idime,inode) = umome(idime,ipoin,1)
           end do
           elden(inode)  = densi(ipoin,icomp)
           elene(inode) = energ(ipoin,icomp)
           elpre(inode)  = press(ipoin,icomp)
           eltem(inode)  = tempe(ipoin,icomp)
           elvis(inode)  = visco(ipoin,icomp)
        end do

        ! hleng and tragl at center of gravity
        call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)
        
        htran= hleng(    1)
        hconv= hleng(ndime)
        hmini= hleng(ndime)  ! hleng(ndime) is the smallest
        hleti = htran*htran           
        
        do igaus=1,pgaus
           call elmder(pnode,ndime,elmar(pelty)%deriv(1,1,igaus),elcod,cartd(1,1,igaus),detjm,xjacm,xjaci)
           dvolu(igaus)=elmar(pelty)%weigp(igaus)*detjm        
           xdens = 0.0_rp
           xener = 0.0_rp 
           xenth = 0.0_rp
           xtemp = 0.0_rp
           xpres = 0.0_rp
           xvisc = 0.0_rp
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
             ! xpres = xpres + elpre(inode) * xshap(igaus)
              xvisc = xvisc + elvis(inode) * xshap(igaus)
              do idime=1,ndime
                 !xvelo(idime) = xvelo(idime) + elvel(idime,inode)*xshap(igaus)
                 xumom(idime) = xumom(idime) + elumo(idime,inode)*xshap(igaus)
                 dvelo      = dvelo + cartd(idime,inode,igaus) * elvel(idime,inode)
              end do
           end do
           velmo = 0.0_rp 
           do idime=1,ndime
              xvelo(idime) = (xumom(idime)+xsube(idime))/(xdens+xsube(ndime+1))
              ! Velocity module
              velmo = velmo + xvelo(idime)*xvelo(idime)
           end do
           if (kfl_unkse_nsa < 10 ) then
              xpres = rgasc_nsa * (xener+xsube(ndime+2)-0.5_rp*(xdens+xsube(ndime+1))*velmo) / cvcoe_nsa
           else if (kfl_unkse_nsa < 20 ) then
              xpres = rgasc_nsa * (xdens+xsube(ndime+1)) * (xtemp+xsube(ndime+2))
           end if

           velmo= sqrt(velmo)
           
           if (velmo > zensa) then           
              !!           ! Compute the characteristic length chale if needed
              !!           if (kfl_hconv_nsa == 1) then
              !!              chave(1:ndime,1) = xvelo(1:ndime)
              !!              call nsa_elmchl(tragl,hleng,chave,chale,ndime,&
              !!                   pnode,hnatu(pelty),1_ip,zensa)
              !!              hconv = chale(1)
              !!           end if
           else
              
              velmo = zensa 
              
           end if
        
           ! Compute the viscosity
           
           call nsa_lawvis(one,one,xvisc,xtemp,rdumy(ndime,ndime))
           
           ! Compute the sound speed
           
           sound=0.0
           !        if (kfl_foreg_nsa == 0) sound= sofac_nsa * sqrt(xtemp)
           if (kfl_foreg_nsa == 0) sound= sqrt(adgam_nsa * xpres / (xdens+xsube(ndime+1)))
           if (kfl_relat_nsa == 1) sound= sqrt((1.0_rp - 1.0_rp / xenth)*(adgam_nsa-1.0_rp))

           
           ! Compute the different time increments (MUY INCOMPLETO TODAVIA...!!!!!)
           
           seige = sound
           veige = velmo
           
           reate = (1.0_rp + rgasc_nsa / cvcoe_nsa)*dvelo
           if (dvelo <= 0.0_rp) then
              reate= -dvelo
              !           reate= 0.0_rp
           end if
           
           if (kfl_relat_nsa == 1) then
              !           tauxi(1,1) = (veige + seige) / (1.0_rp + veige*seige)
              !           tauxi(1,2) = (veige - seige) / (1.0_rp - veige*seige)           
              tauxi(1,1) = abs((veige*(1.0_rp - seige*seige) + seige*(1.0_rp-veige*veige)) &
                   / (1.0_rp - veige*seige*veige*seige))
              tauxi(1,2) = abs((veige*(1.0_rp - seige*seige) - seige*(1.0_rp-veige*veige)) &
                   / (1.0_rp - veige*seige*veige*seige))
              seige = tauxi(1,1)
              
              if (tauxi(1,2) > seige) seige = tauxi(1,2)
              seige= 1.0_rp / seige
              
              veige = 0.0_rp
              reate = 0.0_rp
           end if
           
           tauxi = 0.0_rp
           
           xdtau = xdens
           if (kfl_taufa_nsa(1,1) == 1) tauxi(1,1) = veige/hconv
           if (kfl_taufa_nsa(2,1) == 1) tauxi(1,2) = seige/hmini
           if (kfl_taufa_nsa(3,1) == 1) tauxi(1,3) = reate
           if (kfl_taufa_nsa(4,1) == 1) tauxi(1,4) = 4.0*xvisc/hleti/xdtau
           
           !                tauxi(1,1) = tauxi(1,1)*tauxi(1,1) 
           !                tauxi(1,2) = tauxi(1,2)*tauxi(1,2) 
           !                tauxi(1,3) = tauxi(1,3)*tauxi(1,3) 
           !                tauxi(1,4) = tauxi(1,4)*tauxi(1,4) 
           
           timom = tauxi(1,1) + tauxi(1,2) + tauxi(1,3) + tauxi(1,4)
           !        timom = sqrt(timom)
           
           timom = qufac / timom
           
           tiene= timom
           ticon= timom

           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              
              aomp1 = timom * dvolu(igaus) * xshap(igaus) / vmass(ipoin) 
              aomp2 = ticon * dvolu(igaus) * xshap(igaus) / vmass(ipoin) 
              aomp3 = tiene * dvolu(igaus) * xshap(igaus) / vmass(ipoin) 
              
              !$*OMP      CRITICAL (lock_dtieq_nsa)
              dtieq_nsa(1,ipoin,1) = dtieq_nsa(1,ipoin,1) + aomp1
              if (kfl_foreg_nsa == 0) then
                 dtieq_nsa(2,ipoin,1) = dtieq_nsa(2,ipoin,1) + aomp2
                 dtieq_nsa(3,ipoin,1) = dtieq_nsa(3,ipoin,1) + aomp3
              end if
              !$*OMP      END CRITICAL (lock_dtieq_nsa)
              
           end do

        end do

     end do elementary_loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     elementary_loop_2: do ielem=1,nelem
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do idime=1,ndime
              elcod(idime,inode) = coord(idime,ipoin)
           end do
        end do

        ! hleng and tragl at center of gravity
        call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)        
        htran= hleng(    1)
        hconv= hleng(ndime)
        hmini= hleng(ndime)  ! hleng(ndime) is the smallest
        hleti = htran*htran           

        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           dtieq_nsa(1,ipoin,1) = 0.0_rp
           do idime=1,ndime
              dtieq_nsa(1,ipoin,1) = dtieq_nsa(1,ipoin,1) + conve_nsa(idime) * conve_nsa(idime)    
           end do
           dtieq_nsa(1,ipoin,1) = sqrt(dtieq_nsa(1,ipoin,1)) / hconv
           dtieq_nsa(1,ipoin,1) = dtieq_nsa(1,ipoin,1) + react_nsa
           dtieq_nsa(1,ipoin,1) = dtieq_nsa(1,ipoin,1) + diffu_nsa / hleti
           dtieq_nsa(1,ipoin,1) = 1.0_rp / dtieq_nsa(1,ipoin,1)
           if (kfl_foreg_nsa == 0) then
              dtieq_nsa(2,ipoin,1) = dtieq_nsa(1,ipoin,1)
              dtieq_nsa(3,ipoin,1) = dtieq_nsa(1,ipoin,1)
           end if
        end do
     end do elementary_loop_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     !
     ! Sum up nodal time steps
     !
     call nsa_parall(eight)

     nodal_loop: do ipoin = 1,npoin   

        if (dtieq_nsa(1,ipoin,1) < dtcri_nsa(1)) dtcri_nsa(1) = dtieq_nsa(1,ipoin,1)
        if (dtieq_nsa(2,ipoin,1) < dtcri_nsa(2)) dtcri_nsa(2) = dtieq_nsa(2,ipoin,1)
        if ((ndtdf_nsa == 3) .and. (dtieq_nsa(3,ipoin,1) < dtcri_nsa(3))) dtcri_nsa(3) = dtieq_nsa(3,ipoin,1)

     end do nodal_loop

     dtcri_nsa(1) = dtcri_nsa(1) * safet_nsa
     dtcri_nsa(2) = dtcri_nsa(2) * safet_nsa
     dtcri_nsa(3) = dtcri_nsa(3) * safet_nsa
     
  end if
  !
  ! Look for minimum over whole mesh
  !
  if( IPARALL ) call PAR_MIN(ndtdf_nsa,dtcri_nsa)
  
  do iequa= 1,ndtdf_nsa
     dtinv_nsa  = 1.0_rp/(dtcri_nsa(iequa))  
     if(kfl_timco==1) dtinv=max(dtinv,dtinv_nsa)
  end do

  if( INOTMASTER ) then     
     !
     ! Assign 1/dt
     !
     ! For the first time step, dt^(n-1) = dt^n
     if (itinn(modul)==0) then
        do ipoin = 1,npoin
           dtieq_nsa(1:ndtdf_nsa,ipoin,2) = dtieq_nsa(1:ndtdf_nsa,ipoin,1) 
        end do
     end if

  end if


end subroutine nsa_updtss_cdr
