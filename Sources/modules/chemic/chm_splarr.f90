subroutine chm_splarr()
  !-----------------------------------------------------------------------
  !****f* domain/chm_splarr
  ! NAME
  !    chm_splarr
  ! DESCRIPTION
  !    This routines calculates the diagonal mass matrix using a 
  !    closed rule
  ! OUTPUT
  !    VMASS_CHM(NPOIN) : Diagonal mass matrix
  ! USED BY
  !    Domain
  !*** 
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  use def_chemic
  use mod_memchk
  implicit none
  integer(ip)   :: ipoin,inode,igaus,ielem,izdom,jzdom
  integer(ip)   :: pgaus,pnode,pelty
  integer(4)    :: istat
  real(rp)      :: elcod(ndime,mnode)
  real(rp)      :: gpdet,gpvol,xjacm(9) 
  real(rp)      :: numer(5),denom,volum,xfact,aux

  if( kfl_assem_chm >= 2 .and. INOTMASTER ) then

     !-------------------------------------------------------------------
     !
     ! IDIMA_CHM: Look for diagonal position
     !
     !-------------------------------------------------------------------

     if( solve(1) % kfl_symme == 0 ) then
        allocate(idima_chm(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'IDIMA_CHM','chm_splarr',idima_chm)
        do ipoin = 1,npoin
           izdom = r_dom(ipoin)-1
           do while( izdom < r_dom(ipoin+1) - 1 )
              izdom = izdom +1 
              if( c_dom(izdom) == ipoin ) then
                 jzdom = izdom
                 izdom = r_dom(ipoin+1) -1
              end if
           end do
           idima_chm(ipoin) = jzdom
        end do
     end if

     !-------------------------------------------------------------------
     !
     ! SMATR_CHM: Save matrix
     !
     !-------------------------------------------------------------------

     if( kfl_assem_chm == 3 ) then
        if( kfl_coupl_chm == 0 ) then
           allocate(smatr_chm(nzmat * nclas_chm),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'SMATR_CHM','chm_splarr',smatr_chm)
           allocate(shsid_chm(nzrhs * nclas_chm),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'SHSID_CHM','chm_splarr',shsid_chm)
        else
           call runend('CHM_SLPARR: JACOBI WITH SPLIT MATRICES NOT CODED')
        end if
     end if

     !-------------------------------------------------------------------
     !
     ! VMASS_CHM: Mass matrix
     !
     !-------------------------------------------------------------------

     allocate(vmass_chm(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VMASS_CHM','chm_splarr',vmass_chm)
     !
     ! Loop over elements
     !
     if( ndime == 3 .and. lexis(PYR05) == 0 .and. kfl_naxis == 0 .and. kfl_spher == 0 ) then
        !
        ! 3D: Avoid if and go fast
        !
        !$OMP  PARALLEL DO SCHEDULE (GUIDED) & 
        !$OMP  DEFAULT (NONE)                                                    &
        !$OMP  PRIVATE ( aux, denom, elcod, gpdet, gpvol,                        &
        !$OMP            ielem, igaus, inode, ipoin, numer, pelty, pgaus, pnode, &
        !$OMP            volum, xjacm, xfact                                   ) &
        !$OMP  SHARED  ( coord, elmar, lnods, ltype, kfl_spher, kfl_naxis,       &
        !$OMP            nelem, ngaus, nnode,                                    &
#ifndef NDIMEPAR
        !$OMP            ndime,                                                  &
#endif
        !$OMP            vmass_chm                                             )
        !
        do ielem = 1,nelem
           pelty = ltype(ielem) 
           pnode = nnode(pelty)
           do inode = 1,pnode
              ipoin          = lnods(inode,ielem)
              elcod(1,inode) = coord(1,ipoin)
              elcod(2,inode) = coord(2,ipoin)
              elcod(3,inode) = coord(3,ipoin)
           end do
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              call jacdet(&
                   ndime,pnode,elcod,elmar(pelty)%deric(1,1,inode),&
                   xjacm,gpdet)
              gpvol = elmar(pelty)%weigc(inode)*gpdet
              !$OMP         ATOMIC     
              vmass_chm(ipoin) = vmass_chm(ipoin) + gpvol
           end do
        end do
        !$OMP END PARALLEL DO

     else if( ndime == 2 .and. kfl_naxis == 0 .and. kfl_spher == 0 ) then
        !
        ! 2D: Avoid if and go fast
        !
        !$OMP  PARALLEL DO SCHEDULE (GUIDED) & 
        !$OMP  DEFAULT (NONE)                                                    &
        !$OMP  PRIVATE ( aux, denom, elcod, gpdet, gpvol,                        &
        !$OMP            ielem, igaus, inode, ipoin, numer, pelty, pgaus, pnode, &
        !$OMP            volum, xjacm, xfact                                   ) &
        !$OMP  SHARED  ( coord, elmar, lnods, ltype, kfl_naxis, kfl_spher,       &
        !$OMP            nelem, ngaus, nnode,                                    &
#ifndef NDIMEPAR
        !$OMP            ndime,                                                  &
#endif
        !$OMP            vmass_chm                                             )
        !
        do ielem = 1,nelem
           pelty = ltype(ielem) 
           pnode = nnode(pelty)
           do inode = 1,pnode
              ipoin          = lnods(inode,ielem)
              elcod(1,inode) = coord(1,ipoin)
              elcod(2,inode) = coord(2,ipoin)
           end do
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              call jacdet(&
                   ndime,pnode,elcod,elmar(pelty)%deric(1,1,inode),&
                   xjacm,gpdet)
              gpvol = elmar(pelty)%weigc(inode)*gpdet
              !$OMP         ATOMIC     
              vmass_chm(ipoin) = vmass_chm(ipoin) + gpvol
           end do
        end do
        !$OMP END PARALLEL DO

     else        
        !
        ! Cannot avoid if
        !
        !$OMP  PARALLEL DO SCHEDULE (GUIDED) & 
        !$OMP  DEFAULT (NONE)                                                    &
        !$OMP  PRIVATE ( aux, denom, elcod, gpdet, gpvol,                        &
        !$OMP            ielem, igaus, inode, ipoin, numer, pelty, pgaus, pnode, &
        !$OMP            volum, xjacm, xfact                                   ) &
        !$OMP  SHARED  ( coord, elmar, lnods, ltype, kfl_naxis, kfl_spher,       &
        !$OMP            nelem, ngaus, nnode,                                    &
#ifndef NDIMEPAR
        !$OMP            ndime,                                                  &
#endif
        !$OMP            vmass_chm                                               )
        !
        do ielem = 1,nelem
           pelty = ltype(ielem) 
           pnode = nnode(pelty)
           if( ndime == 1 ) then
              do inode=1,pnode
                 ipoin          = lnods(inode,ielem)
                 elcod(1,inode) = coord(1,ipoin) 
              end do
           else if( ndime == 2 ) then
              do inode=1,pnode
                 ipoin          = lnods(inode,ielem)
                 elcod(1,inode) = coord(1,ipoin)
                 elcod(2,inode) = coord(2,ipoin)
              end do
           else
              do inode = 1,pnode
                 ipoin          = lnods(inode,ielem)
                 elcod(1,inode) = coord(1,ipoin)
                 elcod(2,inode) = coord(2,ipoin)
                 elcod(3,inode) = coord(3,ipoin)
              end do
           end if

           if( pelty == PYR05 ) then
              !
              ! Pyramid element
              !
              pgaus = ngaus(pelty)
              volum = 0.0_rp
              denom = 0.0_rp
              do inode=1,pnode
                 numer(inode)=0.0_rp
                 do igaus=1,pgaus
                    call jacdet(&
                         ndime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                         xjacm,gpdet)
                    volum = volum + elmar(pelty)%weigp(igaus) * gpdet
                    denom = denom + elmar(pelty)%weigp(igaus) * gpdet &
                         *(elmar(pelty)%shape(igaus,inode))**2
                    numer(inode)=numer(inode)+elmar(pelty)%weigp(igaus)*gpdet&
                         *(elmar(pelty)%shape(igaus,inode))**2
                 end do
              end do
              xfact=volum/denom
              do inode=1,pnode
                 ipoin        = lnods(inode,ielem)
                 aux          = numer(inode)*xfact
                 !$OMP         ATOMIC
                 vmass_chm(ipoin) = vmass_chm(ipoin)+aux
              end do

           else
              !
              ! Other elements: loop over Gauss points (which are nodes)
              !
              if( kfl_naxis == 1 ) then

                 do inode=1,pnode
                    ipoin=lnods(inode,ielem)
                    call jacdet(&
                         ndime,pnode,elcod,elmar(pelty)%deric(1,1,inode),&
                         xjacm,gpdet)
                    gpvol=elmar(pelty)%weigc(inode)*gpdet
                    if( elcod(1,inode) == 0.0_rp ) then
                       gpvol = gpvol * twopi * 1.0e-12
                    else
                       gpvol = gpvol * twopi * elcod(1,inode)
                    end if
                    !$OMP         ATOMIC             
                    vmass_chm(ipoin) = vmass_chm(ipoin) + gpvol
                 end do

              else if( kfl_spher == 1 ) then

                 do inode=1,pnode
                    ipoin=lnods(inode,ielem)
                    call jacdet(&
                         ndime,pnode,elcod,elmar(pelty)%deric(1,1,inode),&
                         xjacm,gpdet)
                    gpvol=elmar(pelty)%weigc(inode)*gpdet
                    if( elcod(1,inode) == 0.0_rp ) then
                       gpvol = gpvol * twopi * 1.0e-12
                    else
                       gpvol = gpvol * 2.0_rp * twopi * elcod(1,inode) * elcod(1,inode)
                    end if
                    !$OMP         ATOMIC             
                    vmass_chm(ipoin) = vmass_chm(ipoin) + gpvol
                 end do

              else

                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    call jacdet(&
                         ndime,pnode,elcod,elmar(pelty)%deric(1,1,inode),&
                         xjacm,gpdet)
                    gpvol = elmar(pelty)%weigc(inode)*gpdet
                    !$OMP         ATOMIC             
                    vmass_chm(ipoin) = vmass_chm(ipoin) + gpvol
                 end do

              end if

           end if

        end do
        !$OMP END PARALLEL DO
        
     end if

  end if

end subroutine chm_splarr
