!------------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    chm_partis.f90
!> @author  Guillaume Houzeaux
!> @brief   Coupling with Partis module
!> @details Exchange mass source in energy equation 
!> @}
!------------------------------------------------------------------------

subroutine chm_partis()
  use def_master
  use def_domain
  use mod_parall,  only : commd
  use def_chemic,  only : nspec_chm
  implicit none
  integer(ip) :: ipoin,iclas,iclas_source,idofn,jj
  real(rp)    :: multiplicity_loc

  if( associated(mass_sink) ) then

     jj    = 0 
     iclas_source = 3
     do ipoin = 1,npoin
        if (ipoin <= npoi1) then
           multiplicity_loc = 1.0_rp
        else
           jj = jj + 1
           multiplicity_loc = real(commd % bound_multiplicity(jj),rp)
        endif
        do iclas = 1,nspec_chm
           idofn = (ipoin - 1) * nspec_chm + iclas

           !
           ! Add vapour source:
           !
           if (iclas == iclas_source) &
              rhsid(idofn)  = rhsid(idofn) + mass_sink(ipoin)/multiplicity_loc

           !
           ! Subtract dilution term:
           !
           rhsid(idofn)  = rhsid(idofn) - mass_sink(ipoin)/multiplicity_loc * unkno(idofn) 
        enddo
     end do
  end if

  !real(rp)    :: elmassk(mnode)
  !real(rp)    :: gpvol(mgaus)
  !real(rp)    :: gpmassk(mgaus)
  !real(rp)    :: gpcar(ndime,mnode,mgaus) 
  !real(rp)    :: gphes(ntens,mnode) 
  !real(rp)    :: elcod(ndime,mnode) 
  !integer(ip) :: pnode,pgaus,pelty
  !integer(ip) :: ielem,igaus,inode
  !!
  !! Loop over elements
  !!
  !iclas_source = 3
  !elements: do ielem = 1,nelem
  !   ! 
  !   ! Element properties and dimensions
  !   !
  !   pelty = ltype(ielem)
  !   if( pelty > 0 ) then
  !      pnode = nnode(pelty)
  !      pgaus = ngaus(pelty)
  !      !
  !      ! Gather operations
  !      !
  !      do inode = 1,pnode
  !         ipoin = lnods(inode,ielem)
  !         elcod(1:ndime,inode) = coord(1:ndime,ipoin)
  !      enddo
  !      elmassk(:) = 0.0_rp
  !      if( associated(mass_sink) ) then
  !         do inode = 1,pnode
  !            ipoin = lnods(inode,ielem)
  !            elmassk(inode) = mass_sink(ipoin)
  !         enddo
  !      endif
  !      call elmcar(&
  !           pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
  !           elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
  !           gphes,ielem)
  !         
  !      gpmassk(:) = 0.0_rp
  !      do igaus = 1,pgaus      
  !         do inode = 1,pnode
  !            gpmassk(igaus) = gpmassk(igaus)                   &
  !                          + elmar(pelty) % shape(inode,igaus) &
  !                          * elmassk(inode)              
  !         enddo
  !      enddo

  !      !
  !      ! Add to right hand side
  !      !
  !      do igaus = 1,pgaus     
  !         do inode = 1,pnode
  !            ipoin = lnods(inode,ielem)
  !            do iclas = 1,nclas_chm
  !               idofn = (ipoin - 1) * nclas_chm + iclas
  !               if( kfl_fixno_chm(iclas,ipoin) <= 0 ) then
  !                  !
  !                  ! Add source term:
  !                  !
  !                  if (iclas == iclas_source) &
  !                     rhsid(idofn) = rhsid(idofn) + elmar(pelty) % shape(inode,igaus) * gpvol(igaus) * gpmassk(igaus)  
  !                  !
  !                  ! Subtract dilution term:
  !                  !
  !                  rhsid(idofn) = rhsid(idofn) -  elmar(pelty) % shape(inode,igaus) * gpvol(igaus) * gpmassk(igaus) * unkno(idofn)
  !               end if
  !            enddo
  !         enddo
  !      end do
  !   endif



  !end do elements

end subroutine chm_partis
