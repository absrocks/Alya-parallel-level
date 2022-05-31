!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_bouope.f90
!> @author  Guillaume Houzeaux
!> @brief   Matrix assembly: boundary contribution
!> @details Boundary operations
!!          - Element matrix calculation
!!          - Scatter in global matrix
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_bouope()
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_nastal
  use mod_ker_proper 
  implicit none

  real(rp)    :: elrhs(nevat_nsa)

  real(rp)    :: baloc(ndime,ndime)                       ! Gather  
  real(rp)    :: bovel(ndime,mnodb)     
  real(rp)    :: bovfi(ndime,mnodb)     
  real(rp)    :: bocod(ndime,mnodb)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: elvel(ndime,mnode)

  real(rp)    :: gbcar(ndime,mnode,mgaus)
  real(rp)    :: gpcar(ndime,mnode,mgaus)
  real(rp)    :: gbsur(mgaub),eucta                 ! Values at Gauss points
  real(rp)    :: gbden(mgaub)
  real(rp)    :: gbvis(mgaub)
  real(rp)    :: tract(3),chale(3),chave(3),dummr
  real(rp)    :: gpvis,ustar,tragl(9),hleng(3)

  integer(ip) :: ielem,ipoin,igaus,inode,idime             ! Indices and dimensions
  integer(ip) :: pnode,pgaus,iboun,igaub,inodb
  integer(ip) :: pelty,pmate,pblty,pnodb,pgaub
  integer(ip) :: ievat,jevat,porde,kboun
  integer(ip) :: idofn
  !
  ! Loop over boundaries
  !
  boundaries: do iboun = 1,nboun

     if(  kfl_fixbo_nsa(iboun) ==  3 ) then    ! u.n in weak form
        !
        ! Element properties and dimensions
        !
        pblty = ltypb(iboun) 
        pnodb = nnode(pblty)
        ielem = lelbo(iboun)
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaub = ngaus(pblty) 
        pgaus = ngaus(pelty)
        porde = lorde(pelty)
        pmate = 1
        if( nmate > 1 ) pmate = lmate(ielem)

        if( pmate /= -1 ) then
           !
           ! Initialize
           !
           do ievat = 1,nevat_nsa
              elrhs(ievat) = 0.0_rp
           end do
           !
           ! Gather operations: ELVEL, ELCOD, BOVEL
           !
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 elvel(idime,inode) = veloc(idime,ipoin,1)
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do

           do inodb = 1,pnodb
              ipoin = lnodb(inodb,iboun)
              do idime = 1,ndime
                 bovel(idime,inodb) = veloc(idime,ipoin,1)
              end do
              if( kfl_coupl(ID_NASTAL,ID_ALEFOR) /= 0 ) then 
                 do idime = 1,ndime
                    bovfi(idime,inodb) = velom(idime,ipoin)
                 end do
              else
                 do idime = 1,ndime
                    bovfi(idime,inodb) = 0.0_rp
                 end do
              end if
           end do

           gauss_points: do igaub = 1,pgaub

              tract(1) = 0.0_rp
              tract(2) = 0.0_rp
              tract(3) = 0.0_rp

              bocod = 0.0_rp  !! DMM To be computed in nsa_elmgap

              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&    ! Cartesian derivative
                   bocod,baloc,eucta)                                   ! and Jacobian
              gbsur(igaub) = elmar(pblty) % weigp(igaub) * eucta 
              call chenor(pnode,baloc,bocod,elcod)                      ! Check normal

             if( kfl_fixbo_nsa(iboun) == 3 ) then
                 !
                 ! Wall law: sig.n = - rho* (U*^2) * (u_tan-u_fix_tan)/|u_tan-u_fix_tan|
                 !
                 gbvis = 0.0_rp
                 gbden = 0.0_rp
                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    gbvis(igaub) = gbvis(igaub) + visco(ipoin,1) * elmar(pblty) % shape(inodb,igaub)
                    gbden(igaub) = gbden(igaub) + densi(ipoin,1) * elmar(pblty) % shape(inodb,igaub)
                 end do
                 call nsa_bouwal(&                        
                      pnodb,iboun,lboel(1,iboun),elmar(pblty) % shape(1,igaub),bovel,bovfi,tract,&
                      gbvis(igaub),gbden(igaub),baloc,ustar,rough_dom)
                 ! Computation of RHS

                 do inodb = 1,pnodb                    
                    idofn = (lboel(inodb,iboun) - 1 ) * ndofn_nsa 
                    do idime = 1,ndime
                       idofn = idofn + 1
                       elrhs(idofn) = elrhs(idofn) + tract(idime) * elmar(pblty) % shape(inodb,igaub) * gbsur(igaub)
                    end do
                 end do

              end if

           end do gauss_points
           call assrhs(ndofn_nsa,pnode,lnods(1,ielem),elrhs,rhsid)

        end if

     end if

  end do boundaries

end subroutine nsa_bouope
