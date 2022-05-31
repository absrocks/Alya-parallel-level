subroutine nsa_elmean(pnode,elvel,elpre,elunk,emene,emvel,emve2,emso2,emve3,emesv)
!-----------------------------------------------------------------------
! 
! Elemental mean values of:
! 
!    e                    1       emene
!    c^2                  1       emso2
!    k  e / rho           1       ...
!    k / rho              1
!    mu / rho             1
!    mu / rho^2           1
!    u_i c^2              ndime
!    mu u_i / rho         ndime
!    k  u_i / rho         ndime
!    u_i                  ndime   emvel
!    u_i e                ndime   
!    u_i u_j              ntens   emve2
!    u_i u_j u_j          ntens
!    mu u_i u_j / rho     ntens
!    k  u_i u_j / rho     ntens
!
!-----------------------------------------------------------------------------------------------------------
  use      def_kintyp
  use      def_domain
  use      def_nastal
  implicit none
  integer(ip)         :: inode,pnode,idime,jdime,itens,jtens
  real(rp)            :: &
       elvel(ndime,mnode),elpre(mnode),elunk(ndofn_nsa,mnode,ncomp_nsa),&
       elvis(mnode),emene,emvel,emve2,emso2,emve3,emesv
!  real(rp) ::  emean(6+5*ndime+4*ntens)
  


  emene        = 0.0_rp
  emso2        = 0.0_rp
  emvel        = 0.0_rp
  emve2        = 0.0_rp
  emve3        = 0.0_rp
  emesv        = 0.0_rp

  do inode=1,pnode
     emene        = emene + elunk(ndime+2,inode,1) / elunk(ndime+1,inode,1)
     emso2        = emso2 + adgam_nsa * elpre(inode) / elunk(ndime+1,inode,1)
     emvel        = emvel + elvel(1,inode)
     emve2        = emve2 + elvel(1,inode) * elvel(1,inode)
     emve3        = emve3 + elvel(1,inode) * elvel(1,inode) * elvel(1,inode)
     emesv        = emesv + elvel(1,inode) * (elunk(ndime+2,inode,1) + elpre(inode)) / elunk(ndime+1,inode,1)
  end do

  emene        = emene / real(pnode)
  emso2        = emso2 / real(pnode)
  emvel        = emvel / real(pnode)
  emve2        = emve2 / real(pnode)
  emve3        = emve3 / real(pnode)
  emesv        = emesv / real(pnode)



!  itens= 6 + 5 * ndime
!  emean(1) = elene(1) / elden(1)
!  emean(2) = adgam_nsa * elpre(1) / elden(1)
!  emean(3) = cppra_nsa * elvis(1) * elene(1) / elden(1) / elden(1)
!  emean(4) = cppra_nsa * elvis(1) / elden(1)
!  emean(5) = elvis(1) / elden(1)
!  emean(6) = elvis(1) / elden(1) / elden(1)
!  jtens= 0
!  do idime=1,ndime
!     emean(6+  idime)       = elvel(idime,1) * adgam_nsa * elpre(1) / elden(1)
!     emean(6+  ndime+idime) = elvel(idime,1) * elvis(1) / elden(1)
!     emean(6+2*ndime+idime) = cppra_nsa * elvel(idime,1) * elvis(1) / elden(1)
!     emean(6+3*ndime+idime) = elvel(idime,1)
!     emean(6+4*ndime+idime) = elvel(idime,1) * elene(1) / elden(1)     
!     do jdime=1,idime
!        jtens=jtens+1
!        emean(itens+        jtens) = elvel(idime,1) * elvel(jdime,1)
!        emean(itens+  ntens+jtens) = elvel(idime,1) * elvel(jdime,1) * elvel(jdime,1)
!        emean(itens+2*ntens+jtens) = elvis(1) * elvel(idime,1) * elvel(jdime,1) / elden(1) 
!        emean(itens+3*ntens+jtens) = cppra_nsa * elvis(1) * elvel(idime,1) * elvel(jdime,1) / elden(1) 
!     end do
!  end do
!
!  do inode= 2,pnode     
!     emean(1) = emean(1) +elene(inode) / elden(inode)
!     emean(2) = emean(2) +adgam_nsa * elpre(inode) / elden(inode)
!     emean(3) = emean(3) +cppra_nsa * elvis(inode) * elene(inode) / elden(inode) / elden(inode)
!     emean(4) = emean(4) +cppra_nsa * elvis(inode) / elden(inode)
!     emean(5) = emean(5) +elvis(inode) / elden(inode)
!     emean(6) = emean(6) +elvis(inode) / elden(inode) / elden(inode)
!     jtens= 0
!     do idime=1,ndime
!        emean(6+  idime)       = emean(6+  idime)       + elvel(idime,inode) * adgam_nsa * elpre(inode) / elden(inode)
!        emean(6+  ndime+idime) = emean(6+  ndime+idime) + elvel(idime,inode) * elvis(inode) / elden(inode)
!        emean(6+2*ndime+idime) = emean(6+2*ndime+idime) + cppra_nsa * elvel(idime,inode) * elvis(inode) / elden(inode)
!        emean(6+3*ndime+idime) = emean(6+3*ndime+idime) + elvel(idime,inode)
!        emean(6+4*ndime+idime) = emean(6+4*ndime+idime) + elvel(idime,inode) * elene(inode) / elden(inode)     
!        do jdime=1,idime
!           jtens=jtens+1
!           emean(itens+        jtens) = emean(itens+        jtens) + &
!                elvel(idime,inode) * elvel(jdime,inode)
!           emean(itens+  ntens+jtens) = emean(itens+  ntens+jtens) + &
!                elvel(idime,inode) * elvel(jdime,inode) * elvel(jdime,inode)
!           emean(itens+2*ntens+jtens) = emean(itens+2*ntens+jtens) + &
!                elvis(inode) * elvel(idime,inode) * elvel(jdime,inode) / elden(inode) 
!           emean(itens+3*ntens+jtens) = emean(itens+3*ntens+jtens) + &
!                cppra_nsa * elvis(inode) * elvel(idime,inode) * elvel(jdime,inode) / elden(inode) 
!        end do
!     end do
!     
!  end do

end subroutine nsa_elmean
