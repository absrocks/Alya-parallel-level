!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_copy_mesh.f90
!> @author  Guillaume Houzeaux
!> @date    27/02/2013
!> @brief   Copy mesh in subdomain-wise arrays
!> @details Copy mesh in subdomain-wise arrays
!> @} 
!-----------------------------------------------------------------------
subroutine dod_copy_mesh()

  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme 
  use mod_memory
  use mod_elmgeo
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: isubd,ipoin,ielem,inode,iboun,idime
  integer(ip) :: pnode,pnodb,inodb,izone,kelem,ninve
  integer(ip) :: ipoin_subd,ielem_subd,iboun_subd,ii
  integer(ip) :: lnodb_aux(mnodb),lboel_aux(mnodb)
  real(rp)    :: dummr(3)

  call livinf(0_ip,'COPY MESH ARRAYS TO SUBDOMAIN ARRAYS',0_ip)
  !
  ! Mesh dimensions
  !
  do isubd = 1,nsubd
     subdomain(isubd) % npoin = 0
     subdomain(isubd) % nelem = 0
     subdomain(isubd) % nboun = 0
  end do
  do ielem = 1,nelem
     isubd                    = lsubd_nelem(ielem)
     subdomain(isubd) % nelem = subdomain(isubd) % nelem + 1
     linvp_nelem(ielem)       = subdomain(isubd) % nelem
     pnode                    = nnode(abs(ltype(ielem)))
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        if( lsubd_npoin(ipoin) == 0 ) then
           subdomain(isubd) % npoin = subdomain(isubd) % npoin + 1
           linvp_npoin(ipoin)       = subdomain(isubd) % npoin
        end if
        lsubd_npoin(ipoin) = isubd
     end do
  end do
  do iboun = 1,nboun
     pnodb                    = nnode(abs(ltypb(iboun)))
     ielem                    = lelbo(iboun)
     isubd                    = lsubd_nelem(ielem)
     lsubd_nboun(iboun)       = isubd
     subdomain(isubd) % nboun = subdomain(isubd) % nboun + 1
     linvp_nboun(iboun)       = subdomain(isubd) % nboun
  end do
  !
  ! Allocate memory for subdomain arrays
  !
  call dod_memall(2_ip)
  !
  ! Copy nodal mesh arrays: LNPER, COORD
  !
  do ipoin = 1,npoin
     isubd                                = lsubd_npoin(ipoin)
     ipoin_subd                           = linvp_npoin(ipoin)
     subdomain(isubd) % lnper(ipoin_subd) = ipoin
     do idime = 1,ndime
        subdomain(isubd) % coord(idime,ipoin_subd) = coord(idime,ipoin)
     end do
  end do
  !
  ! Copy element mesh arrays: LNODS, LTYPE, LELCH, LNNOD
  !
  do ielem = 1,nelem
     isubd                                = lsubd_nelem(ielem)
     ielem_subd                           = linvp_nelem(ielem)
     pnode                                = nnode(abs(ltype(ielem)))     
     subdomain(isubd) % leper(ielem_subd) = ielem
     do inode = 1,pnode
        ipoin      = lnods(inode,ielem)
        ipoin_subd = linvp_npoin(ipoin)
        subdomain(isubd) % lnods(inode,ielem_subd) = ipoin_subd        
     end do
     subdomain(isubd) % ltype(ielem_subd) = ltype(ielem)
     subdomain(isubd) % lelch(ielem_subd) = lelch(ielem)
     subdomain(isubd) % lnnod(ielem_subd) = lnnod(ielem)
  end do
  !
  ! Materials
  !
  do ielem = 1,nelem
     pnode = nnode(abs(ltype(ielem)))   
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        lmatn_dod(ipoin) = lmate(ielem) 
     end do
  end do
  !
  ! Zones
  !
  do izone = 1,nzone
     do ielem = 1,nelem
        pnode = nnode(abs(ltype(ielem)))   
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           lpoiz_dod(ipoin) = izone
        end do
     end do
  end do
  !
  ! Copy boundary mesh arrays: LNODB, LBOCH, LTYPB, LBOEL
  !
  do iboun = 1,nboun
     isubd                                = lsubd_nboun(iboun)
     pnodb                                = nnode(abs(ltypb(iboun)))
     ielem                                = lelbo(iboun)
     iboun_subd                           = linvp_nboun(iboun)
     ielem_subd                           = linvp_nelem(ielem)
     subdomain(isubd) % lbper(iboun_subd) = iboun
     do inodb = 1,pnodb
        ipoin                                      = lnodb(inodb,iboun)
        ipoin_subd                                 = linvp_npoin(ipoin)
        subdomain(isubd) % lnodb(inodb,iboun_subd) = ipoin_subd
        subdomain(isubd) % lboel(inodb,iboun_subd) = lboel(inodb,iboun)
        subdomain(isubd) % lelbo(iboun_subd)       = lelbo(inodb,iboun)
     end do
     subdomain(isubd) % ltypb(iboun_subd) = ltypb(iboun)
     subdomain(isubd) % lboch(iboun_subd) = lboch(iboun)
     subdomain(isubd) % lelbo(iboun_subd) = ielem_subd
     ninve = 0_ip
     call  elmgeo_bounor(&
       1_ip,1_ip,ndime,pnodb,mnode,subdomain(isubd) % lnodb(1:pnodb,iboun_subd),subdomain(isubd) %ltypb(iboun_subd:),&
       subdomain(isubd) % lelbo(iboun_subd),subdomain(isubd) % ltype,subdomain(isubd) % lnods,nnode(1:),subdomain(isubd) % coord,ninve,dummr)
     if(ninve/=0_ip)then
        do inodb=1,pnodb
           lnodb_aux(inodb) = subdomain(isubd)%lnodb(inodb,iboun_subd)
           lboel_aux(inodb) = subdomain(isubd)%lboel(inodb,iboun_subd)
        end do
        ii = 0
        do inodb=2,pnodb
           subdomain(isubd)%lnodb(inodb,iboun_subd) = lnodb_aux(pnodb-ii)
           subdomain(isubd)%lboel(inodb,iboun_subd) = lboel_aux(pnodb-ii)
           ii = ii + 1
       end do
     end if
  end do
  !
  ! Deallocate inverse permutation arrays
  !
  call dod_memall(-5_ip)

end subroutine dod_copy_mesh
