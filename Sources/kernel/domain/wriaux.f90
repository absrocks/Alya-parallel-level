! must be aclled at the end of Turnon

!-----------------------------------------------------------------------
!> @addtogroup DomainInput
!> @{
!> @file    wriaux.f90
!> @author  Herbert Owen 
!> @date    27/04/2014
!> @brief   Write a modified mesh
!> @details In this first case a 2d mesh from a 3d one
!> @} 
!-----------------------------------------------------------------------
subroutine wriaux()
  use def_kintyp
  use def_parame
  use def_master 
  use def_domain
  use def_inpout
  use mod_iofile
  use mod_memory
  use def_elmtyp
  implicit none
  integer(ip)           :: ielem,ipoin,jpoin,inode,idime,iimbo
  integer(ip)           :: iboun,inodb,ielty,ktype,dummi,kfl_defau
  integer(ip)           :: iskew,jskew,jdime,imate,kelem,nlimi
  integer(ip)           :: iblty,knode,kfl_gidsk,kfl_dontr,kpoin
  integer(ip)           :: kfl_binme,ipara,iperi,imast,kfl_icgns
  integer(ip)           :: pelty,ifiel,izone,kfl_ifmat,kfl_elino
  integer(ip)           :: isubd
  real(rp)              :: dummr
  character(20)         :: chnod

  integer(ip),  pointer :: lnnhh(:)
  integer(ip)           :: led(2,4)
  integer(ip),parameter :: notherbo=3
  integer(ip)           :: lotherbo(notherbo)
  integer(ip)           :: i,j,kount,nnpoin,pblty,pnodb
  
  integer(ip)           :: codextract,nnelem,iob,nnboun

  if( INOTSLAVE ) then
     codextract = 4_ip  ! atocar a mano
     lotherbo(1)=1
     lotherbo(2)=2
     lotherbo(3)=3
     !
     ! write type and obtain lnnhh
     !
     write(777,*)'TYPES'
     nnelem = 0_ip
     nnpoin = 0_ip
     nullify(lnnhh)
     call memory_alloca(memor_dom,'LNNHH','wriaux',lnnhh,npoin)
     lnnhh = 0_ip
     do iboun = 1,nboun
        if ((kfl_codbo(iboun) == codextract) .and. coord(2,lnodb(1,iboun)) > -0.5_rp) then
           nnelem = nnelem + 1_ip
           write(777,*) nnelem,ltypb(iboun)
           pblty = abs(ltypb(iboun))
           pnodb = nnode(pblty)
           do inodb=1,pnodb
              if ( lnnhh(lnodb(inodb,iboun)) == 0 ) then
                 nnpoin = nnpoin + 1_ip
                 lnnhh(lnodb(inodb,iboun)) = nnpoin
              end if
           end do
        end if
     end do
     write(777,*)'END_TYPES'

     write(777,*)'ELEMENTS'
     nnelem = 0_ip
     do iboun = 1,nboun
        if ((kfl_codbo(iboun) == codextract) .and. coord(2,lnodb(1,iboun)) > -0.5_rp) then
           nnelem = nnelem + 1_ip
           pblty = abs(ltypb(iboun))
           pnodb = nnode(pblty)
!          write(777,*) nnelem,(lnnhh(lnodb(inodb,iboun)),inodb=1,pnodb)
           write(777,*) nnelem,lnnhh(lnodb(1,iboun)),lnnhh(lnodb(4,iboun)),lnnhh(lnodb(3,iboun)),lnnhh(lnodb(2,iboun))
        end if
     end do
     write(777,*)'END_ELEMENTS'
     write(*,*)'!!!!!!!!!!!!!!!END_ELEMENTS'
     write(*,*)'size(lnnhh),npoin',size(lnnhh,KIND=ip),npoin


     !ojo no van a salir ordnados, importa?????
     write(777,*)'COORDINATES   NOT_SORTED'
     do ipoin = 1,npoin
        if (lnnhh(ipoin) /= 0_ip ) write(777,*) lnnhh(ipoin),coord(1,ipoin),coord(3,ipoin)
     end do
     write(777,*)'END_COORDINATES'
     write(*,*)'!!!!!!!!!!!END_COORDINATES'
     flush(6)
     flush(777)


     led(1,1) = 1
     led(2,1) = 2
     led(1,2) = 2
     led(2,2) = 3
     led(1,3) = 3
     led(2,3) = 4
     led(1,4) = 4
     led(2,4) = 1



     write(777,*)'BOUNDARIES, ELEMENT'
     nnelem = 0_ip
     nnboun = 0_ip
     do iboun = 1,nboun
        if ((kfl_codbo(iboun) == codextract) .and. coord(2,lnodb(1,iboun)) > -0.5_rp) then
           nnelem = nnelem + 1_ip
           pblty = abs(ltypb(iboun))
           pnodb = nnode(pblty)
           do iob = 1,notherbo
              do inodb=1,pnodb   ! acually he we use pnodb=nedge for 2d
                 kount=0
                 do i=1,2
                    do j=1,3
                       if (kfl_codno(j,lnodb(led(i,inodb),iboun)) == lotherbo(iob)) kount=kount+1
                    end do
                 end do
                 if (kount==2) then
                    nnboun = nnboun+1
                    write(777,*) nnboun,lnnhh(lnodb(led(1,inodb),iboun)),lnnhh(lnodb(led(2,inodb),iboun)),nnelem
                    write(778,*) nnboun,lotherbo(iob)
                 end if
              end do
           end do
        end if
     end do
     write(777,*)'END_BOUNDARIES'

     call memory_deallo(memor_dom,'LNNHH','wriaux',lnnhh)


  end if   !notslave
  print*,'SALIOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
end subroutine wriaux
