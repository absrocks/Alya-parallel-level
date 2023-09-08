subroutine got_elmgat(&
     order,itask,pnode,lnods,elvdr,elcdr,elcod,elvel,eldif)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmgat
  ! NAME 
  !    got_elmgat
  ! DESCRIPTION
  !    Gather operations
  ! INPUT
  !    VDROP ... Global droplet droplet velocity
  !    CDROP ... Global droplet concentration
  !    COORD ... Global coordinates
  !    VELOC ... Global Air droplet velocity
  ! OUTPUT
  !    ELVDR ... Element droplet velocity
  !    ELCDR ... Element droplet concentration
  !    ELCOD ... Element coordinates
  !    ELVEL ... Element air droplet velocity
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,npoin,coord
  use def_gotita, only     :  kfl_timei_got,diffm_got,kfl_diffu_got,&
       &                      kfl_probl_got
  use def_master, only     :  cdrop,vdrop,unkno
  implicit none
  integer(ip), intent(in)  :: order,itask,pnode
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: elvdr(ndime,pnode,2),elcdr(pnode,2)
  real(rp),    intent(out) :: elcod(ndime,pnode),elvel(ndime,pnode)
  real(rp),    intent(out) :: eldif(pnode)
  integer(ip)              :: inode,ipoin,idime

  select case(order)

  case(1)
     !
     ! Current droplet velocity, droplet concentration and coordinates
     !
     if(itask==1) then
        !
        ! Monolithic
        !
        do inode=1,pnode
           ipoin=lnods(inode)
           do idime=1,ndime
              elcod(idime,inode)   = coord(idime,ipoin)
           end do
        end do
        if(kfl_probl_got<=2) then
           do inode=1,pnode
              ipoin=lnods(inode)
              do idime=1,ndime
                 elvdr(idime,inode,1) = vdrop(idime,ipoin,1)
              end do
           end do
        end if
        if(kfl_probl_got/=2) then
           do inode=1,pnode
              ipoin=lnods(inode)
              elcdr(inode,1) = cdrop(ipoin,1)
              do idime=1,ndime
                 elvdr(idime,inode,1) = vdrop(idime,ipoin,1)
              end do
           end do
        end if
     else
        !
        ! Block Gauss-Seidel: VDROP and CDROP are in UNKNO
        !
        do inode=1,pnode
           ipoin=lnods(inode)
           do idime=1,ndime
              elvdr(idime,inode,1) = unkno((ipoin-1)*ndime+idime)
              elcod(idime,inode)   = coord(idime,ipoin)
           end do
           elcdr(inode,1)          = unkno(npoin*ndime+ipoin)
        end do
     end if
     !
     ! Velocity
     !
     if(kfl_probl_got<=2) then
        do inode=1,pnode
           ipoin=lnods(inode)
           call got_veloci(ipoin,elcod(1,inode),elvel(1,inode))
        end do
     end if
     !
     ! Time integration
     !
     if(kfl_timei_got==1) then
        if(kfl_probl_got==1) then
           do inode=1,pnode
              ipoin=lnods(inode)
              do idime=1,ndime
                 elvdr(idime,inode,2) = vdrop(idime,ipoin,3)
              end do
              elcdr(inode,2) = cdrop(ipoin,3)
           end do
        else if(kfl_probl_got==2) then
           do inode=1,pnode
              ipoin=lnods(inode)
              do idime=1,ndime
                 elvdr(idime,inode,2) = vdrop(idime,ipoin,3)
              end do
           end do
        else if(kfl_probl_got==3) then
           do inode=1,pnode
              ipoin=lnods(inode)
              elcdr(inode,2) = cdrop(ipoin,3)
           end do
        end if
     end if
     !
     ! Smoothed diffusion
     !
     if(kfl_diffu_got==2.and.kfl_probl_got<=2) then
        do inode=1,pnode
           ipoin=lnods(inode)
           eldif(inode)=diffm_got(ipoin)
        end do
     end if

  case(2)
     !
     ! Gather for time step calculations
     !
     if(kfl_probl_got==1) then
        do inode=1,pnode
           ipoin=lnods(inode)
           do idime=1,ndime
              elvdr(idime,inode,1) = vdrop(idime,ipoin,1)
              elcod(idime,inode)   = coord(idime,ipoin)
           end do
           elcdr(inode,1) = cdrop(ipoin,1)
           call got_veloci(ipoin,elcod(1,inode),elvel(1,inode))
        end do
     else if(kfl_probl_got==2) then
        do inode=1,pnode
           ipoin=lnods(inode)
           do idime=1,ndime
              elvdr(idime,inode,1) = vdrop(idime,ipoin,1)
              elcod(idime,inode)   = coord(idime,ipoin)
           end do
           call got_veloci(ipoin,elcod(1,inode),elvel(1,inode))
        end do
     else if(kfl_probl_got==3) then
        do inode=1,pnode
           ipoin=lnods(inode)
           do idime=1,ndime
              elvdr(idime,inode,1) = vdrop(idime,ipoin,1)
              elcod(idime,inode)   = coord(idime,ipoin)
           end do
           elcdr(inode,1) = cdrop(ipoin,1)
        end do
     end if

     if(kfl_probl_got<=2) then
        if(kfl_diffu_got==2) then
           do inode=1,pnode
              ipoin=lnods(inode)
              eldif(inode)=diffm_got(ipoin)
           end do
        end if
     end if

  case(3)
     !
     ! Gather for diffusion calculations
     !
     do inode=1,pnode
        ipoin=lnods(inode)
        do idime=1,ndime
           elvdr(idime,inode,1) = vdrop(idime,ipoin,1)
           elcod(idime,inode)   = coord(idime,ipoin)
        end do
        elcdr(inode,1) = cdrop(ipoin,1)
     end do

  end select

end subroutine got_elmgat
