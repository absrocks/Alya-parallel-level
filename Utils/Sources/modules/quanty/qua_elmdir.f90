subroutine qua_elmdir(itask,pnode,lnods,elmat,elrhs)
  !------------------------------------------------------------------------
  !****f* Quanty/qua_elmdir
  ! NAME 
  !    qua_elmdir
  ! DESCRIPTION
  !    This routine prescribes the boundary conditions for the equation. 
  ! USES
  ! USED BY
  !    qua_elmope
  !    qua_bouope
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_quanty, only       :  bvess_qua,bvessH_qua,kfl_fixno_qua
  use def_quanty, only       :  ncomp_eig,kfl_dftgs_qua,kfl_alele_qua 
  use def_master, only       :  eigen
  implicit none
  integer(ip), intent(in)    :: pnode,itask
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)                :: inode,ipoin,jnode,ii
  real(rp)                   :: adiag

  if( itask == 0 ) then  
     !
     ! Caso KhonSham
     !
     do inode = 1,pnode
        ipoin = lnods(inode)
        if(  kfl_fixno_qua(1,ipoin) == 1 .or. &
             kfl_fixno_qua(1,ipoin) == 4 .or. &
             kfl_fixno_qua(1,ipoin) == 5 ) then
           adiag  = elmat(inode,inode)
           do jnode = 1,pnode
              elmat(inode,jnode) = 0.0_rp
              elmat(jnode,inode) = 0.0_rp
           end do
           elmat(inode,inode) = adiag
        end if
     end do

  else  
     ! 
     ! Caso Poisson
     !
     do inode = 1,pnode
        ipoin = lnods(inode)
        if(  kfl_fixno_qua(1,ipoin) == 1 .or. &
             kfl_fixno_qua(1,ipoin) == 4 .or. &
             kfl_fixno_qua(1,ipoin) == 5 ) then
           adiag  = elmat(inode,inode)
           do jnode = 1,pnode
              elmat(inode,jnode) = 0.0_rp
              elrhs(jnode)       = elrhs(jnode)-elmat(jnode,inode)*bvessH_qua(ipoin,1)
              elmat(jnode,inode) = 0.0_rp
           end do
           elmat(inode,inode) = adiag
           elrhs(inode)       = adiag*bvessH_qua(ipoin,1)
        end if
     end do

  endif


  !  agrego las condiciones de contorno a los eigenvectors OJO guillaume ?

  !  if(kfl_dftgs_qua /=0 .or. kfl_alele_qua/=0) then

  !    do ii=1,ncomp_eig
  !      do inode = 1,pnode
  !        ipoin = lnods(inode)

  !	    if(  kfl_fixno_qua(1,ipoin) == 1 .or. &
  !          kfl_fixno_qua(1,ipoin) == 4 .or. &
  !          kfl_fixno_qua(1,ipoin) == 5 ) then

  !			  eigen((ii-1)*pnode + ipoin) = 0.0_rp

  !	    end if

  !      end do
  !    enddo

  !  endif

end subroutine qua_elmdir
