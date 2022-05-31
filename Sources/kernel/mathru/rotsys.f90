!------------------------------------------------------------------------
!> @addtogroup Mathru
!> @{
!> @file    rotsys.f90
!> @date    10/12/2013
!> @author  Mariano Vazquez
!> @brief   Rotate elmat and elrhs
!> @details Rotate elmat and elrhs
!> @}
!------------------------------------------------------------------------
subroutine rotsys(&
     itask,imodi,pnode,pdofn,pevat,elmat,elrhs,rotma)
  !------------------------------------------------------------------------
  !
  ! This routine modifies a matrix the components of which are
  ! submatrices by rotating them using a matrix R ( = ROTMA) as follows:
  !
  !                |       A_11  ....      A_1i  R   ....      A_1n   |
  !                |   ............................................   |
  !     ELMAT <--  |   R^t A_i1  ....  R^t A_ii  R   ....  R^t A_in   |
  !                |   ............................................   |
  !                |       A_n1  ....      A_ni  R   ....      A_nn   |
  !
  ! where i ( = IMODI) is a given position. Also, the i-subvector of the
  ! given vector BVECT is multiplied by R^t. 
  !
  !------------------------------------------------------------------------
  use def_parame
  implicit none
  integer(ip), intent(in)    :: pnode,pevat,pdofn,itask
  real(rp),    intent(in)    :: rotma(pdofn,pdofn)
  real(rp),    intent(out)   :: elmat(pevat,pevat)
  real(rp),    intent(out)   :: elrhs(pevat)
  real(rp)                   :: worma(pevat,pevat)
  integer(ip)                :: imodi,itot0,jtot0,inode,jnode,itotp,jtotp
  integer(ip)                :: idofn,jdofn,kdime,itott,jtott,ktott,pauxi
  !
  ! Modifies column number IMODI of ELMAT ( A_j,imodi <-- A_j,imodi R )
  !
  
  if (itask .le. 1) then
     jtot0=(imodi-1)*pdofn
     do inode=1,pnode
        itot0=(inode-1)*pdofn
        do idofn=1,pdofn
           itott=itot0+idofn
           do jdofn=1,pdofn
              worma(idofn,jdofn)=0.0_rp
              do kdime=1,pdofn
                 ktott=jtot0+kdime
                 worma(idofn,jdofn)=worma(idofn,jdofn)&
                      +elmat(itott,ktott)*rotma(kdime,jdofn)
              end do
           end do
        end do
        do idofn=1,pdofn
           itott=itot0+idofn
           do jdofn=1,pdofn
              jtott=jtot0+jdofn
              elmat(itott,jtott)=worma(idofn,jdofn)
           end do
        end do
     end do
     !
     ! Modifies row number IMODI of ELMAT ( A_imodi,j <-- R^t A_imodi,j )
     !
     itot0=(imodi-1)*pdofn
     do jnode=1,pnode
        jtot0=(jnode-1)*pdofn
        do idofn=1,pdofn
           do jdofn=1,pdofn
              jtott=jtot0+jdofn
              worma(idofn,jdofn)=0.0_rp
              do kdime=1,pdofn
                 ktott=itot0+kdime
                 worma(idofn,jdofn)=worma(idofn,jdofn)&
                      +rotma(kdime,idofn)*elmat(ktott,jtott)
              end do
           end do
        end do
        do idofn=1,pdofn
           itott=itot0+idofn
           do jdofn=1,pdofn
              jtott=jtot0+jdofn
              elmat(itott,jtott)=worma(idofn,jdofn)
           end do
        end do
     end do
  end if
  !
  ! Modifies the IMODI subvector of the vector ELRHS
  !
  if (itask .le. 2) then
     itot0=(imodi-1)*pdofn
     do idofn=1,pdofn
        worma(idofn,1)=0.0_rp
        do kdime=1,pdofn
           ktott=itot0+kdime
           worma(idofn,1)=worma(idofn,1)&
                +rotma(kdime,idofn)*elrhs(ktott)
        end do
     end do
     do idofn=1,pdofn
        itott=itot0+idofn
        elrhs(itott)=worma(idofn,1)
     end do
  end if
  
end subroutine rotsys
