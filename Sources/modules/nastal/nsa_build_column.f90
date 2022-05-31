subroutine nsa_build_column
  !-----------------------------------------------------------------------
  !****f* nastal/nsa_build_column
  ! NAME 
  !    nsa_build_columb
  ! DESCRIPTION
  !    This routine builds the connectivity matrix/array of the columns
  !    needed by the Kessler microphysics.
  !
  !    Simone Marras (SM), Dec 2011, After Jim Kelly, NPS, March 2011
  !
  !
  ! NOTE: FOR NOW THIS ROUTINE WILL ONLY WORK IN SERIAL!!!
  !-----------------------------------------------------------------------

  use      def_master
  use      def_domain
  use      def_nastal

  implicit none

  integer(ip)   :: icol, iz, ipoin, iicol, ii, ie, k, ngl, nz, ncol
  real(rp)      :: x, y, eps, x1, y1
  logical isnew
 
  ngl = 1_ip
  eps = 1e-5_rp
  nz  = nelz_nsa + 1_ip                     !Number of points in the vertical
  ncol = ncol_nsa                           !Number of columns
  
  !Build nnode_column_nsa
  ipoin = 0
  do iz = 1,nz
     do icol = 1,ncol
        ipoin = ipoin + 1
        node_column_nsa(icol,iz) = ipoin
     end do
  end do
  
!!$  !Store NODE_COLUMN Data in Element-wise Fashion
!!$  do icol=1,ncol
!!$     ii=0
!!$     do ie=1,nelz_nsa
!!$        ii = ii+1
!!$        !if (ii == 0) &
!!$        !     ii = 1
!!$        
!!$        ipoin = node_column_nsa(icol,ii)
!!$        intma_column_nsa(icol,ie) = ipoin
!!$        ii=ii-2
!!$     end do !ie
!!$  end do !icol
  
end subroutine nsa_build_column
