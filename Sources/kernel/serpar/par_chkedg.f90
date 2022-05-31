!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_chkedg.f90
!> @date    24/10/2013
!> @author  Guillaume Houzeaux
!> @brief   Check if an edge is mine
!> @details Given two nodes, PARD1 and PARD2, check if the edge 
!>          PARD1-PARD2 is shared by a neighbor. If this is the
!>          case, the one who owns the edge if the subdomain with
!>          lowest rank. Then, PARD3 = 1, otherwise, PARD3 = 0
!> @} 
!-----------------------------------------------------------------------
subroutine par_chkedg()
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_parall, only :  commd
  implicit none
  integer(ip) :: which_dom,knode
  integer(ip) :: ineig,dom_i,ini,pardk

  pard3 = 1
  if( ISLAVE ) then
     !
     ! Loop over the neighbors to check if they have edge PARD1-PARD2
     !
     which_dom = kfl_paral
     do ineig = 1,commd % nneig
        knode = 0
        dom_i = commd % neights(ineig)
        ini   = commd % bound_size(ineig)
        do while( ini < commd % bound_size(ineig+1)-1 )
           pardk = commd % bound_perm(ini)
           if( pardk == pard1 .or. pardk == pard2 ) knode = knode + 1
           ini   = ini + 1
        end do
        if( knode == 2 ) then
           !
           ! My neigbor DOM_I owns this edge too, decide who owns it
           !
           which_dom = min(which_dom,dom_i)           
        end if
     end do
     if( which_dom /= kfl_paral ) pard3 = 0
  end if

end subroutine par_chkedg
