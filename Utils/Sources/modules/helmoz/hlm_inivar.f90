subroutine hlm_inivar(itask)

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_inivar.f90
  ! NAME 
  !    hlm_inivar
  ! DESCRIPTION
  !    This routine initializes some variables.
  ! USES
  ! USED BY
  !    hlm_turnon
  !-----------------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain
  use def_helmoz
  use def_solver
  use def_kermod
  use mod_iofile

  implicit none

  integer(ip), intent(in) :: itask

  integer(ip)             :: ielem,iedge,idime
  integer(ip)             :: ii,jj,kk,nn,ipoin,jpoin,iedgg
  integer(ip)             :: clnod1(nmlsi_hlm)
  real(rp)                :: clnod2(nmlsi_hlm) 
  integer(ip)             :: clsite1(1)
  real(rp)                :: clsite2(1) 
  real(rp)                :: tpoin(3)
  real(rp)                :: cpu_closnod1,cpu_closnod2

  select case (itask)

  case ( 0_ip )
     !
     ! Postprocess variables
     !
     postp(1) % wopos (1,1) = 'VECPO'
     postp(1) % wopos (1,2) = 'SCAPO'
     postp(1) % wopos (1,3) = 'ELEFI'
     postp(1) % wopos (1,4) = 'MAGFI'
     postp(1) % wopos (1,5) = 'DIFFJ'
     postp(1) % wopos (1,6) = 'DESIG'
     postp(1) % wopos (1,7) = 'CORRE'
     postp(1) % wopos (1,8) = 'DESCD'
     postp(1) % wopos (1,9) = 'FIXIT'

     postp(1) % wopos (2,1) = 'VECTX'
     postp(1) % wopos (2,2) = 'SCALX'
     postp(1) % wopos (2,3) = 'VECTX'
     postp(1) % wopos (2,4) = 'VECTX'
     postp(1) % wopos (2,5) = 'SCALA'
     postp(1) % wopos (2,6) = 'SCALA'
     postp(1) % wopos (2,7) = 'SCALA'
     postp(1) % wopos (2,8) = 'SCALA'
     postp(1) % wopos (2,9) = 'SCALA'

  case ( 1_ip ) 
     !
     ! Solver
     !
     call soldef(-1_ip)
     solve(1) % wprob     = 'MAXWELL'       ! Equation name: the gauged coupled vector-scalar potential formulation of Maxwell's equations
     solve(1) % kfl_solve = 1               ! Output flag
     !solve(1) % ndofn     = nequs_hlm       ! Number of degrees of freedom Vladimir
     solve(1) % ndofn     = 1               ! Number of degrees of freedom

!!$     if( INOTSLAVE ) then ! Vladimir
!!$        !write(*,*) 'Strategy: ', ielse(8)
!!$        call cputim(cpu_closnod1)
!!$        nn = nmlsi_hlm                     !Number of closest mesh nodes to a site needed for the MLSI
!!$        !write(*,*) 'nmlsi_hlm: ', nn
!!$        do ii = 1,nsite_hlm
!!$           !write(*,*) 'closnod: ', ii
!!$           !Test point
!!$           tpoin(1) = site_hlm(1,ii)
!!$           tpoin(2) = site_hlm(2,ii)
!!$           tpoin(3) = site_hlm(3,ii)
!!$           !write(*,*) 'tpoin(1): ', tpoin(1)
!!$           !write(*,*) 'tpoin(2): ', tpoin(2)
!!$           !write(*,*) 'tpoin(3): ', tpoin(3)
!!$           !Find N closest nodes of the mesh to a test point
!!$           call hlm_closnod(tpoin,nn,clnod1,clnod2)
!!$           kk = nn * (ii - 1_ip)
!!$           do jj = 1,nn
!!$              clnod1_hlm(kk+jj) = clnod1(jj)
!!$              clnod2_hlm(kk+jj) = clnod2(jj)
!!$              clcoor_hlm(1,kk+jj) = coord(1,clnod1(jj))
!!$              clcoor_hlm(2,kk+jj) = coord(2,clnod1(jj))
!!$              clcoor_hlm(3,kk+jj) = coord(3,clnod1(jj))
!!$           enddo
!!$        enddo
!!$        call cputim(cpu_closnod2)
!!$        !write(*,*)'Time to find the closest nodes: ', cpu_closnod2-cpu_closnod1
!!$        nn = 1_ip                         !Number of closest mesh node to a site
!!$        do ii = 1,nsite_hlm
!!$           !Test point
!!$           tpoin(1) = site_hlm(1,ii)
!!$           tpoin(2) = site_hlm(2,ii)
!!$           tpoin(3) = site_hlm(3,ii)
!!$           !Find 1 closest nodes of the mesh to a site
!!$           call hlm_closnod(tpoin,nn,clsite1,clsite2)
!!$           kk = nn * (ii - 1_ip)
!!$           do jj = 1,nn
!!$              clsite1_hlm(kk+jj) = clsite1(jj)
!!$              clsite_hlm(1,kk+jj) = coord(1,clsite1(jj))
!!$              clsite_hlm(2,kk+jj) = coord(2,clsite1(jj))
!!$              clsite_hlm(3,kk+jj) = coord(3,clsite1(jj))
!!$              clsite2_hlm(kk+jj) = clsite2(jj)
!!$           end do
!!$        end do
!!$     end if
     
  case ( 2_ip )

     if (INOTMASTER) then
        solve_sol(1)%limpo => kfl_fixno_hlm(1,1:)       
        call cregro()             
     endif

  case ( 3_ip )
     !
     ! Sign of the edges
     !
     if( kfl_edges_hlm == 1 ) then
        do ielem = 1,nelem
           do iedge = 1,meshe(ndivi) % lnned(ielem)
              iedgg = meshe(ndivi) % ledgs(iedge,ielem)
              ipoin = lninv_loc(meshe(ndivi) % edge_to_node(1,iedgg))
              jpoin = lninv_loc(meshe(ndivi) % edge_to_node(2,iedgg)) 
              if( ipoin > jpoin ) then
                 sign_edges_hlm(iedge,ielem) =  1.0_rp 
              else
                 sign_edges_hlm(iedge,ielem) = -1.0_rp 
              end if
           end do
        end do
        do iedgg = 1,meshe(ndivi) % nedge
           ipoin = meshe(ndivi) % edge_to_node(1,iedgg)
           jpoin = meshe(ndivi) % edge_to_node(2,iedgg)
           length_edges_hlm(iedgg) = 0.0_rp
           do idime = 1,ndime
              length_edges_hlm(iedgg) = length_edges_hlm(iedgg) + ( coord(idime,ipoin)-coord(idime,jpoin) )**2
           end do
           length_edges_hlm(iedgg) = sqrt( length_edges_hlm(iedgg) )
        end do
     end if

  end select

end subroutine hlm_inivar
