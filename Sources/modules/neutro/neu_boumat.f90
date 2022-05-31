!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_boumat.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Element gather
!> @details Element gather for the boundary elements
!> @}
!-----------------------------------------------------------------------

subroutine neu_boumat(pnodb,pgaub,pnode,lboel,gbsha,gbsur,baloc,borad,bvnat,elmat,elrhs)

  use def_kintyp, only :  ip,rp
  use def_parame, only :  pi,Stefan_Boltzmann
  use def_domain, only :  ndime
  use def_neutro, only :  ADR_neu
  use def_neutro, only :  num_energies_neu
  use def_neutro, only :  num_directions_neu
  use def_neutro, only :  current_energy_neu
  use def_neutro, only :  current_direction_neu
  use def_neutro, only :  direc_neu
  use def_neutro, only :  weigd_neu
  use def_neutro, only :  scattering_neu
  use def_neutro, only :  nitsche_neu
  use def_neutro, only :  kfl_smobo_neu
  implicit none
  integer(ip), intent(in)  :: pnodb !> Number of nodes in the  boundary element
  integer(ip), intent(in)  :: pgaub !> Number of gauss points in the boundary element
  integer(ip), intent(in)  :: pnode !> Number of nodes in the element (Domain)
  integer(ip), intent(in)  :: lboel(pnodb)  ! connectivity matrix in the contour (relates number of nodes with elements)
  real(rp),    intent(in)  :: gbsha(pnodb,pgaub)  !> test/shape function, it is saved in elmar % shape
  real(rp),    intent(in)  :: gbsur(pgaub) !> surface , jacobian times weight of the Gauss point
  real(rp),    intent(in)  :: baloc(ndime,ndime,pgaub)   !> Local directions, normal to boundary element (baloc(1:ndime, ndime))  and two tangencial (baloc(1:ndime,1) and baloc(1:ndime,2) for 3D 
  real(rp),    intent(in)  :: borad(num_energies_neu,num_directions_neu,pnodb) !> Boundary radiation for each energy, direction and boundary node
  real(rp),    intent(in)  :: bvnat !> Natural boundary conditions
  real(rp),    intent(out) :: elmat(pnode,pnode) !> Element matrix
  real(rp),    intent(out) :: elrhs(pnode) !> Element right-hand side
  integer(ip)              :: idofn,jdofn,jdir,out !> Indexes
  integer(ip)              :: outgo(num_directions_neu) !> Outward direction
  integer(ip)              :: igaub,inodb,inode !> Indexes for boundary gauss point, boundary node and element node
  real(rp)                 :: btemp!,bolt4!> Used to be the radiation at the boundary (boundary temperature)
  real(rp)                 :: ndots(num_directions_neu),refle !> producto escalar de direccion con normal hacia fuera (n dot s)
  real(rp)                 :: tract,right(pnode)!> variables intermedias para computar condiciones de contorno
  real(rp)                 :: work1,dotpr,mmwat !> variables intermedias
  real(rp)                 :: gprad(num_energies_neu,num_directions_neu)!> Radiacion interpolada a los puntos de gauss.

  real(rp) :: emiss,bsour !>Boundary SOURce

  refle = 0.0_rp
  emiss = 1.0_rp
  bsour = 100000.0_rp
  igaub = 1
  !
  ! Exiting direction, no boundary condition
  !
  dotpr = dot_product(direc_neu(1:ndime,current_direction_neu),baloc(1:ndime,ndime,igaub))
  if( dotpr > 0.0_rp ) return  ! if outgoing direction then return

  gauss: do igaub = 1,pgaub  ! impose b.c.

     ndots(current_direction_neu) = dotpr
     !btemp        = btemp_neu(iboun, igaub)
     !bolt4        = btemp*btemp*btemp*btemp*boltz   ! For us this is bvnat
     !bolt4        = bvnat   ! For us this is bvnat
     !bolt4 = (tempe**4.0_rp)*Stefan_Boltzmann !>We no longer need to calculate the radiation
     !
     ! if nodal matrices, add contibution of nitsche's smooth term to matrix
     !
     do inodb = 1,pnodb
        inode              = lboel(inodb) !> We find to which node the boundary point corresponds to
        mmwat              = -ndots(current_direction_neu) * nitsche_neu * gbsha(inodb,igaub) !> Gauss matrix contribution
        elmat(inode,inode) = elmat(inode,inode) + mmwat * gbsur(igaub)!> We add the contribution to the elemental matrix
     end do
     !
     ! rhs  emmission contribution: u_in
     !
     !tract = -ndots(current_direction_neu) * nitsche_neu * bvnat   ! emiss*bolt4/pi
     tract = -ndots(current_direction_neu) * nitsche_neu * emiss * bsour / pi  
     right =  0.0_rp
     !
     !   If reflection: rho/pi*(...) : DE MOMENTO = ZERO
     !
     if( refle > 0.0_rp ) then
        !
        ! store exiting directions : this ro/2*pi * int_beta wb*ub* (s.n) (beta is outflow)
        ! outgo = list of directions that go out
        !
        out = 0
        do idofn = 1,num_directions_neu
           dotpr = dot_product(direc_neu(1:ndime,idofn),baloc(1:ndime,ndime,igaub))
           if( dotpr > 0.0_rp ) then    ! exiting directions
              ndots(idofn) = dotpr      ! n·s (s:direction of propagation)
              out          = out + 1    ! This stores the number of directions pointing outwards
              outgo(out)   = idofn       ! store exiting directions
           endif
        end do

        if( kfl_smobo_neu == 0 ) then   ! =1 always smooth
           !
           ! Nitsche's boundary conditions (reflective terms)
           !
           gprad = 0.0_rp
           work1 = 0.0_rp
           do jdir = 1,out
              jdofn = outgo(jdir)
              do inodb = 1,pnodb
                 gprad(current_energy_neu,jdofn) = gprad(current_energy_neu,jdofn) + borad(current_energy_neu,jdofn,inodb) * gbsha(inodb,igaub)
              end do
              work1  = work1 + weigd_neu(jdofn) * ndots(jdofn) * gprad(current_energy_neu,jdofn)
           end do
           tract = tract - work1 * refle / pi * ndots(current_direction_neu) * nitsche_neu
        else
           !
           ! Smooth nitsche's boundary conditions
           !
           do inodb = 1,pnodb
              gprad = 0.0_rp
              do jdir = 1,out
                 jdofn        = outgo(jdir)
                 right(inodb) = right(inodb) + weigd_neu(jdofn) * ndots(jdofn) * borad(current_energy_neu,jdofn,inodb) ! integral n.s u ds we do not multiply by gpsur or gpsha because we integrate over the directions, and this integral is already done with the weights of the directions
              end do
              right(inodb) = - right(inodb) * refle / pi * ndots(current_direction_neu) * nitsche_neu      ! integral por reflectividad/ pi nitsche *ndots z
           end do
        end if
     end if
     !
     ! Assembly
     !
     do inodb = 1,pnodb
        inode = lboel(inodb)
        elrhs(inode) = elrhs(inode) + gbsha(inodb,igaub) * (tract+right(inodb)) * gbsur(igaub) ! elemental rhs  vector  due to emissivity (tract) and reflectivity (right)
     end do

  end do gauss

end subroutine neu_boumat
