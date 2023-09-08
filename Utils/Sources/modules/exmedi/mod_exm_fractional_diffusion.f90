!-----------------------------------------------------------------------
!
!> @addtogroup Exmedi
!> @{
!> @name    Fractional difussion for electrophisiology
!> @file    mod_exm_fractional_diffusion.f90
!> @author  Alfonso Santiago
!> @brief   Fractional diffusion 
!> @details Fractional diffusion for electrophysiology. Based on
!>          Nicole Cusimano's mathematical developement and 
!>          Mathlab code.
!> @{
!
!-----------------------------------------------------------------------
module mod_exm_fractional_diffusion

!  use def_master
  use def_kintyp, only : ip, rp, lg
  use def_parame, only : pi

  implicit none
  real(rp),  parameter :: epsil = epsilon(1.0_rp)

  real(rp), save                                :: MSeigenmax, &
                                                   MSeigenmin

  real(rp), save, dimension(:), allocatable  :: quadr_nodes(:)
  real(rp), save, dimension(:), allocatable     :: quadr_weigh(:)
 
  private

  interface jacobi_elliptic_function
    module procedure  jacobi_elliptic_scalar, &
                      jacobi_elliptic_matrix
  end interface


   public :: elliptic_integral
   public :: jacobi_elliptic_function

   contains
    !-----------------------------------------------------------------------
    !> 
    !> @author  Alfonso Santiago
    !> @date    2019-FEB-18
    !> @brief   Initialisation for the fractional diffusion
    !> @details This subroutine initialises the fractional diffusion
    !>          for electrophysiology. The objective is to obtain
    !>                 MSeigenmax  : the largest eigenvalue of M{^1}S
    !>                 MSeigenmin  : the smalles eigenvalue of M{-1}S
    !>                 quadr_nodes : the quadrature nodes
    !>                 quadr_weigh : the quadrature weights
    !>         
    !>         
    !-----------------------------------------------------------------------
    subroutine initialisation(nintp)

        implicit none
        integer(ip), intent(in)                 :: nintp ! Number of INTegration Points
        complex(rp), dimension(:), allocatable  :: midpn
        real(rp), dimension(:), allocatable     :: sn, cn, dn, dxidt
        integer(ip)                             :: k, i
        real(rp)                                :: ellipK, ellipKp

        allocate(quadr_nodes(nintp))
        allocate(quadr_weigh(nintp))
        allocate(midpn(nintp))
        allocate(sn(nintp))
        allocate(cn(nintp))
        allocate(dn(nintp))
        allocate(dxidt(nintp))


        ! These two values should be obtained for the M^{-1}S matrix
        MSeigenmax=1.0_rp
        MSeigenmin=0.0_rp
        
        k=(sqrt(MSeigenmax/MSeigenmin)-1.0_rp)/(sqrt(MSeigenmax/MSeigenmin)-1.0_rp)
        
        call elliptic_integral(-log(real(k))/pi,ellipK,ellipKp)


        do i=1,nintp
          midpn(i)=nintp-0.5_rp-i
          midpn(i) = 0.5_rp*cmplx(0.0,1.0,kind=rp)*ellipKp - ellipK + midpn(i)*2.0_rp*K/nintp
        enddo


        do  i=1,nintp
           call jacobi_elliptic_function(midpn(i), -log(real(k))/pi, sn(i), cn(i), dn(i))
           quadr_nodes(i) = sqrt(MSeigenmax*MSeigenmin)*(1/k+sn(i))*(1/k-sn(i))

           dxidt(i)=cn(i)*dn(i)/(1/k-sn(i))**2.0_rp
           quadr_weigh(i)=f(quadr_nodes(i))*dxidt(i)
        enddo

        contains
            function f(z) result(r)
              implicit none
              real(rp), intent(in)    :: z
              real(rp)                :: r
              real(rp)                :: dt, coeff

              !! dt is not defined!!
              !! What is coefff?
              call runend("dt and coeff not defined")
              r = 1.0_rp/(1.0_rp+dt*coeff*z)

            endfunction

    end subroutine initialisation
    !-----------------------------------------------------------------------
    !> 
    !> @author  Alfonso Santiago
    !> @date    2019-JAN-30
    !> @brief   Complete eliptic integral of the first kind, with complement
    !> @details Returns the value of the complete elliptic integral of the 
    !>          first kind, evaluated at M=exp(-2*pi*L), 0<L<inf
    !>          Also returns the result for the complementary parameter 1-M
    !>          which is useful when M<EPS. Even when M<1e-6.
    !>          Recall that the elliptic modulus k is related to the parameter
    !>          M by M=k**2
    !>         
    !>         This subrutine uses the method of the arithmetic-geometric 
    !>         mean described in 17.6 of M. Abramowitz and IA Stegun
    !>         "handbook of mathematical functions", Dover ,1965.
    !>         
    !>         This subroutine was adapted  from Toby Driscoll's
    !>         ellipkkp.m.        
    !>         
    !>         
    !-----------------------------------------------------------------------
    subroutine elliptic_integral(L,K,Kp)

    implicit none
    real(rp), intent(in)       :: L !Must be a scalar
    real(rp), intent(out)      :: K, Kp
    real(rp)                   :: m

      
     K=0.0_rp
     Kp=0.0_rp

    if(L .gt. 10.0_rp) then
    !If m=exp(e,-2*pi*L) is extremely small, use 0(m) approximations
      K=pi/2.0_rp
      Kp=pi*L+log(4.0_rp)

    else
    ! If  not, do the proper computations

      m=exp(-2.0_rp*pi*L)
      call computeK(m,K,'GETKK')

      call computeK(m,Kp,'GETKP')

    endif
    

    contains
     
        subroutine computeK(m_in,K_out,option)
            implicit none
            real(rp), intent(in)  :: m_in
            real(rp), intent(out) :: K_out
            character(len=5)      :: option
            real(rp)              :: a0, b0,s0, a1, b1, c1 , mm, w1, m_comp
            real(rp),parameter    :: eps=1e-10
            integer(ip)           ::i1

            !! Initialisation depending on the case
            if(option .eq. 'GETKK') then
              a0=1.0_rp
              b0=sqrt(1.0_rp-m_in)
              s0=m_in
            elseif(option .eq. 'GETKP') then
              a0=1.0_rp
              b0=sqrt(m_in)
              s0=1.0_rp-m_in
            else
              call runend ('MOD_EXMEDI: SUBRU ELLIPTIC: OPTION NO DEFINED')
            endif

            i1=0_ip
            mm=1.0_rp
            K_out=0.0_rp

            !! Loop independent of the case
            !
            do while(mm .gt. eps)
              a1=(a0+b0)/2.0_rp
              b1=sqrt(a0*b0)
              c1=(a0-b0)/2.0_rp
              i1=i1+1_ip
              w1=(2.0_rp**i1)*(c1**2.0_rp)

              !! Code adaptation as w1 is not a matrix/vector
              !
              !mm=maxval(maxval(w1))
              mm=abs(w1)
              !
              !!
              s0=s0+w1
              a0=a1
              b0=b1
            end do
            
            K_out = pi/(2*a1)
            
            !! Ending, dependent on the case
            !
            if(option .eq. 'GETKK') then
              m_comp=1.0_rp
            elseif(option .eq. 'GETKP') then
              m_comp=0.0_rp
            else
              call runend ('MOD_EXMEDI: SUBRU ELLIPTIC: OPTION NO DEFINED')
            endif
            
            !! Code adaptation as m_in is not a vector
            !! And we dont have find and isempty functions
            !
            !im = find(m_in==m_comp)
            !if( the vector IM is not empty ) then
              !K_out(im) = K_out(im)*huge(K_out(im))
            !endif
            if(m_in .eq. m_comp) then
              K_out = K_out*huge(K_out)
            endif
            !
            !!

        end subroutine computeK

    end subroutine elliptic_integral

    !-----------------------------------------------------------------------
    !> 
    !> @author  Alfonso Santiago
    !> @date    2019-JAN-30
    !> @brief   Jacobi Elliptic functions for complex arguments
    !> @details Returns the values of the Jacobi elliptic functions evaluated
    !>          at complex argument U, and parameter M=exp(-2*pi*L), 0<L<inf.
    !>          Recall that M=k**2 where k is thee elliptic modulus.
    !>          U may be a matrix, L must be a scalar. The entries of U
    !>          are expected to lie within the rectangle |Re(u)|<K,
    !>          0<Im(U)<Kp, where K and Kp were computed from the
    !>          elliptic integral.
    !>         
    !>         This algorithm is the descending Landen transformation
    !>         described in L. Howell's PhD thesis from MIT. Additional
    !>         formulas from Gradshteyn and Ryzhik, 5th ed., and  
    !>         Abbramowitz and Stegun
    !>         
    !>         This subroutine was adapted  from Toby Driscoll's
    !>         ellipjc.m.        
    !>         
    !>         
    !-----------------------------------------------------------------------

    recursive subroutine jacobi_elliptic_matrix(U, L , sn, cn, dn, flag)
      implicit none
      complex(rp), intent(inout)         :: U(:,:) ! U may be a matrix (of complex values)
      real(rp), intent(in)               :: L      ! Must be a scalar
      logical(lg), optional, intent(in)  :: flag
     
      real(rp), intent(out)              :: sn(:,:), cn(:,:), dn(:,:) !No allocatable needed

      complex(rp), allocatable           :: v(:,:)
      real(rp), allocatable              :: sn1(:,:),cn1(:,:),dn1(:,:), denom(:,:)

      real(rp),parameter                 :: eps=1e-6
      real(rp)                           :: m
      real(rp)                           :: K,Kp
      logical(lg), allocatable           :: high(:,:)
      integer(ip)                        :: ud1,ud2,i,j
      real(rp)                           :: kappa, x, mu


      sn=0.0_rp
      cn=0.0_rp
      dn=0.0_rp

      ud1=size(U,1)
      ud2=size(U,2)
      allocate(high (ud1,ud2))
      allocate(sn1 (ud1,ud2))
      allocate(cn1 (ud1,ud2))
      allocate(dn1 (ud1,ud2))
      allocate(v (ud1,ud2))
      allocate(denom (ud1,ud2))
      high=.false.

      if(present(flag))then
          high=.false.
          m=L
      else
          call elliptic_integral(L,K,Kp)
             
          do i = 1, ud1
             do j = 1, ud2
               if(aimag(u(i,j)).gt.Kp/2.0_rp) then
                 high(i,j) = .true.
                 u(i,j)=cmplx(0.0,1.0,kind=rp)*Kp-u(i,j)
               endif
             enddo
          enddo
          m=exp(-2.0_rp*pi*L)

      endif


      if (m.lt.4.0_rp*eps) then
          do i = 1, ud1
             do j = 1, ud2
               sn=sin(u(i,j))+m/4.0_rp*(sin(u(i,j))*cos(u(i,j))-u)*cos(u(i,j))
               cn=cos(u(i,j))+m/4.0_rp*(-sin(u(i,j))*cos(u(i,j))+u)*sin(u(i,j))
               dn=1.0_rp+m/4*(cos(u(i,j))**2.0_rp-sin(u(i,j))**2.0_rp-1.0_rp)
             enddo
          enddo

      else
          if (m.gt.1e-3) then
            kappa = (1-sqrt(1-m))/(1+sqrt(1-m))
          else
            x=m/4.0_rp
            kappa = 132.0_rp*x**6.0_rp + 42.0_rp*x**5.0_rp + &  
                     14.0_rp*x**4.0_rp +  5.0_rp*x**3.0_rp + &  
                     2.0_rp*x**2.0_rp +  1.0_rp*x**1.0_rp
          endif
        
          mu=kappa**2.0_rp

          do i = 1, ud1
             do j = 1, ud2
                v(i,j)=u(i,j)/(1.0_rp+kappa)
             enddo
          enddo

          call jacobi_elliptic_matrix(v, mu , sn1, cn1, dn1, .true.)

          do i = 1, ud1
             do j = 1, ud2
                  denom(i,j) = (1.0_rp+kappa*sn1(i,j)**2.0_rp)
                  sn(i,j) = (1.0_rp+kappa)*sn1(i,j)/denom(i,j)
                  cn(i,j) = cn1(i,j)*dn1(i,j)/denom(i,j)
                  dn(i,j) = (1.0_rp-kappa*sn1(i,j)**2)/denom(i,j)
             enddo
          enddo
      endif


      do i = 1, ud1
         do j = 1, ud2
           if(high(i,j))then
             dn(i,j) = cmplx(0.0,1.0,kind=rp)*cn(i,j)/sn(i,j)
             cn(i,j) = cmplx(0.0,1.0,kind=rp)*dn(i,j)/(sqrt(m)*sn(i,j))
             sn(i,j) = -1.0_rp/sqrt(m)*sn(i,j)
           endif
         enddo
      enddo
      




    end subroutine jacobi_elliptic_matrix

    !-----------------------------------------------------------------------
    !> 
    !>
    !>  SCALAR VERSION OF THE SUBROUTINE
    !>
    !> @author  Alfonso Santiago
    !> @date    2019-JAN-30
    !> @brief   Jacobi Elliptic functions for complex arguments
    !> @details Returns the values of the Jacobi elliptic functions evaluated
    !>          at complex argument U, and parameter M=exp(-2*pi*L), 0<L<inf.
    !>          Recall that M=k**2 where k is thee elliptic modulus.
    !>          U may be a matrix, L must be a scalar. The entries of U
    !>          are expected to lie within the rectangle |Re(u)|<K,
    !>          0<Im(U)<Kp, where K and Kp were computed from the
    !>          elliptic integral.
    !>         
    !>         This algorithm is the descending Landen transformation
    !>         described in L. Howell's PhD thesis from MIT. Additional
    !>         formulas from Gradshteyn and Ryzhik, 5th ed., and  
    !>         Abbramowitz and Stegun
    !>         
    !>         This subroutine was adapted  from Toby Driscoll's
    !>         ellipjc.m.        
    !>         
    !>         
    !-----------------------------------------------------------------------

    recursive subroutine jacobi_elliptic_scalar(U, L , sn, cn, dn, flag)
      implicit none
      complex(rp), intent(inout)         :: U ! U is a complex scalar
      real(rp), intent(in)               :: L      ! Must be a scalar
      logical(lg), optional, intent(in)  :: flag
     
      real(rp), intent(out)              :: sn, cn, dn
      complex(rp)                        :: v
      real(rp)                           :: sn1,cn1,dn1, denom

      real(rp),parameter                 :: eps=1e-10
      real(rp)                           :: m
      real(rp)                           :: K,Kp
      logical(lg)                        :: high
      real(rp)                           :: kappa, x, mu


      sn=0.0_rp
      cn=0.0_rp
      dn=0.0_rp
      high=.false.

      if(present(flag))then
          high=.false.
          m=L
      else
          call elliptic_integral(L,K,Kp)
             
          if(aimag(u).gt.Kp/2.0_rp) then
            high = .true.
            u=cmplx(0.0,1.0,kind=rp)*Kp-u
          endif
          m=exp(-2.0_rp*pi*L)
      endif


      if (m.lt.4.0_rp*eps) then
          sn=sin(u)+m/4.0_rp*(sin(u)*cos(u)-u)*cos(u)
          cn=cos(u)+m/4.0_rp*(-sin(u)*cos(u)+u)*sin(u)
          dn=1.0_rp+m/4*(cos(u)**2.0_rp-sin(u)**2.0_rp-1.0_rp)

      else
          if (m.gt.1e-3) then
            kappa = (1-sqrt(1-m))/(1+sqrt(1-m))
          else
            x=m/4.0_rp
            kappa = 132.0_rp*x**6.0_rp + 42.0_rp*x**5.0_rp + &  
                     14.0_rp*x**4.0_rp +  5.0_rp*x**3.0_rp + &  
                     2.0_rp*x**2.0_rp +  1.0_rp*x**1.0_rp
          endif
        
          mu=kappa**2.0_rp
          v=u/(1.0_rp+kappa)

          call jacobi_elliptic_scalar(v, mu , sn1, cn1, dn1, .true.)

          denom = (1.0_rp+kappa*sn1**2.0_rp)
          sn = (1.0_rp+kappa)*sn1/denom
          cn = cn1*dn1/denom
          dn = (1.0_rp-kappa*sn1**2)/denom
      endif


      if(high)then
        dn = cmplx(0.0,1.0,kind=rp)*cn/sn
        cn = cmplx(0.0,1.0,kind=rp)*dn/(sqrt(m)*sn)
        sn = -1.0_rp/sqrt(m)*sn
      endif

    end subroutine jacobi_elliptic_scalar
end module mod_exm_fractional_diffusion
!> @}
