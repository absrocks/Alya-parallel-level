!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_addarr.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine computes some additional arrays
!> @details This routine computes some additional arrays
!> @} 
!-----------------------------------------------------------------------
subroutine exm_addarr()
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_addarr
  ! NAME 
  !    exm_addarr
  ! DESCRIPTION
  !   subroutine that computes the conductivity tensor for 2D and 3D 
  ! USES
  !    exm_...
  ! USED BY
  !    exm_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_exmedi
  use      mod_maths,   only: maths_normalize_vector
  implicit none
  integer(ip) :: ipoin,ielem,imate,pelty,pnode,inode
  integer(ip) :: idime,igaus
  real(rp)    :: xauxi
  real(rp)    :: fib_aux(ndime), fib_aux_ortho1(ndime), fib_aux_ortho2(ndime) !! Auxiliary fiber fields
  real(rp)    :: xref_fiber(ndime,3) ! Fiber at the reference frame for the fiber fields
  real(rp)    :: gdiff(3)        ! Conductivity of the fiber
  real(rp)    :: xbalo(3,3)
  real(rp)    :: xnorm(3),xlong,xtra1,xtra2,elcod(ndime,mnode)
  real(rp)    :: elfib(ndime,mnode)
  real(rp)    :: elshe(ndime,mnode),elnor(ndime,mnode) ! Fiber field for orthotropic model
  real(rp)    :: cartd(ndime,mnode)
  real(rp)    :: dvolu,detjm
  real(rp)    :: xjaci(ndime,ndime),xjacm(ndime,ndime),xshap
  logical(lg) :: ortho_model

  ! Initialize arrays
  fib_aux_ortho2(:)= 0.0_rp
  fib_aux(:)= 0.0_rp

  !
  ! Evaluation of orthotropic model
  !
  ortho_model=.false.
  if( (modor_exm(1).ne.0_ip) .and. (modor_exm(2) .ne. 0_ip) ) ortho_model=.true.


  if( INOTEMPTY ) then
     !--------------------------------------------------------!  
     !--------------< INITIALISE FIBER FIELDS >---------------!
     !    
     ! For the isotropic part
     !
     if( modfi_exm < 0 ) then
        fiber_exm => xfiel(-modfi_exm) % a(:,:,1)
     else if ( modfi_exm > 0 ) then
        fib_aux(:)= 0.0_rp
        if(modfi_exm.eq.1_ip) then
          ! X ALIGNED
          fib_aux(1) = 1.0_rp
        elseif(modfi_exm.eq.2_ip) then
          ! Y ALIGNED
          fib_aux(2) = 1.0_rp
        elseif(modfi_exm.eq.3_ip) then
          ! Z ALIGNED
          fib_aux(3) = 1.0_rp
        else
          call runend('EXM_ADDARR: FIBER PREDEFINED OPTION NOT RECOGNISED')
        endif
     else
       call runend('EXM_ADDARR: FIBERS REQUIRED FOR ELECTROPHYSIOLOGY')
     end if
     !
     ! Sheet and normal vectors
     !
     if ( ortho_model ) then
        ! Sheet fibre field
        if(modor_exm(1)<0_ip)then
          sheet_exm => xfiel(-modor_exm(1)) % a(:,:,1)
        elseif(modor_exm(1)>0_ip)then
          fib_aux_ortho1(:)= 0.0_rp
          if(modor_exm(1).eq.1_ip) then
            ! X ALIGNED
            fib_aux_ortho1(1) = 1.0_rp
          elseif(modor_exm(1).eq.2_ip) then
            ! Y ALIGNED
            fib_aux_ortho1(2) = 1.0_rp
          elseif(modor_exm(1).eq.3_ip) then
            ! Z ALIGNED
            fib_aux_ortho1(3) = 1.0_rp
          else
            call runend('EXM_ADDARR: FIBER PREDEFINED OPTION NOT RECOGNISED')
          endif
        endif

        ! Normal fibre field
        if(modor_exm(2)<0_ip)then
            normal_exm => xfiel(-modor_exm(2)) % a(:,:,1)
        elseif(modor_exm(2)>0_ip)then
          fib_aux_ortho2(:)= 0.0_rp
          if(modor_exm(2).eq.1_ip) then
            ! X ALIGNED
            fib_aux_ortho2(1) = 1.0_rp
          elseif(modor_exm(2).eq.2_ip) then
            ! Y ALIGNED
            fib_aux_ortho2(2) = 1.0_rp
          elseif(modor_exm(2).eq.3_ip) then
            ! Z ALIGNED
            fib_aux_ortho2(3) = 1.0_rp
          else
            call runend('EXM_ADDARR: FIBER PREDEFINED OPTION NOT RECOGNISED')
          endif


        endif
     end if
     !----------///////////////////////////////---------------!  
     !----------< END INITIALISE FIBER FIELDS >---------------!

     !--------------------------------------------------------!  
     !--------------< INITIALISE CELL FIELDS >----------------!
     !
     ! CEDIF_EXM: Cell types
     !
     if( kfl_heter_exm == 1_ip ) then
        celty_exm => xfiel(-modce_exm) % a(:,:,1)
     end if
     if( kfl_atbhe_exm == 1_ip ) then
        atbhe_exm => xfiel(-modab_exm) % a(:,:,1)
     end if
     !----------///////////////////////////////---------------!  
     !----------< END INITIALISE CELL FIELDS >----------------!


     !--------------------------------------------------------!  
     !--------------< DIFFUSION TENSOR CEDIF >----------------!

     cedif_exm(:,:,:) = 0.0_rp
 
     elements: do ielem = 1,nelem

     !------------------------------------------!  
     !------< Gather operations >---------------!
 
      pelty = ltype(ielem)
      
      if( pelty > 0 ) then
         pnode = nnode(pelty)
         imate = lmate(ielem)
 
         do inode=1,pnode
            ipoin= lnods(inode,ielem)
            elcod(1:ndime,inode)=coord(1:ndime,ipoin)
         end do
 
         ! Element gather for fibre direction
         if (modfi_exm < 0) then
            do inode=1,pnode
               ipoin= lnods(inode,ielem)
               do idime=1,ndime
                  elfib(idime,inode)=fiber_exm(idime,ipoin)
               end do
            end do
         else if (modfi_exm > 0) then
            do inode=1,pnode
               elfib(1:ndime,inode)=fib_aux(1:ndime)                 
            end do
         end if
 
         if ( ortho_model ) then
            ! Element gather for shear fibre field
            if(modor_exm(1)<0_ip)then
              do inode=1,pnode
                 ipoin= lnods(inode,ielem)
                 do idime=1,ndime
                    elshe(idime,inode)=sheet_exm(idime,ipoin)
                 end do
              end do
            elseif(modor_exm(1)>0_ip)then
              do inode=1,pnode
                 elshe(1:ndime,inode)=fib_aux_ortho1(1:ndime)                 
              end do
            endif
            ! Element gather for normal fibre field
            if(modor_exm(2)<0_ip)then
              do inode=1,pnode
                 ipoin= lnods(inode,ielem)
                 do idime=1,ndime
                    elnor(idime,inode)=normal_exm(idime,ipoin)
                 end do
              end do
            elseif(modor_exm(2)>0_ip)then
              do inode=1,pnode
                 elnor(1:ndime,inode)=fib_aux_ortho2(1:ndime)                 
              end do
            endif
 
         end if
 
     !------< END gather operations >-----------!
     !-------\\\\\\\\\\\\\\\\\\\\\\-------------!  
 
         gauss_points: do igaus=1,ngaus(pelty)
 
            ! Cartesian derivatives and Jacobian
            call elmder(pnode,ndime,elmar(pelty)%deriv(1,1,igaus),elcod,cartd,detjm,xjacm,xjaci)
            dvolu=elmar(pelty)%weigp(igaus)*detjm
 
            xref_fiber=0.0_rp
            do inode=1,pnode
               xshap= elmar(pelty)%shape(inode,igaus)
               do idime= 1,ndime
                  xref_fiber(idime,1)= xref_fiber(idime,1)+elfib(idime,inode)*xshap
 
                  if ( ortho_model ) then
                     xref_fiber(idime,2)= xref_fiber(idime,2)+elshe(idime,inode)*xshap
                     xref_fiber(idime,3)= xref_fiber(idime,3)+elnor(idime,inode)*xshap
                  end if
               end do
            end do
            
            gdiff(1) = gdiff_exm(1,1,imate)  ! longitudinal intra
            gdiff(2) = gdiff_exm(1,2,imate)  ! cross wise 2
            if(ndime.eq.3) gdiff(3) = gdiff_exm(1,3,imate)  ! cross wise 3
 
 
             xauxi=sqrt(dot_product(xref_fiber(1:ndime,1),xref_fiber(1:ndime,1)))
 
            xbalo(1,1) = xref_fiber(1,1)/xauxi
            xbalo(2,1) = xref_fiber(2,1)/xauxi
            if (ndime==3) xbalo(3,1) = xref_fiber(3,1)/xauxi
 
 
            !
            ! Code for 2D problem
            ! WARNING: this code is not functional
            !
            if( ndime == 2 ) then
               call runend('EXM_ADDARR: EXMEDI NOT PREPARED FOR 2D PROBLEMS')
               ! Evaluate perpendicular vector (:,2) cross-fibre vector c
               xbalo(1,2)=   - xbalo(2,1)
               xbalo(2,2)=     xbalo(1,1)
               
               do inode=1,pnode
                  ipoin= lnods(inode,ielem)
                  xshap= elmar(pelty)%shape(inode,igaus)
                  ! Evaluate xbalo (a) * diffusion tensor (g) * xbalo (a)
                  ! f - fibre, c- cross-fibre
 
                  ! g_f * a_fx**2 + g_c * a_cx**2
                  cedif_exm(1,1,ipoin) = cedif_exm(1,1,ipoin) + &
                        xshap * (gdiff(1)*xbalo(1,1)*xbalo(1,1) &
                        + gdiff(2)*xbalo(1,2)*xbalo(1,2) ) * dvolu / vmass(ipoin)
 
                  ! g_f * a_fy**2 + g_c * a_cy**2
                  cedif_exm(2,2,ipoin) = cedif_exm(2,2,ipoin) + &
                        xshap * (gdiff(1)*xbalo(2,1)*xbalo(2,1) &
                       + gdiff(2)*xbalo(2,2)*xbalo(2,2) ) * dvolu / vmass(ipoin)
 
                  ! g_f * a_fx*a_fy + g_c * a_cx*a_cy
                  cedif_exm(1,2,ipoin) = cedif_exm(1,2,ipoin) + &
                        xshap * (gdiff(1)*xbalo(1,1)*xbalo(2,1) &
                        + gdiff(2)*xbalo(1,2)*xbalo(2,2) ) * dvolu / vmass(ipoin)
                  
                  cedif_exm(2,1,ipoin) = cedif_exm(1,2,ipoin)
                  
               end do
 
            else if (ndime==3) then
 
                  do inode= 1,pnode
 
                     ipoin= lnods(inode,ielem)
                     xshap= elmar(pelty)%shape(inode,igaus)
 
                     if (.not. ortho_model) then
                       !
                       ! Transversely isotropic diffusion
                       !
                       ! Evaluate arbitrary cross-fibre vector xbalo(:,2)
                       xnorm(1)=0.0_rp
                       xnorm(2)=1.0_rp
                       xnorm(3)=0.0_rp
 
                       xbalo(1,2) = xnorm(2)*xbalo(3,1)-xnorm(3)*xbalo(2,1)
                       xbalo(2,2) = xnorm(3)*xbalo(1,1)-xnorm(1)*xbalo(3,1)
                       xbalo(3,2) = xnorm(1)*xbalo(2,1)-xnorm(2)*xbalo(1,1)
 
                       xauxi=sqrt(xbalo(1,2)*xbalo(1,2)+xbalo(2,2)*xbalo(2,2)+xbalo(3,2)*xbalo(3,2))
 
                       if (xauxi < 0.000001_rp) then
                          xnorm(1) = 1.0_rp
                          xnorm(2) = 0.0_rp
                          xnorm(3) = 0.0_rp
 
                          xbalo(1,2) = xnorm(2)*xbalo(3,1)-xnorm(3)*xbalo(2,1)
                          xbalo(2,2) = xnorm(3)*xbalo(1,1)-xnorm(1)*xbalo(3,1)
                          xbalo(3,2) = xnorm(1)*xbalo(2,1)-xnorm(2)*xbalo(1,1)
 
                          xauxi = sqrt(xbalo(1,2)*xbalo(1,2)+xbalo(2,2)*xbalo(2,2)+xbalo(3,2)*xbalo(3,2))
                       end if
 
                       xbalo(1,2) = xbalo(1,2)/xauxi
                       xbalo(2,2) = xbalo(2,2)/xauxi
                       xbalo(3,2) = xbalo(3,2)/xauxi
 
                       ! axial vector
                       xbalo(1,3) = xbalo(2,1)*xbalo(3,2)-xbalo(3,1)*xbalo(2,2)
                       xbalo(2,3) = xbalo(3,1)*xbalo(1,2)-xbalo(1,1)*xbalo(3,2)
                       xbalo(3,3) = xbalo(1,1)*xbalo(2,2)-xbalo(2,1)*xbalo(1,2)
 
                       xlong= gdiff(1)
                       xtra1= gdiff(2)
                       xtra2= xtra1
 
                  else
                     !
                     ! Orthotropic diffusion
                     !
                     ! Normalise fiber vector xbalo(:,1)
                     xauxi=sqrt(dot_product(xref_fiber(1:ndime,1),xref_fiber(1:ndime,1)))

                     xbalo(1,1) = xref_fiber(1,1)/xauxi
                     xbalo(2,1) = xref_fiber(1,2)/xauxi
                     xbalo(3,1) = xref_fiber(1,3)/xauxi

                     ! Sheet vector is read-in
                     xauxi=sqrt(dot_product(xref_fiber(1:ndime,2),xref_fiber(1:ndime,2)))
                     xbalo(1,2) = xref_fiber(2,1)/xauxi
                     xbalo(2,2) = xref_fiber(2,2)/xauxi
                     xbalo(3,2) = xref_fiber(2,3)/xauxi

                     ! Sheet-normal vector is read-in
                     xauxi=sqrt(dot_product(xref_fiber(1:ndime,3),xref_fiber(1:ndime,3)))
                     xbalo(1,3) = xref_fiber(3,1)/xauxi
                     xbalo(2,3) = xref_fiber(3,2)/xauxi
                     xbalo(3,3) = xref_fiber(3,3)/xauxi
 
                     xlong= gdiff(1)
                     xtra1= gdiff(2)
                     xtra2= gdiff(3)
 
                  end if
 
                    ! Evaluate diffusion tensor
                    cedif_exm(1,1,ipoin) = cedif_exm(1,1,ipoin) + &
                         xshap*(xlong*xbalo(1,1)*xbalo(1,1) + &
                         xtra1*xbalo(1,2)*xbalo(1,2) + &
                         xtra2*xbalo(1,3)*xbalo(1,3)) * dvolu / vmass(ipoin)                
                    cedif_exm(2,2,ipoin) = cedif_exm(2,2,ipoin) +&
                         xshap*(xlong*xbalo(2,1)*xbalo(2,1) + &
                         xtra1*xbalo(2,2)*xbalo(2,2) + &
                         xtra2*xbalo(2,3)*xbalo(2,3)) * dvolu / vmass(ipoin) 
                    cedif_exm(3,3,ipoin) = cedif_exm(3,3,ipoin) + &
                         xshap*(xlong*xbalo(3,1)*xbalo(3,1) + &
                         xtra1*xbalo(3,2)*xbalo(3,2) + &
                         xtra2*xbalo(3,3)*xbalo(3,3)) * dvolu / vmass(ipoin)              
                    cedif_exm(1,2,ipoin) = cedif_exm(1,2,ipoin) + &
                         xshap*(xlong*xbalo(1,1)*xbalo(2,1) + &
                         xtra1*xbalo(1,2)*xbalo(2,2) + &
                         xtra2*xbalo(1,3)*xbalo(2,3)) * dvolu / vmass(ipoin)              
                    cedif_exm(1,3,ipoin) = cedif_exm(1,3,ipoin) + &
                         xshap*(xlong*xbalo(1,1)*xbalo(3,1) + &
                         xtra1*xbalo(1,2)*xbalo(3,2) + &
                         xtra2*xbalo(1,3)*xbalo(3,3)) * dvolu / vmass(ipoin)              
                    cedif_exm(2,3,ipoin) = cedif_exm(2,3,ipoin) + &
                         xshap*(xlong*xbalo(2,1)*xbalo(3,1) + &
                         xtra1*xbalo(2,2)*xbalo(3,2) + &
                         xtra2*xbalo(2,3)*xbalo(3,3)) * dvolu / vmass(ipoin)              
 
                    cedif_exm(3,2,ipoin) = cedif_exm(2,3,ipoin)
                    cedif_exm(3,1,ipoin) = cedif_exm(1,3,ipoin)
                    cedif_exm(2,1,ipoin) = cedif_exm(1,2,ipoin)
                  
               end do
            end if
           end do gauss_points
           
        end if
     end do elements
 
    call rhsmod(ndime*ndime,cedif_exm)

  end if  

  
end subroutine exm_addarr
