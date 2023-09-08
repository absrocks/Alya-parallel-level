subroutine opt_upddes
  !------------------------------------------------------------------------
  !****f* Optsol/opt_upddes
  ! NAME
  !    opt_upddes
  ! DESCRIPTION
  ! OUTPUT
  ! USED BY
  !    Optsol
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_optsol
  use      def_inpout
  use      def_solver
  implicit none

  integer(ip) :: indvars, isInside
  real(rp) :: maxabsdescdir, minabsdesign, minabsdescdir

  scalefactor   =0.0_rp
  maxabsdescdir =-1.0E+30_rp
  minabsdescdir = 1.0E+30_rp
  minabsdesign  = 1.0E+30_rp

  do indvars=1,kfl_ndvars_opt

     isInside=diffj_isInside(indvars)
     if(isInside>0_ip .and. descdir(indvars)/=0.0_rp )then
        maxabsdescdir = max(abs(descdir(indvars)),maxabsdescdir)
        minabsdescdir = min(abs(descdir(indvars)),minabsdescdir)
        minabsdesign  = min(abs(design_vars(indvars)),minabsdesign)
     end if
  end do

  call pararr('MAX',0_ip,1_ip,maxabsdescdir) 
  call pararr('MIN',0_ip,1_ip,minabsdescdir) 
  call pararr('MIN',0_ip,1_ip,minabsdesign) 

  scalefactor = kfl_scale_opt*(minabsdesign/ (0.5_rp*(minabsdescdir+maxabsdescdir)))

  design_vars_tmp(:) = 0.0_rp

  if(kfl_first_opt==1) then
     design_vars_prev(:) = 0.0_rp
  else 
     if(kfl_curlin_opt==1)then 
        design_vars_prev(:) = design_vars(:)
     end if
  end if

  if(kfl_curlin_opt==1) then
     do indvars=1,kfl_ndvars_opt
        design_vars_tmp(indvars) = design_vars(indvars) + stepj * scalefactor * descdir(indvars)
     end do
     diffj_prev(:)   = diffj(:)
     descdir_prev(:) = descdir(:)
  else
     do indvars=1,kfl_ndvars_opt
        design_vars_tmp(indvars) = design_vars_prev(indvars) + stepj * scalefactor * descdir_prev(indvars)
     end do
  end if
 
end subroutine opt_upddes   
