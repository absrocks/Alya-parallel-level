subroutine maxmax(ndim1,ndim2,ndim3,vecto,vemax)
  !------------------------------------------------------------------------
  !****f* mathru/maxmax
  ! NAME 
  !    minmax
  ! DESCRIPTION
  !    Compute the maximum of a vector 
  ! USES
  ! USED BY
  !    Modules: *_cvgunk
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master
  implicit none
  integer(ip), intent(in)  :: ndim1,ndim2,ndim3
  real(rp),    intent(in)  :: vecto(ndim1,ndim2)
  real(rp),    intent(out) :: vemax
  integer(ip)              :: idim1,idim2,ndim4
  real(rp)                 :: uvalu
  real(rp),    target      :: vtmax(1)

  vemax=-1.0e30

  if(kfl_paral==-1.or.(kfl_paral==0.and.kfl_ptask==0)) then
     ! 
     ! Sequential case
     !
     if(ndim1==1) then
        do idim2=1,ndim2
           if(vecto(1,idim2)>vemax) vemax=vecto(1,idim2)
        end do

     else if(ndim3<0) then
        ndim4=-ndim3
        do idim2=1,ndim2
           if(vecto(ndim4,idim2)>vemax) vemax=vecto(ndim4,idim2)
        end do

     else
        ndim4=min(ndim1,ndim3)
        do idim2=1,ndim2
           uvalu=0.0_rp
           do idim1=1,ndim4
              uvalu=uvalu+vecto(idim1,idim2)*vecto(idim1,idim2)
           end do
           uvalu=sqrt(uvalu)
           if(uvalu>vemax) vemax=uvalu
        end do
     end if

  else
     !
     ! Parallel case
     !
     if(kfl_paral>=1) then
        if(ndim1==1) then
           do idim2=1,ndim2
              if(vecto(1,idim2)>vemax) vemax=vecto(1,idim2)
           end do
        else
           ndim4=min(ndim1,ndim3)
           do idim2=1,ndim2
              uvalu=0.0_rp
              do idim1=1,ndim4
                 uvalu=uvalu+vecto(idim1,idim2)*vecto(idim1,idim2)
              end do
              uvalu=sqrt(uvalu)
              if(uvalu>vemax) vemax=uvalu
           end do
        end if
     end if
     !
     ! Maximum
     ! 
     vtmax(1) =  vemax
     nparr    =  1
     parre    => vtmax
     call par_operat(2_ip)
     vemax    =  vtmax(1)
    
  end if

end subroutine maxmax
