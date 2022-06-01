!------------------------------------------------------------------------
!> @addtogroup alya2pos
!> @{
!> @file    gidres_gp.f90
!> @author  Gerard Guillamet
!> @brief   Gauss Point Result file using GiD format
!> @details This routine write Gauss Point results using GiD format
!           - Uses:
!             alya2pos
!> @}
!------------------------------------------------------------------------
subroutine gidres_gp(iel1,iel2,result_name,step_value,ndime,pdime, &
     nelem_2,mgaus,ngaus,lexis,ltype,leinv,lelch,flag_coh,geve4)

  use def_kintyp, only         : ip,rp,lg,cenam,nnode
  use def_elmtyp

  implicit none

  integer(ip),   intent(in)   :: iel1,iel2,nelem_2
  integer(ip),   intent(in)   :: ndime,pdime,mgaus
  integer(ip),   intent(in)   :: lexis(*)
  integer(ip),   intent(in)   :: ltype(*)
  integer(ip),   intent(in)   :: leinv(*)
  integer(ip),   intent(in)   :: ngaus(*)
  integer(ip),   intent(in)   :: lelch(*)
  real(rp),      intent(in)   :: geve4(ndime,mgaus,nelem_2)
  real(rp),      intent(in)   :: step_value
  character(5),  intent(in)   :: result_name
  logical(lg),   intent(in)   :: flag_coh

  integer(ip)                 :: ielty,ielem
  integer(ip)                 :: igaus
  character(13)               :: elemt

  !
  ! Write Values (Result block)
  !
  do ielty = iel1,iel2                 ! Loop over all element types available
     if( lexis(ielty) /= 0 ) then      ! Element types identified

        if( pdime == 1 ) then          ! R3P
           if( flag_coh ) then
              write(11,2) result_name,'ALYA',step_value,'Scalar','GP_'//trim(cenam(ielty))//'_COH'
           else
              write(11,2) result_name,'ALYA',step_value,'Scalar','GP_'//trim(cenam(ielty))
           end if

        else                           ! R3PVE
           write(11,2) result_name,'ALYA',step_value,'Vector','GP_'//trim(cenam(ielty))
           write(11,4) result_name//'_X,'//result_name//'_Y,'//result_name//'_Z'
        end if

        write(11,1) 'Values'
        if( flag_coh ) then            ! Cohesive elements
           do ielem = 1,nelem_2
              if( abs(ltype(ielem)) == ielty .and. lelch(ielem) == 17_ip) then
                 write(11,5) leinv(ielem),geve4(1:pdime,1,ielem)
                 do igaus = 2,ngaus(ielty)/2_ip
                    write(11,7) geve4(1:pdime,igaus,ielem)
                 end do
                 do igaus = 1,ngaus(ielty)/2_ip
                    write(11,7) geve4(1:pdime,igaus,ielem)
                 end do
              end if
           end do
        else
           do ielem = 1,nelem_2
              if( abs(ltype(ielem)) == ielty) then
                 write(11,5) leinv(ielem),geve4(1:pdime,1,ielem)
                 do igaus = 2,ngaus(ielty)
                    write(11,7) geve4(1:pdime,igaus,ielem)
                 end do
              end if
           end do
        end if
        write(11,1) 'End Values'
        write(11,1) ' '

     end if
  end do
  !
  ! GiD writting formats
  !
1 format(a)
2 format('Result ',a,' ',a,' ',e15.8,' ',a,' OnGaussPoints ',a)
4 format('ComponentNames ',a,a,a)
5 format(i11, 3(1x,e16.8E3))
7 format(3(1x,e16.8E3))

end subroutine gidres_gp
