!------------------------------------------------------------------------
!> @addtogroup alya2pos
!> @{
!> @file    gidres_he.f90
!> @author  Gerard Guillamet
!> @brief   Gauss Point Header using GiD format
!> @details This routine write Gauss Point results using GiD format
!           - Uses:
!             alya2pos
!> @}
!------------------------------------------------------------------------
subroutine gidres_he(iel1,iel2,ndime,ngaus,lexis,flag_coh)

  use def_kintyp, only         : ip,rp,lg,cenam,nnode
  use def_elmtyp

  implicit none

  integer(ip),   intent(in)   :: iel1,iel2
  integer(ip),   intent(in)   :: ndime
  integer(ip),   intent(in)   :: lexis(*)
  integer(ip),   intent(in)   :: ngaus(*)
  logical(lg),   intent(in)   :: flag_coh

  integer(ip)                 :: ielty
  character(13)               :: elemt

  !
  ! Gauss Point Header
  !
  do ielty = iel1,iel2                ! Loop over all element types available in the library
     if( lexis(ielty) /= 0_ip ) then  ! Element types identified
        !
        ! Element type and number of gauss points
        !
        if( ndime == 2_ip ) then
           if( nnode(ielty) == 3_ip .or. nnode(ielty) == 6_ip .or. nnode(ielty) == 7_ip ) then
              elemt = 'Triangle'
           else
              elemt = 'Quadrilateral'
           end if
        else
           if( nnode(ielty) == 4_ip .or. nnode(ielty) == 10_ip ) then
              elemt = 'Tetrahedra'
           else if( nnode(ielty) == 5_ip ) then
              elemt = 'Pyramid'
           else if( nnode(ielty) == 8_ip .or. nnode(ielty) == 20_ip .or. nnode(ielty) == 27_ip ) then
              elemt = 'Hexahedra'
           else if( nnode(ielty) == 6_ip .or. nnode(ielty) == 15_ip ) then
              elemt = 'Prism'
           end if
        end if

        write(11,1) 'GaussPoints '//'GP_'//trim(cenam(ielty))//' Elemtype '//trim(elemt)
        write(11,6) 'Number of Gauss Points: ',ngaus(ielty)
        !
        ! Natural coordinates
        !
        write(11,1) 'Natural Coordinates: given'
        if( ielty == TRI03 ) then
           if( ngaus(ielty) == 1_ip ) then
              write(11,*)  ' 0.33333333333333331       0.33333333333333331'
           else if( ngaus(ielty) == 3_ip ) then
              write(11,*)  ' 0.66666666666666663       0.16666666666666666'
              write(11,*)  ' 0.16666666666666666       0.66666666666666663'
              write(11,*)  ' 0.16666666666666666       0.16666666666666666'
           else
              call runend('TRI3: NOT CODED')
              stop
           end if
           !
           ! QUA04
           !
        else if( ielty == QUA04 ) then
           write(11,*)  '-0.57735026918962595      -0.57735026918962595'
           write(11,*)  '-0.57735026918962595       0.57735026918962595'
           write(11,*)  ' 0.57735026918962595      -0.57735026918962595'
           write(11,*)  ' 0.57735026918962595       0.57735026918962595'
        else if( ielty == TET04 ) then
           write(11,*)  ' 0.13819660112501050       0.13819660112501050       0.13819660112501050'
           write(11,*)  ' 0.58541019662496852       0.13819660112501050       0.13819660112501050'
           write(11,*)  ' 0.13819660112501050       0.58541019662496852       0.13819660112501050'
           write(11,*)  ' 0.13819660112501050       0.13819660112501050       0.58541019662496852'
        else if( ielty == PYR05 ) then
           write(11,*)  '-0.58423739467217717      -0.58423739467217717      -0.66666666666666663'
           write(11,*)  ' 0.58423739467217717      -0.58423739467217717      -0.66666666666666663'
           write(11,*)  ' 0.58423739467217717       0.58423739467217717      -0.66666666666666663'
           write(11,*)  '-0.58423739467217717       0.58423739467217717      -0.66666666666666663'
           write(11,*)  ' 0.0000000000000000        0.0000000000000000        0.40000000000000002'
           !
           ! HEX08
           !
        else if( ielty == HEX08 ) then
           write(11,*)  ' -0.57735026918962595 -0.57735026918962595 -0.57735026918962595'
           write(11,*)  '  0.57735026918962595 -0.57735026918962595 -0.57735026918962595'
           write(11,*)  '  0.57735026918962595  0.57735026918962595 -0.57735026918962595'
           write(11,*)  ' -0.57735026918962595  0.57735026918962595 -0.57735026918962595'
           write(11,*)  ' -0.57735026918962595 -0.57735026918962595  0.57735026918962595'
           write(11,*)  '  0.57735026918962595 -0.57735026918962595  0.57735026918962595'
           write(11,*)  '  0.57735026918962595  0.57735026918962595  0.57735026918962595'
           write(11,*)  ' -0.57735026918962595  0.57735026918962595  0.57735026918962595'
        end if
        write(11,1) 'End GaussPoints'
        write(11,1) ' '
     end if
  end do
  !
  ! Gauss Point Header for Cohesive elements
  !
  if( flag_coh ) then
     write(11,1) 'GaussPoints '//'GP_HEX08_COH '//'Elemtype '//'Hexahedra  '//'coh'
     write(11,6) 'Number of Gauss Points: ',8_ip
     write(11,1) 'Natural Coordinates: given'
     write(11,*) ' -1.0 -1.0 -0.001'
     write(11,*) '  1.0 -1.0 -0.001'
     write(11,*) '  1.0  1.0 -0.001'
     write(11,*) ' -1.0  1.0 -0.001'
     write(11,*) ' -1.0 -1.0  0.001'
     write(11,*) '  1.0 -1.0  0.001'
     write(11,*) '  1.0  1.0  0.001'
     write(11,*) ' -1.0  1.0  0.001'
     write(11,1) 'End GaussPoints'
     write(11,1) ' '
  end if

  !
  ! GiD writting formats
  !
1 format(a)
6 format(a,1x,i2)

end subroutine gidres_he
