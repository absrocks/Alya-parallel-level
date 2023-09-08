module mod_por_interp

contains

  !------------------------------------------------------------------------
  !> @addtogroup Porous 
  !> @{
  !> @file    por_gtabl2.f90
  !> @date    13/05/2013
  !> @author  Herbert Owen
  !> @brief   Get Pbh from table
  !> @details Get Pbh from table
  !> @} 
  !------------------------------------------------------------------------

  subroutine por_interp(ntime,nposi,coord,cutim,xvalu,table_value,table_coord)
    use def_kintyp, only       :  ip,rp
    implicit none 
    integer(ip), intent(in)           :: ntime                       !> Number of time steps
    integer(ip), intent(in)           :: nposi                       !> Number of positions
    real(rp),    intent(in)           :: coord                       !> Current coordinate
    real(rp),    intent(in)           :: cutim                       !> Current time
    real(rp),    intent(out)          :: xvalu                       !> current time
    real(rp),    intent(in)           :: table_value(ntime,nposi+1)  !> Table time values
    real(rp),    intent(in), optional :: table_coord(nposi)          !> table coord values
    integer(ip)                       :: itim1,itim2,itime,jtime
    integer(ip)                       :: ipos1,ipos2,iposi,jposi
    real(rp)                          :: rtime,rcoor
    real(rp)                          :: x1,x2,t1,t2,deltat,deltax
    real(rp)                          :: v_x1_t1,v_x1_t2
    real(rp)                          :: v_x2_t1,v_x2_t2
    real(rp)                          :: v_x_t1,v_x_t2
    !
    ! Time interpolation
    !
    if( cutim <= table_value(1,1) ) then
       itim1 = 1
       itim2 = 1
    else if( cutim >= table_value(ntime,1) ) then
       itim1 = ntime
       itim2 = ntime
    else
       itime = 0
       jtime = 0
       do while( itime < ntime )
          itime = itime + 1
          rtime = table_value(itime,1)
          if( rtime > cutim ) then
             jtime = itime
             itime = ntime
          end if
       end do
       if( jtime == 0 ) then
          call runend('POR_INTERP: PROBLEM 1')
       else if( jtime == 1 ) then
          call runend('POR_INTERP: PROBLEM 2')
       else
          itim1 = jtime - 1
          itim2 = jtime 
       end if
    end if
    !
    ! Coordinate interpolation
    !
    if( .not. present(table_coord) ) then
       ipos1 = 1
       ipos2 = 1
    else
       if( coord <= table_coord(1) ) then
          ipos1 = 1
          ipos2 = 1
       else if( coord >= table_coord(nposi) ) then
          ipos1 = nposi
          ipos2 = nposi
       else
          iposi = 0
          jposi = 0
          do while( iposi < nposi )
             iposi = iposi + 1
             rcoor = table_coord(iposi)
             if( rcoor > coord ) then
                jposi = iposi
                iposi = nposi
             end if
          end do
          if( jposi == 0 ) then
             call runend('POR_INTERP: PROBLEM 1')
          else if( jposi == 1 ) then
             call runend('POR_INTERP: PROBLEM 2')
          else
             ipos1 = jposi - 1
             ipos2 = jposi 
          end if
       end if
       x1 = table_coord(ipos1)
       x2 = table_coord(ipos2)

    end if
    !     
    !    /|\
    !     |
    !  t2 o     +------+
    !     |     |   x  |
    !     |     |      |
    !  t1 o     +------+
    !     |
    !     +-----o------o----->
    !          x1     x2
    !
    t1 = table_value(itim1,1)
    t2 = table_value(itim2,1)
    !
    ! Interpolate v1 and v2
    !
    if( ipos1 == ipos2 ) then
       deltax = 1.0_rp
       x1     = coord
    else
       deltax = x2-x1
    end if
    if( itim1 == itim2 ) then
       deltat = 1.0_rp
    else
       deltat = t2-t1
    end if
    !
    ! Interpolate spatially at time t1 and t2
    !
    v_x1_t1 = table_value(itim1,ipos1+1)
    v_x1_t2 = table_value(itim2,ipos1+1)
    v_x2_t1 = table_value(itim1,ipos2+1)
    v_x2_t2 = table_value(itim2,ipos2+1)

    v_x_t1  = v_x1_t1 + (coord-x1)/deltax*(v_x2_t1-v_x1_t1)
    v_x_t2  = v_x1_t2 + (coord-x1)/deltax*(v_x2_t2-v_x1_t2)
    !
    ! Time interpolation
    !
    xvalu   = v_x_t1  + (cutim-t1)/deltat*(v_x_t2-v_x_t1)

  end subroutine por_interp

end module mod_por_interp
