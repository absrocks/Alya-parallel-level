subroutine ibm_openfi(itask)
  !------------------------------------------------------------------------
  !****f* master/ibm_openfi
  ! NAME 
  !    ibm_openfi
  ! DESCRIPTION
  !    This subroutine 
  ! OUTPUT
  ! USES
  ! USED BY
  !    Turnon
  !    Turnof
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_postpr
  use def_immbou
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask
  character(150)          :: fil_resib_tmp
  character(150)          :: fil_outpu_ibm

  if( INOTSLAVE ) then

     select case ( itask )
        
     case ( 0_ip )

        !----------------------------------------------------------------
        !
        ! Open IB result file          
        !        
        !----------------------------------------------------------------

        if( nimbo > 0 ) then
           if (kfl_naked==0) then
              call GET_ENVIRONMENT_VARIABLE('FOR2114',fil_outpu_ibm)               
           else if (kfl_naked==1) then
              fil_outpu_ibm = adjustl(trim(namda))//'-IB.res'
           end if
           call iofile(zero,lun_outpu_ibm,fil_outpu_ibm,'IB RESULTS')
        end if

     case ( 1_ip )

        !----------------------------------------------------------------
        !
        ! Open IB mesh file
        !     
        !----------------------------------------------------------------

        if (kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR2110',fil_mshib)
           call GET_ENVIRONMENT_VARIABLE('FOR2111',fil_resib)
           call GET_ENVIRONMENT_VARIABLE('FOR2112',fil_mshi2)
           call GET_ENVIRONMENT_VARIABLE('FOR2113',fil_resi2)
        else
           fil_mshib = adjustl(trim(namda))
           fil_resib = adjustl(trim(namda))
           fil_mshi2 = adjustl(trim(namda))
           fil_resi2 = adjustl(trim(namda))
        end if

        if( kfl_outfo == 1 ) then

           fil_mshib = adjustl(trim(fil_mshib))//'-IB.post.msh'
           call iofile(zero,lun_mshib,fil_mshib,'IB MESH')

        else if( kfl_outfo == 3000 ) then

           fil_mshib = adjustl(trim(fil_mshib))//'-IB.post-bin.msh'
           call iofile(zero,lun_mshib,fil_mshib,'IB BIN MESH','unknown','unformatted')

        else if( kfl_outfo == 30 ) then

           fil_mshib = adjustl(trim(fil_mshib))//'-IB.msh.vu'
           call iofile(zero,lun_mshib,fil_mshib,'IB VU MESH')
           fil_mshi2 = adjustl(trim(fil_mshi2))//'-IB.msh.vubin'
           call iofile(zero,lun_mshi2,fil_mshi2,'IB VU BIN MESH','unknown','unformatted')

        end if

       if( kfl_outfo == 1 ) then

           fil_resib = adjustl(trim(fil_resib))//'-IB.post.res'

        else if( kfl_outfo == 3000 ) then

           fil_resib = adjustl(trim(fil_resib))//'-IB.post-bin.res'

        else if( kfl_outfo == 30 ) then

           fil_resib = adjustl(trim(fil_resib))//'-IB.res'
           fil_resi2 = adjustl(trim(fil_resi2))//'-IB.res'

        end if

     case ( 2_ip )
 
        !----------------------------------------------------------------
        !
        ! Open IB postprocess
        !
        !----------------------------------------------------------------

        if( kfl_outfo == 1 ) then

           call iofile(zero,lun_resib,fil_resib,'IB POSTPROCESS')

        else if( kfl_outfo == 3000 ) then

           fil_resib_tmp = fil_resib
           call appnum(ittim,fil_resib_tmp)
           call iofile(zero,lun_resib,fil_resib_tmp,'IB BIN POSTPROCESS','unknown','unformatted')

        else if( kfl_outfo == 30 ) then

           fil_resib_tmp = adjustl(trim(fil_resib))
           call appnum(ittim,fil_resib_tmp)
           fil_resib_tmp = trim(fil_resib_tmp)//'.vu'
           call iofile(zero,lun_resib,fil_resib_tmp,'IB VU POSTPROCESS')

           fil_resi2 = adjustl(trim(fil_resib))
           call appnum(ittim,fil_resi2)
           fil_resi2 = trim(fil_resi2)//'.vubin'
           call iofile(zero,lun_resi2,fil_resi2,'IB VU BIN POSTPROCESS','unknown','unformatted')

        end if

     case ( 3_ip )

        !----------------------------------------------------------------
        !
        ! Close IB mesh
        !
        !----------------------------------------------------------------

        if( kfl_outfo == 1 .or. kfl_outfo == 3000 .or. kfl_outfo == 30 ) then
           call iofile(two,lun_mshib,' ','IB MESH')
        end if
        if( kfl_outfo == 30 ) then
           continue ! files are closed by geomvu
        end if

     case ( 4_ip )
 
        !----------------------------------------------------------------
        !
        ! Close IB postprocess
        !
        !----------------------------------------------------------------

        if( kfl_outfo == 1 .or. kfl_outfo == 3000 ) then
           call iofile(two,lun_resib,' ','IB POSTPROCESS') 
        else if( kfl_outfo == 30 ) then
           call iofile(two,lun_resib,' ','IB VU POSTPROCESS') 
           call iofile(two,lun_resi2,' ','IB VU BIN POSTPROCESS') 
        end if


     end select

  end if


end subroutine ibm_openfi
