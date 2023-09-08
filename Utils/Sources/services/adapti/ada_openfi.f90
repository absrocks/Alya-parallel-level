subroutine ada_openfi(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_openfi
  ! NAME
  !    par_openfi
  ! DESCRIPTION
  !    This subroutine gets ALL the file names to be used by Parall 
  !    service of Alya in two
  !    possible ways and them open them:
  !***
  !-----------------------------------------------------------------------

  use      def_kintyp
  use      def_parame
  use      def_master
  use      def_domain
  use      def_adapti
  use      mod_iofile
  
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iunit
  character(150)          :: fil_pdata_ada,fil_outpu_ada

  if(kfl_paral<=0) then

     select case(itask)

     case(1)
        if (kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR5901',fil_pdata_ada)
           call GET_ENVIRONMENT_VARIABLE('FOR5902',fil_outpu_ada)
        else if (kfl_naked==1) then
           fil_pdata_ada = adjustl(trim(namda))//'.'//exser(servi)//'.dat'
           fil_outpu_ada = adjustl(trim(namda))//'.'//exser(servi)//'.log'
           fil_oudom_ada = adjustl(trim(namda))//'-adapti.mgi.0-post.msh'
           fil_ougid_ada = adjustl(trim(namda))//'-adapti.post.msh'
           fil_ouerr_ada = adjustl(trim(namda))//'-adapti.post.res'
           fil_algeo_ada = adjustl(trim(namda))//'-a.dom.geo'
           fil_aldom_ada = adjustl(trim(namda))//'-a.dom.dat'
           fil_alfix_ada = adjustl(trim(namda))//'-a.fix'
           fil_alfix_ada = adjustl(trim(namda))//'-a.fix'
           fil_chpoi_ada = adjustl(trim(namda))//'-adapti.chk.out'
           fil_oumgi_ada = &
                adjustl(trim(namda))//'-adapti.mgi' ! this will be corrected in ada_outsau
        end if

        call iofile(zero,lun_pdata_ada,fil_pdata_ada,'ADAPTI DATA'    ,'old')
        call iofile(zero,lun_outpu_ada,fil_outpu_ada,'ADAPTI LOG FILE','unknown')

     case(2)
        !
        ! Output files which will be inputs for the mesh generator
        !
        call iofile(zero,lun_oumgi_ada,fil_oumgi_ada,'ADAPTI INPUT TO MG FILE','unknown')

     case(3)
        !
        ! Output files
        !
        call iofile(zero,lun_ougid_ada,fil_ougid_ada,'ADAPTI OUTGID FILE','unknown')
        call iofile(zero,lun_ouerr_ada,fil_ouerr_ada,'ADAPTI OUTERR FILE','unknown')
        call iofile(zero,lun_oudom_ada,fil_oudom_ada,'ADAPTI DOMAIN OUTPUT','unknown')
        call iofile(zero,lun_aldom_ada,fil_aldom_ada,'ADAPTI ALYA DOMDAT FILE','unknown')
        call iofile(zero,lun_algeo_ada,fil_algeo_ada,'ADAPTI ALYA DOMGEO FILE','unknown')
!        call iofile(zero,lun_alfix_ada,fil_alfix_ada,'ADAPTI ALYA FIXITY FILE','unknown')

        if (kfl_chpoi_ada > 0) &
             call iofile(&
             zero,lun_chpoi_ada,fil_chpoi_ada,'ADAPTI INTERPOLATED CHECKPOINT','unknown','unformatted')

     case(4)
        !
        ! Output file to append
        !
        call iofile(zero,lun_ougid_ada,fil_ougid_ada,'ADAPTI OUTGID FILE','old','formatted','append')

     case(10)
        !
        ! Close data file
        !
        call iofile(two,lun_pdata_ada,' ','DATA')
     case(11)
        !
        ! Close output files
        !
        call iofile(two,lun_outpu_ada,' ','LOG FILE')
        call iofile(two,lun_algeo_ada,' ','ALYA DOMGEO FILE')
        call iofile(two,lun_aldom_ada,' ','ALYA DOMDAT FILE')
!        call iofile(two,lun_alfix_ada,' ','ALYA FIXITY FILE')
        if (kfl_oudom_ada == 1) call iofile(two,lun_oudom_ada,' ','OUTDOM FILE')
        if (kfl_ougid_ada == 1) call iofile(two,lun_ougid_ada,' ','OUTGID FILE')
     end select

  end if

end subroutine ada_openfi
