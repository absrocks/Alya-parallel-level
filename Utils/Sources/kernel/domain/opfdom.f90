subroutine opfdom(itask)
  !-----------------------------------------------------------------------
  !****f* Domain/opfdom
  ! NAME
  !    opfdom
  ! DESCRIPTION
  !    This subroutine gets ALL the file names to be used by Alya in two
  !    possible ways and them open them:
  ! 
  !    1. Recalling them from the environment, when Alya is launched
  !    encapsulated in a shell script, or
  ! 
  !    2. Composing the names out of the problem name which is given as argument
  !    when the binary file Alya is launched "naked". 
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_inpout
  use mod_iofile
  use def_postpr
  implicit none
  integer(ip), intent(in) :: itask
  character(150)          :: fil_pdata_dom,fil_elsta_dom
  character(150)          :: fil_elmsh_dom,fil_elres_dom

  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if( INOTSLAVE ) then

     select case (itask)

     case(1_ip)
        !
        ! Get file names:
        !
        ! If kfl_naked=0 -->> encapsulated, then get names from the environment (DEFAULT value)
        ! If kfl_naked=1 -->> naked, then compose the names
        !    
        if (kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR021',fil_pdata_dom) 
           call GET_ENVIRONMENT_VARIABLE('FOR025',fil_elsta_dom)
           call GET_ENVIRONMENT_VARIABLE('FOR026',fil_elmsh_dom)
           call GET_ENVIRONMENT_VARIABLE('FOR027',fil_elres_dom)

        else if (kfl_naked==1) then
           fil_pdata_dom = adjustl(trim(namda))//'.dom.dat'
           fil_elsta_dom = adjustl(trim(namda))//'.els.log'
           fil_elmsh_dom = adjustl(trim(namda))//'.els.post.msh'
           fil_elres_dom = adjustl(trim(namda))//'.els.post.res'
        end if
        !
        ! Open file
        !
        call iofile(zero,lun_pdata_dom,fil_pdata_dom,'DOMAIN','old')
        !
        ! Open postprocess files
        !
        if (kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR030',fil_pos00) 
           call GET_ENVIRONMENT_VARIABLE('FOR031',fil_pos01)
           call GET_ENVIRONMENT_VARIABLE('FOR032',fil_pos02)           
        else if (kfl_naked==1) then
           fil_pos00 = adjustl(trim(namda))//'.post.alyapar'
           fil_pos01 = adjustl(trim(namda))//'.post.alyafil'
           fil_pos02 = adjustl(trim(namda))//'.post.alyalog'           
        end if
        call iofile(zero,lun_pos00,fil_pos00,'POSTPROCESS PARALL INFO')
        !
        ! Define unit opening option if this is a restart run for alyafil
        !
        if(kfl_rstar==2) then 
           statu='old'
           forma='formatted'
           posit='append'
        else
           statu='unknown'
           forma='formatted'
           posit='asis'
        end if
        call iofile(zero,lun_pos01,fil_pos01,'POSTPROCESS FILE NAMES',statu,forma,posit)
        call iofile(zero,lun_pos02,fil_pos02,'POSTPROCESS LOG FILE')

     case(2_ip)

     case(3_ip)
        !
        ! Close domain data file
        !
        call iofile(two,lun_pdata_dom,'','DOMAIN')
        lispa = 0
        lisda = lun_pdata      ! Recover data file
        lisre = lun_outpu      ! Recover results file

     case(4_ip)

     case(5_ip)
        !
        ! Open Elsest files
        !
        if(ielse( 7)/=0) call iofile(zero,lun_elsta_dom,fil_elsta_dom,'ELSEST STATISTICS')
        if(ielse(12)/=0) call iofile(zero,lun_elmsh_dom,fil_elmsh_dom,'ELSEST MESH')
        if(ielse(13)/=0) call iofile(zero,lun_elres_dom,fil_elres_dom,'ELSEST RESULTS')

     end select

  end if

end subroutine opfdom
