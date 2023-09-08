module mod_opfpos

    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    use def_inpout
    use mod_iofile
    use def_postpr
    use def_mpio

    implicit none

    public                                     :: opfpos_name, opfposvx, postpr_intto8

    contains

    subroutine opfpos_name(name,extension,mesh,TAG1,TAG2)
        use def_mpio, only                     :  mpio_ext
        implicit none
        character(5)  ,  intent(in)            :: name
        integer(ip)   ,  intent(in), optional  :: TAG1
        integer(ip)   ,  intent(in), optional  :: TAG2
        character(150)                         :: filsa
        character(*),    intent(in)            :: extension
        logical       ,  intent(out)           :: mesh
        character(150)                         :: wtags
        character(20)                          :: wtag1
        character(20)                          :: wtag2

        mesh=.false.


         !-------------------------------------------------------------------
         !
         ! Open postprocess file
         !
         !-------------------------------------------------------------------
        !
        ! Time step identifier
        !
        nunam_pos=postpr_intto8(ittim)

        !
        ! Compose file name: example.VELOC-123456.h5
        !
        filsa = adjustl(trim(namda))
        !
        ! File name
        !
        filsa = trim(filsa)//'-'//trim(name)
        wtags = ''
        if( present(TAG1) ) then
            if (TAG1 /= -1) then
                wtags = trim(wtags)//'.'//postpr_intto8(TAG1)
            end if
        end if
        if( present(TAG2) ) then
            if (TAG2 /= -1) then
                wtags = trim(wtags)//'.'//postpr_intto8(TAG2)
            end if
        end if
        filsa = trim(filsa)//trim(wtags)

        if(name == 'COORD' .or. &
           name == 'LNODS' .or. &
           name == 'LTYPE' .or. &
           name == 'LNINV' .or. &
           name == 'LELCH' .or. &
           name == 'LNODB' .or. &
           name == 'LTYPB' .or. &
           name == 'LBOCH' .or. &
           name == 'LNOCH' .or. &
           name == 'LBINV' .or. &
           name == 'LELBO' .or. &
           name == 'LESUB' .or. &
           name == 'LMATE' .or. &
           name == 'LMAST' .or. &
           name == 'LEINV' .or. &
           name == 'CODBO' .or. &
           name == 'CODNO' .or. &
           name == 'LESET' .or. &
           name == 'LBSET' .or. &
           name == 'XFIEL' ) then
             !
             ! Open geometry arrays file
             !
           fil_postp = trim(filsa)//trim(post_ext)//trim(extension)
        if (kfl_mpio_export==1) then
            fil_postp = trim(filsa)//trim(extension)
        end if
           mesh=.true.

        else

            if( kfl_reawr == 0 ) then
            !
            ! Open postprocess file name
            !
                fil_postp = trim(filsa)//'-'//adjustl(trim(nunam_pos))//trim(post_ext)//trim(extension)

            else
                fil_postp = trim(filsa)//'.rst'!//trim(extension)
                if (extension==mpio_ext) then
                    fil_postp = trim(fil_postp)//trim(extension)
                end if
                if( kfl_reawr == 2) then
                !
                ! Open restart file name for writing
                !
                    if( kfl_rsfil == 1 ) call appnum(ittim,fil_postp)
                end if
            end if
        end if

    end subroutine

    subroutine opfposvx()
      !-----------------------------------------------------------------------
      !****f* Domain/opfpos
      ! NAME
      !    opfpos
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

      implicit none
      character(150) :: filsa
      character(200) :: messa

      !-------------------------------------------------------------------
      !
      ! Open voxel file
      !
      !-------------------------------------------------------------------

      if( INOTSLAVE ) then
         nunam_pos=postpr_intto8(ittim)
         filsa = adjustl(trim(namda))
         fil_posvx = trim(filsa)//'-'//trim(wopos_pos(1))//'-'//adjustl(trim(nunam_pos))//'.bvox'
         open (unit=lun_posvx, FILE = trim(fil_posvx) ,ACTION='WRITE', &
              & STATUS = 'UNKNOWN', form='unformatted', access='stream' )
      end if

    end subroutine


  function postpr_intto8(v) result(v_out)

    integer(ip),  intent(in) :: v
    character(8)             :: v_out

    write(v_out,'(a)') 'XXXXXXXX'
    if(      v < 10      ) then
       write(v_out,'(a,i1)') '0000000',v
    else if( v < 100      ) then
       write(v_out,'(a,i2)') '000000',v
    else if( v < 1000     ) then
       write(v_out,'(a,i3)') '00000',v
    else if( v < 10000    ) then
       write(v_out,'(a,i4)') '0000',v
    else if( v < 100000   ) then
       write(v_out,'(a,i5)') '000',v
    else if( v < 1000000  ) then
       write(v_out,'(a,i6)') '00',v
    else if( v < 10000000 ) then
       write(v_out,'(a,i7)') '0',v
    else if( v < 100000000 ) then
       write(v_out,'(i8)') v
    end if

  end function

end module
