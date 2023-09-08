subroutine exm_outset()
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_outset
  ! NAME 
  !    exm_outset
  ! DESCRIPTION
  !    Compute and write results on sets:
  !    - Element, boundary and node sets:
  !      1. INTRA
  !      2. EXTRA
  !      3. RECOV
  ! USES
  ! USED BY
  !    exm_output
  !***
  !----------------------------------------------------------------------- 
  use def_parame
  use def_master
  use def_domain
  use def_exmedi
  use mod_iofile
  implicit none
  integer(ip) :: inset,ipoin,dummi

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  if( maxval(postp(1) % npp_setsn) > 0 ) then

     if( INOTMASTER ) then

        do inset = 1,nnset

           if( lnsec(inset) /= 0 ) then
              ipoin = lnsec(inset)
              if (kfl_cellmod(nodemat(ipoin)) == 1) then                  
                 if( postp(1) % npp_setsn(1) /= 0 ) &
                      vnset(1,inset) =  elmag(ipoin,1) * ( poref_fhn_exm(2) - poref_fhn_exm(1) ) + poref_fhn_exm(1)  
              else
                 if( postp(1) % npp_setsn(1) /= 0 ) &
                      vnset(1,inset) =  elmag(ipoin,1)
              end if
           end if 

        end do
     end if

     call posdef(23_ip,dummi)

  end if

end subroutine exm_outset
