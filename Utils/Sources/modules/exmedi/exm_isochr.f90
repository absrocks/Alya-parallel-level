subroutine exm_isochr(icomp)
  !-----------------------------------------------------------------------
  !****f* exmedi/exm_isochr
  ! NAME 
  !    exm_isochr
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    exm_isochr
  !***
  !-----------------------------------------------------------------------

  use      def_parame
  use      def_master
  use      def_domain
  use      def_postpr
  use      def_exmedi
  use      mod_communications, only : PAR_MAX

  implicit none
  integer(ip), intent(in)  :: icomp
  integer(ip)              :: ipoin
  real(rp)                 :: xauxi
  logical                  :: auto_isochrones
  integer(1)               :: upstroke_or_downstroke
  integer(4)               :: save_isochrones, istat
  real(rp), allocatable    :: fisoc_temp(:)

  if (thiso_exm(1)>0.0_rp) then    ! isochrones trigger         


          save_isochrones = 0
          auto_isochrones = thiso_exm(3)>0.0_rp

          !the deafult behaviour is that it saves only the first upstroke
          !if auto_isochrones are set, both upstrokes and downstrokes are recorded. 
          !The file of isochrones is saved when the same point gets two status changes (upstroke, then downstroke).
          !In the latter case, the state of the isochrones is saved before applying the change.
          !I.e. if point 100 has an upstroke recorded and then downstroke happens, 
          !the file is saved and then the downstroke is recorded.
          if( .NOT. auto_isochrones ) then
              if( INOTMASTER) then

                do ipoin= 1,npoin

                    !!        if (kfl_cellmod(nodemat(ipoin)) == CELL_FITZHUGH_EXMEDI) then 
                    !           xauxi = elmag(ipoin,icomp) * ( poref_fhn_exm(2) - poref_fhn_exm(1) ) + poref_fhn_exm(1)        

                    xauxi = elmag(ipoin,icomp) 
                   
                    if (kwave_exm(ipoin) == 0) then !expecting an upstroke
                       if (xauxi > thiso_exm(2)) then ! isochrones threshold
                          fisoc(ipoin) = cutim
                          kwave_exm(ipoin) = 1
                       end if
                    end if
                    
                    !!        end if
                 end do !ipoin

              end if !INOTMASTER

          else !.NOT. auto_isochrones

             if( INOTMASTER) then
                 allocate( fisoc_temp(npoin), STAT=istat )
                 if (istat /= 0) then
                    call runend("exm_isochr.f90: failed to allocate fisoc_temp(npoin)")
                 end if

                 fisoc_temp = fisoc

                 do ipoin= 1,npoin

                    !!        if (kfl_cellmod(nodemat(ipoin)) == CELL_FITZHUGH_EXMEDI) then 
                    !           xauxi = elmag(ipoin,icomp) * ( poref_fhn_exm(2) - poref_fhn_exm(1) ) + poref_fhn_exm(1)        

                    xauxi = elmag(ipoin,icomp) 

                    upstroke_or_downstroke = 0
                
                    
                    
                    !save upstrokes with positive time and downstrokes with negative time. Every other point will have 0, not -1            
                    if (kwave_exm(ipoin) == 0) then !expecting an upstroke
                       if (xauxi > thiso_exm(2)) then ! do we have an upstroke?
                          kwave_exm(ipoin) = 1
                          upstroke_or_downstroke = 1 
                          !print *, 'Upstroke detected ipoin: ', ipoin, ' xauxi = ',xauxi
                       end if
                    else if (kwave_exm(ipoin) == 1) then !expecint a downstroke
                       if (xauxi < thiso_exm(2)) then ! do we have a downstroke?
                          kwave_exm(ipoin) = 0
                          upstroke_or_downstroke = -1
                          !print *, 'Downstroke detected ipoin: ', ipoin, ' xauxi = ',xauxi
                       end if
                    end if


                    if (upstroke_or_downstroke /=0 ) then

                        if( .NOT. isoch_modified(ipoin) ) then
                            isoch_modified(ipoin) = .TRUE. !modified
                            fisoc_temp(ipoin) = cutim*upstroke_or_downstroke
                        else
                            !save isochrones
                            !print *,"Saving isochrones"
                            !kfl_ivari(1) = 26
                            !call exm_outvar(26)
                            save_isochrones = 1
                            
                            isoch_modified = .FALSE. !reset modified flag
                            isoch_modified(ipoin) = .TRUE. !set modified just for this point
                            fisoc_temp(ipoin) = cutim*upstroke_or_downstroke
                        end if

                    end if !upstroke_or_downstroke

                 end do !ipoin

             end if ! INOTMASTER
          

             if (auto_isochrones .AND. IPARALL) then
                call PAR_MAX(save_isochrones,'IN MY CODE')     
             end if


             if (auto_isochrones .AND. save_isochrones>0) then
                kfl_ivari(1) = 26_ip !ISOCH
                call exm_outvar(26_ip)
             end if
         
             if( auto_isochrones .AND. INOTMASTER ) then
                 fisoc = fisoc_temp
                 deallocate(fisoc_temp)
             end if

          end if !.NOT. auto_isochrones

  end if  ! isochrones trigger   


end subroutine exm_isochr
