subroutine ibm_boxwal(bobox,iwaib,ifind)

  use def_kintyp, only    :  ip,rp
  use def_master
  use def_domain
  use def_immbou
  use mod_kdtree
  implicit none

  integer(ip), intent(in)  :: iwaib
  real(rp),    intent(in)  :: bobox(3,2)
  integer(ip), intent(out) :: ifind

  integer(ip)              :: indst, curr, icond, idime
  real(rp)                 :: itemp,dismi1,dismi2


  
  ifind = 0_ip
  indst = 1_ip              
  twall_ibm(iwaib) % stru2(indst) = 1_ip                 
  do while( indst > 0_ip  )
     curr  = twall_ibm(iwaib) % stru2(indst)
     indst = indst - 1_ip
     !
     ! Minimum distance
     !       
     icond = 0_ip
     do idime=1,ndime
        if ( twall_ibm(iwaib) % sabox(1,idime,curr) < bobox(idime,2) .and. twall_ibm(iwaib) % sabox(2,idime,curr) > bobox(idime,1) ) then
           icond=icond+1_ip
        end if
     end do
     if (icond == ndime) then                 
        if( twall_ibm(iwaib) % blink(curr) < 0_ip ) then
           ifind = 1_ip                      
        else
           dismi1 = 0.0_rp
           dismi2 = 0.0_rp
           do idime = 1,ndime
              itemp  = max(0.0_rp,twall_ibm(iwaib) % sabox(1,idime,twall_ibm(iwaib) % blink(curr)) - bobox(idime,2))  
              itemp  = max(itemp,  bobox(idime,1) - twall_ibm(iwaib) % sabox(2,idime,twall_ibm(iwaib) % blink(curr))) 
              dismi1 = dismi1 + itemp * itemp
              
              itemp  = max(0.0_rp,twall_ibm(iwaib) % sabox(1,idime,twall_ibm(iwaib) % blink(curr)+1_ip) - bobox(idime,2))  
              itemp  = max(itemp,  bobox(idime,1) - twall_ibm(iwaib) % sabox(2,idime,twall_ibm(iwaib) % blink(curr)+1_ip))  
              dismi2 = dismi2 + itemp * itemp
           end do
           
           indst = indst + 2_ip                          
           if (dismi1 > dismi2) then
              twall_ibm(iwaib) % stru2(indst-1) = twall_ibm(iwaib) % blink(curr)
              twall_ibm(iwaib) % stru2(indst)   = twall_ibm(iwaib) % blink(curr) + 1_ip
           else
              twall_ibm(iwaib) % stru2(indst-1) = twall_ibm(iwaib) % blink(curr) + 1_ip
              twall_ibm(iwaib) % stru2(indst)   = twall_ibm(iwaib) % blink(curr)
           end if
        end if
     end if
  end do
  
end subroutine ibm_boxwal
