subroutine rad_bcntoe
!-----------------------------------------------------------------------
!        
! This routine transforms the boundary conditions on nodes of 
! Neumann(Robin) type to Neumann(Robin) conditions on boundaries.
!
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_radiat
  use      mod_memchk
  implicit none
  integer(ip) :: ipoin,pnodb,iboun,inodb,iffix
  integer(ip) :: knodb(10)            ! Number of bc appearance
  real(rp)    :: vafix(10),value

  return

  do iboun=1,nboun

     ! Count number of appearance for each fixity on iboun
     ! e.g. knodb(1)= number of nodes with emmisivity
     !      knodb(2)= number of nodes with reflectivity applied
     pnodb=nnode(ltypb(iboun))
     knodb=0
     vafix=0.0_rp
     do inodb=1,pnodb
        ipoin=lnodb(inodb,iboun)
        iffix=kfl_fixno_rad(1,ipoin)
        if(iffix>0) then                           
           knodb(iffix)=knodb(iffix)+1
           vafix(iffix)=vafix(iffix)+bvess_rad(ipoin,1)
        end if
     end do
     where (knodb/=0) vafix=vafix/real(knodb)

     if(knodb(1)>=1) then
        !
        ! Emmisivity
        !
        kfl_fixbo_rad(iboun)=1
        do inodb=1,pnodb
           ipoin=lnodb(inodb,iboun)
           if(kfl_fixno_rad(1,ipoin)==1) then
              value=bvess_rad(ipoin,1)
           else
              value=vafix(1)
           end if
           bvnat_rad(1,iboun,1)=value
        end do
     end if

     if(knodb(2)>=1) then
        !
        ! Reflectivity
        !
        kfl_fixbo_rad(iboun)=2
        do inodb=1,pnodb
           ipoin=lnodb(inodb,iboun)
           if(kfl_fixno_rad(1,ipoin)==2) then
              value=bvess_rad(ipoin,1)
           else
              value=vafix(2)
           end if
           bvnat_rad(1,iboun,1)=1.0_rp-value ! epsilon = 1-rho
        end do
     end if

  end do

end subroutine rad_bcntoe
