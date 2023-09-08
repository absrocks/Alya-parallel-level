subroutine par_senset()
  !------------------------------------------------------------------------
  !****f* Parall/par_senset
  ! NAME
  !    par_senset
  ! DESCRIPTION
  !    Send/receive sets and witness points
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !------------------------------------------------------------------------
  use      def_kintyp
  use      def_parame
  use      def_parall
  use      def_domain
  use      def_master
  use      mod_memchk
  implicit none  
  integer(ip)          :: jpoin,ipoin,inset,domai,ii,kpoin
  integer(ip)          :: knset,jelem,jboun,iboun,ielem,indice0,jj
  integer(ip), pointer :: lnsec_loc(:)
  integer(ip), target  :: nnset_loc(1)

  if( IMASTER ) then
     
     if( neset > 0 ) then

        !----------------------------------------------------------------
        !
        ! Element sets
        !
        !----------------------------------------------------------------

        call memgen(1_ip,nelem,0_ip)
        do ielem = 1,nelem
           jelem        = leper_par(ielem)
           gisca(jelem) = leset(ielem)
        end do
        do kfl_desti_par = 1,npart_par
           npari =  nelem_par(kfl_desti_par)
           parin => gisca(leind_par(kfl_desti_par):)
           strin =  'GISCA'
           call par_sendin()
        end do
        call memgen(3_ip,nelem,0_ip)
        call memose(-1_ip)
        
     end if

     if( nbset > 0 ) then

        !----------------------------------------------------------------
        !
        ! Boundary sets
        !
        !----------------------------------------------------------------

        call memgen(1_ip,nboun,0_ip)
        do iboun = 1,nboun
           jboun        = lbper_par(iboun)
           gisca(jboun) = lbset(iboun)
        end do 
        do kfl_desti_par = 1,npart_par
           npari =  nboun_par(kfl_desti_par)
           parin => gisca(lbind_par(kfl_desti_par):)
           strin =  'LBSET'
           call par_sendin()
        end do 
        call memgen(3_ip,nboun,0_ip)
        call memose(-2_ip)

     end if

     if( nnset > 0 ) then
        
        !----------------------------------------------------------------
        !
        ! Node sets
        !
        !----------------------------------------------------------------

        call memgen(1_ip,npoin_total,0_ip)
        do ipoin = 1,npoin_total
           jpoin        = lninv_loc(ipoin)
           gisca(ipoin) = lnset(jpoin)
        end do
        indice0 = 1
        do kfl_desti_par = 1,npart_par
           npari =  npoin_par(kfl_desti_par)
           parin => gisca(indice0:)
           strin =  'LNSET'
           call par_sendin()
           indice0 = indice0 + npoin_par(kfl_desti_par)
        end do
        call memgen(3_ip,npoin,0_ip)
        call memose(-12_ip)

     end if
     
  else if( ISLAVE ) then

     kfl_desti_par = 0

     if( neset > 0 ) then
        !
        ! Element sets
        !
        call memose( 1_ip)
        npari =  nelem
        parin => leset(:)
        call par_receiv()
     end if

     if( nbset > 0 ) then
        !
        ! Boundary sets
        !
        call memose( 2_ip)
        npari =  nboun
        parin => lbset(:)
        call par_receiv()
     end if

     if( nnset > 0 ) then
        !
        ! Node sets 
        !
        call memose(12_ip)
        npari =  npoin
        parin => lnset(:)
        call par_receiv()
     end if

  end if

end subroutine par_senset
