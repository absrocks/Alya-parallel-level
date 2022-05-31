subroutine lelpo2()
  !-----------------------------------------------------------------------
  !****f* domain/lelpo2
  ! NAME 
  !    lelpo2
  ! DESCRIPTION
  !    This subroutine computed extended graph
  ! USES
  !    chm_updunk
  ! USED BY
  !    endste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_graphs
  implicit none
  integer(ip) :: mepoi_2
  integer(ip) :: medge_2
  integer(ip) :: nedge_2

  if( INOTEMPTY ) then

     if( kfl_lele2 == 1 .and. kfl_lelp2 == 1 ) then
        !
        ! PELPO, LELPO
        ! LELEL, PELEL
        !
        call graphs_elepoi(&
             npoin_2,nelem_2,mnode,lnods,lnnod,mepoi_2,&
             pelpo_2,lelpo_2,&
             PELPO_NAME='PELPO_2',LELPO_NAME='LELPO_2')
        call graphs_eleele(&
             nelem_2,npoin_2,mnode,mepoi_2,lnods,lnnod,&
             pelpo_2,lelpo_2,nedge_2,medge_2,pelel_2,lelel_2)

     else if( kfl_lele2 == 1 .and. kfl_lelp2 == 0 ) then
        !
        ! PELEL, LELEL: deallocate memory for PELPO AND LELPO
        !
        call graphs_elepoi(&
             npoin_2,nelem_2,mnode,lnods,lnnod,mepoi_2,&
             pelpo_2,lelpo_2,&
             PELPO_NAME='PELPO_2',LELPO_NAME='LELPO_2')
        call graphs_eleele(&
             nelem_2,npoin_2,mnode,mepoi_2,lnods,lnnod,&
             pelpo_2,lelpo_2,nedge_2,medge_2,pelel_2,lelel_2)
        call graphs_elepoi_deallocate(&
             pelpo_2,lelpo_2,PELPO_NAME='PELPO_2',LELPO_NAME='LELPO_2')

     else if( kfl_lelp2 == 1 ) then
        !
        ! PELPO, LELPO
        !
        call graphs_elepoi(&
             npoin_2,nelem_2,mnode,lnods,lnnod,mepoi_2,&
             pelpo_2,lelpo_2,&
             PELPO_NAME='PELPO_2',LELPO_NAME='LELPO_2')

     end if


  end if

end subroutine lelpo2
