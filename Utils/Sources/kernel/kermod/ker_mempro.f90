subroutine ker_mempro(itask)
  use def_master
  use def_parame
  use def_kermod
  use def_domain
  use mod_ker_proper  
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( 1_ip )
     !
     ! Nullify laws
     !
     nullify( densi_ker % wlaws )
     nullify( densi_ker % ilaws )
     nullify( densi_ker % rlaws )
     nullify( densi_ker % value_default )
     nullify( densi_ker % update )

     nullify( visco_ker % wlaws )
     nullify( visco_ker % ilaws )
     nullify( visco_ker % rlaws )
     nullify( visco_ker % value_default )
     nullify( visco_ker % update )

     nullify( poros_ker % wlaws )
     nullify( poros_ker % ilaws )
     nullify( poros_ker % rlaws )
     nullify( poros_ker % value_default )
     nullify( poros_ker % update )

     nullify( condu_ker % wlaws )
     nullify( condu_ker % ilaws )
     nullify( condu_ker % rlaws )
     nullify( condu_ker % value_default )
     nullify( condu_ker % update )

     nullify( sphea_ker % wlaws )
     nullify( sphea_ker % ilaws )
     nullify( sphea_ker % rlaws )
     nullify( sphea_ker % value_default )
     nullify( sphea_ker % update)

     nullify( dummy_ker % wlaws )
     nullify( dummy_ker % ilaws )
     nullify( dummy_ker % rlaws )
     nullify( dummy_ker % value_default )
     nullify( dummy_ker % update )

     nullify( turmu_ker % wlaws )
     nullify( turmu_ker % ilaws )
     nullify( turmu_ker % rlaws )
     nullify( turmu_ker % value_default )
     nullify( turmu_ker % update )
     
     nullify( absor_ker % wlaws )
     nullify( absor_ker % ilaws )
     nullify( absor_ker % rlaws )
     nullify( absor_ker % value_default )
     nullify( absor_ker % update )
     
     nullify( scatt_ker % wlaws )
     nullify( scatt_ker % ilaws )
     nullify( scatt_ker % rlaws )
     nullify( scatt_ker % value_default )
     nullify( scatt_ker % update )

     nullify( mixin_ker % wlaws )
     nullify( mixin_ker % ilaws )
     nullify( mixin_ker % rlaws )
     nullify( mixin_ker % value_default )
     nullify( mixin_ker % update )      
     !
     !
     ! Allocate laws
     !
     allocate( densi_ker % wlaws(nmate) )
     allocate( densi_ker % ilaws(nmate) )
     allocate( densi_ker % rlaws(mlapa_ker,nmate) )
     allocate( densi_ker % value_default(nmate) )
     allocate( densi_ker % update(2,nmate) )

     allocate( visco_ker % wlaws(nmate) )
     allocate( visco_ker % ilaws(nmate) )
     allocate( visco_ker % rlaws(mlapa_ker,nmate) )
     allocate( visco_ker % value_default(nmate) )
     allocate( visco_ker % update(2,nmate) )

     allocate( poros_ker % wlaws(nmate) )
     allocate( poros_ker % ilaws(nmate) )
     allocate( poros_ker % rlaws(mlapa_ker,nmate) )
     allocate( poros_ker % value_default(nmate) )
     allocate( poros_ker % update(2,nmate) )

     allocate( condu_ker % wlaws(nmate) )
     allocate( condu_ker % ilaws(nmate) )
     allocate( condu_ker % rlaws(mlapa_ker,nmate) )
     allocate( condu_ker % value_default(nmate) )
     allocate( condu_ker % update(2,nmate) )

     allocate( sphea_ker % wlaws(nmate) )
     allocate( sphea_ker % ilaws(nmate) )
     allocate( sphea_ker % rlaws(mlapa_ker,nmate) )
     allocate( sphea_ker % value_default(nmate) )
     allocate( sphea_ker % update(2,nmate) )

     allocate( dummy_ker % wlaws(nmate) )
     allocate( dummy_ker % ilaws(nmate) )
     allocate( dummy_ker % rlaws(mlapa_ker,nmate) )
     allocate( dummy_ker % value_default(nmate) )
     allocate( dummy_ker % update(2,nmate) )

     allocate( turmu_ker % wlaws(nmate) )
     allocate( turmu_ker % ilaws(nmate) )
     allocate( turmu_ker % rlaws(mlapa_ker,nmate) )
     allocate( turmu_ker % value_default(nmate) )
     allocate( turmu_ker % update(2,nmate) )
     
     allocate( absor_ker % wlaws(nmate) )
     allocate( absor_ker % ilaws(nmate) )
     allocate( absor_ker % rlaws(mlapa_ker,nmate) )
     allocate( absor_ker % value_default(nmate) )
     allocate( absor_ker % update(2,nmate) )
     
     allocate( scatt_ker % wlaws(nmate) )
     allocate( scatt_ker % ilaws(nmate) )
     allocate( scatt_ker % rlaws(mlapa_ker,nmate) )
     allocate( scatt_ker % value_default(nmate) )
     allocate( scatt_ker % update(2,nmate) )

     allocate( mixin_ker % wlaws(nmate) )
     allocate( mixin_ker % ilaws(nmate) )
     allocate( mixin_ker % rlaws(mlapa_ker,nmate) )
     allocate( mixin_ker % value_default(nmate) )
     allocate( mixin_ker % update(2,nmate) )
     !
     ! Put default values to negative value to identfy if the default
     ! exists or not
     ! 
     densi_ker % value_default(1:nmate) = -1.0_rp
     visco_ker % value_default(1:nmate) = -1.0_rp
     poros_ker % value_default(1:nmate) = -1.0_rp
     condu_ker % value_default(1:nmate) = -1.0_rp
     sphea_ker % value_default(1:nmate) = -1.0_rp
     dummy_ker % value_default(1:nmate) = -1.0_rp
     turmu_ker % value_default(1:nmate) = -1.0_rp
     absor_ker % value_default(1:nmate) = -1.0_rp
     scatt_ker % value_default(1:nmate) = -1.0_rp
     mixin_ker % value_default(1:nmate) = -1.0_rp
     !
     ! initialization of where to update to default values 
     ! 
     densi_ker % update(1, 1:nmate) = ITASK_DOITER
     visco_ker % update(1, 1:nmate) = ITASK_DOITER
     poros_ker % update(1, 1:nmate) = ITASK_DOITER
     condu_ker % update(1, 1:nmate) = ITASK_DOITER
     sphea_ker % update(1, 1:nmate) = ITASK_DOITER
     dummy_ker % update(1, 1:nmate) = ITASK_DOITER
     turmu_ker % update(1, 1:nmate) = ITASK_DOITER

     densi_ker % update(2, 1:nmate) = ITASK_ENDSTE
     visco_ker % update(2, 1:nmate) = ITASK_ENDSTE
     poros_ker % update(2, 1:nmate) = ITASK_ENDSTE
     condu_ker % update(2, 1:nmate) = ITASK_ENDSTE
     sphea_ker % update(2, 1:nmate) = ITASK_ENDSTE
     dummy_ker % update(2, 1:nmate) = ITASK_ENDSTE
     turmu_ker % update(2, 1:nmate) = ITASK_ENDSTE
     absor_ker % update(2, 1:nmate) = ITASK_ENDSTE
     scatt_ker % update(2, 1:nmate) = ITASK_ENDSTE
     mixin_ker % update(2, 1:nmate) = ITASK_ENDSTE

  case ( 2_ip ) 

     if( INOTMASTER ) then
        !
        ! Law name vs law number
        !
        if( densi_ker % kfl_exist == 1 ) call ker_wiprop(densi_ker)
        if( visco_ker % kfl_exist == 1 ) call ker_wiprop(visco_ker)
        if( poros_ker % kfl_exist == 1 ) call ker_wiprop(poros_ker)
        if( condu_ker % kfl_exist == 1 ) call ker_wiprop(condu_ker)
        if( sphea_ker % kfl_exist == 1 ) call ker_wiprop(sphea_ker)
        if( dummy_ker % kfl_exist == 1 ) call ker_wiprop(dummy_ker)
        if( turmu_ker % kfl_exist == 1 ) call ker_wiprop(turmu_ker)
        if( absor_ker % kfl_exist == 1 ) call ker_wiprop(absor_ker)
        if( scatt_ker % kfl_exist == 1 ) call ker_wiprop(scatt_ker)
        if( mixin_ker % kfl_exist == 1 ) call ker_wiprop(mixin_ker)
        !
        ! Allocate memory for properties according to law 
        !
        if( densi_ker % kfl_exist == 1 ) call ker_allpro(densi_ker)
        if( visco_ker % kfl_exist == 1 ) call ker_allpro(visco_ker)
        if( poros_ker % kfl_exist == 1 ) call ker_allpro(poros_ker)
        if( condu_ker % kfl_exist == 1 ) call ker_allpro(condu_ker)
        if( sphea_ker % kfl_exist == 1 ) call ker_allpro(sphea_ker)
        if( dummy_ker % kfl_exist == 1 ) call ker_allpro(dummy_ker)
        if( turmu_ker % kfl_exist == 1 ) call ker_allpro(turmu_ker)
        if( absor_ker % kfl_exist == 1 ) call ker_allpro(absor_ker)
        if( scatt_ker % kfl_exist == 1 ) call ker_allpro(scatt_ker)
        if( mixin_ker % kfl_exist == 1 ) call ker_allpro(mixin_ker)

     end if
     
  end select

end subroutine ker_mempro
 
