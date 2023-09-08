subroutine qua_begite
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_begite
  ! NAME 
  !    qua_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the Schrodinger
  !    equation
  ! USES
  !    qua_tistep
  !    qua_updbcs
  !    qua_inisol
  !    qua_updunk
  ! USED BY
  !    qua_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_quanty
  use mod_memchk
  use mod_outfor, only : outfor
  use mod_messages, only : livinf
  implicit none

  character*20    :: fname
  integer(ip)     :: kato,kk,kespe,cont_espe,nespec,kocup,nspin_tot
  integer(4)      :: istat
  !
  ! Initializations
  !
  print*,'Begite: start=',kfl_paral,kfl_dftgs_qua,kfl_alele_qua
  kfl_goite_qua = 1 
  itinn(modul)  = 0
  if(itcou==1) call qua_tistep()
  call livinf(15_ip,' ',modul)
  !
  ! Set up the solver parameters for the schrodinger equation
  !
  eigen_sol => eigen_qua ! Eigen Solver type
  solve_sol => solve_qua ! Algeb Solver type

  if( kfl_dftgs_qua == 1 .or. kfl_alele_qua == 1 ) then
     !
     !  ahora, con toda la info, determino ocupaciones, estados a calcular, etc
     !  numero de estados a calcular --> nestates
     !
     nestates = 0
     nspin_tot= 0 
     do kato = 1,atomo_qua(1) % natoms
        nestates = nestates + especie_qua(atomo_qua(1)%Tiespe(kato))%valencia
        nspin_tot= nspin_tot+ atomo_qua(1)%spin(kato)
     end do
     !
     !  numero de ocupacion          --> noccupa(natoms)
     !
     if( mod(nestates,2) == 0) then
        nestates=nestates/2
        ! ocupacion!
        allocate(noccupa(nestates),stat=istat)
        call memchk(zero,istat,memdi,'noccu','qua_3despa',noccupa)
        do kocup = 1,nestates
           noccupa(kocup) = 2
        end do

     else

        write(lun_outpu_qua,*) 'ODD Number of states. Warning! alya will make unpolarized calculi'  
        nestates = (nestates+1)/2
        ! ocupacion!
        allocate(noccupa(nestates),stat=istat)
        call memchk(zero,istat,memdi,'noccu','qua_3despa',noccupa)
        do kocup = 1,nestates-1
           noccupa(kocup)=2
        end do
        noccupa(nestates)=1

     endif

     !
     ! now I will set the eigensolver. Number de estates to be calculated.
     !
     eigen_qua(1)%neiva = nestates 
     eigen_sol(1)%neiva = nestates 


     if( INOTMASTER ) then

        do kk = 1,npoin
           rhoon(kk,1) = 0.0_rp
        end do
        !
        !  alloco la funcion potencial del pseudo que luego uso en la eq de Schrodinger
        !
        allocate(v_pot_ps(npoin),stat=istat)          
        call memchk(zero,istat,memdi,'vpotp','qua_begite',v_pot_ps)
        !
        ! alloco tambien de paso los potenciales de hartree y XC
        !
        allocate(v_hartree(npoin),stat=istat)          
        call memchk(zero,istat,memdi,'vhart','qua_begite',v_hartree)
        allocate(v_xc(npoin),stat=istat)          
        call memchk(zero,istat,memdi,'vpxc ','qua_begite',v_xc)
        !
        !  realloco vectores para el Eigen solver
        !
        deallocate(eigen)
        deallocate(eigva)
        allocate(eigen(nestates*npoin),stat=istat)
        call memchk(zero,istat,memdi,'EIGEN','qua_begite',eigen)
        allocate(eigva(nestates),stat=istat)
        call memchk(zero,istat,memdi,'EIGVA','qua_begite',eigva)

        allocate(NL_denom(lmaximo+1,2*lmaximo+1),stat=istat)
        call memchk(zero,istat,memdi,'NL_denom','qua_begite',NL_denom)
        allocate(NonLoc(lmaximo+1,2*lmaximo+1,npoin),stat=istat)
        call memchk(zero,istat,memdi,'NonLoc','qua_begite',NonLoc)

        Nonloc   = 0
        Nl_denom = 0
        do kk = 1,npoin
           v_pot_ps(kk)  = 0.0_rp
           v_hartree(kk) = 0.0_rp
           v_xc(kk)      = 0.0_rp
        end do
        !
        !  armo las funciones de onda localizadas en cada atomo. Solo para verificacion y graficos
        !
        print*,'Begite: all electron=',kfl_paral

        if( kfl_alele_qua == 1 ) then
           !
           ! si es all electron contruyo la densidad total + saco armonicos esfericos 
           ! y boludeces varias para chequear
           !
           do kato = 1,atomo_qua(1)%natoms
              call all_electron(kato,atomo_qua(1)%spin(kato),atomo_qua(1)%espe(kato), & 
                   atomo_qua(1)%coorx(kato),atomo_qua(1)%coory(kato),atomo_qua(1)%coorz(kato),atomo_qua(1)%natoms)
           end do


        else
           !
           ! en este caso debo armar el pseudopotencial general expandido sobre todo el dominio 
           !
           do kato = 1,atomo_qua(1)%natoms
              nespec = atomo_qua(1)%Tiespe(kato)
              call arma_pp_3D(kato,nespec,atomo_qua(1)%espe(kato))
           end do

        endif
        !    if(ncomp_qua /= nestates) then
        !    call outfor(3_ip,lun_outpu_qua,'Warning! Nmro of Eigenstates to be calculated need to be updated')  
        ncomp_qua = nestates    ! esta se usa para clacular autovaloer en Schrodinger!!!

        allocate(phion(npoin,ncomp_qua),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'PHION','qua_memall',phion)

        ! quizas deberia allocar phion aca! bueno esto queda provisorio
        !     endif
     end if

  else  

     if( INOTMASTER ) then
        allocate(phion(npoin,ncomp_qua),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'PHION','qua_memall',phion)
     end if

  endif

  print*,'Begite: end=',kfl_paral

end subroutine qua_begite
    
