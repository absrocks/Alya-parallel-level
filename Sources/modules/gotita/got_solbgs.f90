subroutine got_solbgs()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_solbgs
  ! NAME 
  !    got_solite
  ! DESCRIPTION
  !    This routine solves an iteration for the incompressible NS equations
  !    using a monolitic scheme.
  !    Velocity unknowns are in UNKNO(1:NDIME*NPOIN)
  !    Alpha    unknowns are in UNKNO(NDIME*NPOIN+1:) 
  ! USES
  !    got_matrix
  !    Soldir
  !    Solite
  !    got_rotunk
  ! USED BY
  !    got_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  implicit none
  integer(ip)          :: itbgs,kfl_gobgs_got,idofn,icomp,ndofn,idime,ipoin
  integer(ip), save    :: ipass=0
  real(rp)             :: rgsve,rgsal,rgsve_old,rgsal_old,dt
  real(rp)             :: rgsve_ini,rgsal_ini,rels1_got
  real(rp),    pointer :: unkno_old(:)

  kfl_gobgs_got = 1
  itbgs         = 0
  rels1_got     = 1.0_rp-relgs_got
  ndofn         = ndime*npoin
  rgsve         = 0.0_rp
  rgsal         = 0.0_rp

  do while(kfl_gobgs_got==1)

     itbgs=itbgs+1
     rgsve_old=rgsve
     rgsal_old=rgsal
     !-------------------------------------------------------------------
     !
     ! Solve the momentum equations
     ! 
     !-------------------------------------------------------------------
     if(kfl_paral/=0.and.(relgs_got/=1.0_rp.or.mibgs_got/=1)) then
        do idofn=1,ndime*npoin
           unbgs_got(idofn)=unkno(idofn)
        end do
     end if
     call got_inisol(2_ip)
     !if(kfl_paral<=0) write(lun_solve,105) itbgs
     call got_matrix(2_ip)               
     call solver(rhsid,unkno,amatr,pmatr)

     if(relgs_got==1.0_rp.and.mibgs_got==1) then
        rgsve=1.0_rp
     else 
        if(kfl_paral/=0) then
           do icomp=1,ndime*npoin
              unkno(icomp) = relgs_got*unkno(icomp)&
                   &       + rels1_got*unbgs_got(icomp)
           end do
        end if
        call residu(kfl_normc_got,ndime,ndime,unkno,&
             unbgs_got,one,one,ndime,1.0_rp,rgsve)
     end if
     !-------------------------------------------------------------------
     !
     ! Solve the continuity equation
     !
     !-------------------------------------------------------------------
     if(kfl_paral/=0.and.(relgs_got/=1.0_rp.or.mibgs_got/=1)) then
        do idofn=1,npoin
           unbgs_got(idofn)=unkno(ndofn+idofn)
        end do
     end if
     call got_inisol(3_ip)
     !if(kfl_paral<=0) write(lun_solve,104) itbgs
     call got_matrix(3_ip)
     if(kfl_paral/=0) then
        unkno_old => unkno 
        unkno     => unkno(ndime*npoin+1:)
        call solver(rhsid,unkno,amatr,pmatr)
        unkno     => unkno_old
     else
        call solver(rhsid,unkno,amatr,pmatr)
     end if

     if(relgs_got==1.0_rp.and.mibgs_got==1) then
        rgsal=1.0_rp
     else 
        if(kfl_paral/=0) then
           do icomp=1,npoin
              idofn=ndofn+icomp
              unkno(idofn) = relgs_got*unkno(idofn)&
                   &       + rels1_got*unbgs_got(icomp)
           end do
        end if
        call residu(kfl_normc_got,one,one,unkno(ndime*npoin+1),&
             unbgs_got,one,one,ndime,1.0_rp,rgsal)
     end if
     !-------------------------------------------------------------------
     !
     ! Check convergence of the Block Gauss-Seidel 
     ! 
     !-------------------------------------------------------------------
     if(rgsve<tobgs_got.or.itbgs>=mibgs_got) kfl_gobgs_got = 0 
     if(itbgs==1) then
        rgsve_ini=rgsve
        rgsal_ini=rgsal
     end if
     !if(kfl_blogs_got==1.and.kfl_paral<=0) then
     !   if(ipass==0) then
     !      ipass=1
     !      write(lun_blogs_got,100)
     !   end if
     !   write(lun_blogs_got,101)&
     !        ittim,itcou,itinn(modul),itbgs,cutim,rgsve,rgsal
     !end if

  end do
  !
  ! Output BGS convergence
  !
  if(kfl_paral==0) then
     !write(lun_solve,200)
     if(rgsve_old>0.0_rp) then
        !write(lun_solve,202) &
        !     rgsve_ini,rgsve,&
        !     log10(rgsve_old/rgsve),itbgs
     else
        !write(lun_solve,202) &
        !     rgsve_ini,rgsve,&
        !     1.0_rp,itbgs       
     end if
     !write(lun_solve,201)
     if(rgsal_old>0.0_rp) then
        !write(lun_solve,202) &
        !     rgsal_ini,rgsal,&
        !     log10(rgsal_old/rgsal),itbgs
     else
        !write(lun_solve,202) &
        !     rgsal_ini,rgsal,&
        !     1.0_rp,itbgs       
     end if
  end if
  !
  ! Formats
  !
100 format('# --| ALYA NASTIN: Block Gauss-Seidel convergence  ' ,/,&
       & '# --| Columns displayed:' ,/,&
       & '# --| 1. Time Step         2. Global Iteration  3. Inner Iteration   ',/,&
       & '# --| 4. Block GS Iter.    5. Current time      6. Velocity          ',/,& 
       & '# --| 7. Alpha             ',//,&
       & '# ','          1','          2','          3',&
       &      '          4','             5','             6','             7') 
101 format(4x,i9,2x,i9,2x,i9,2x,i9,10(2x,e12.6))
104 format(/,&
       10x,'    CONTINUITY EQ, BGS ITERATION=',i3,/,&
       10x,'    -----------------------------')
105 format(/,&
       10x,'    MOMENTUM EQS,  BGS ITERATION=',i3,/,&
       10x,'    -----------------------------')
200 format(/,&
       10x,'GAUSS-SEIDEL CONVERGENCE, MOMENTUM EQS:')
201 format(/,&
       10x,'GAUSS-SEIDEL CONVERGENCE, CONTINUITY EQ:')
202 format(&
       10x,'    INITIAL RESIDUAL NORM  : ',e12.5,/,&
       10x,'    FINAL RESIDUAL NORM    : ',e12.5,/,&
       10x,'    CONVERGENCE RATE       : ',e12.5,/,&
       10x,'    NUMBER OF ITERATIONS   : ',i12)

end subroutine got_solbgs
