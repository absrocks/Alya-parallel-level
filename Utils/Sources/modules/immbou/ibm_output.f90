subroutine ibm_output()
  !------------------------------------------------------------------------
  !****f* Immbou/ibm_output
  ! NAME 
  !    ibm_output
  ! DESCRIPTION
  !    Output and postprocess of solution
  ! USES
  ! USED BY
  !    ibm_iniunk
  !    ibm_endite
  !    ibm_endste
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use mod_memchk
  use mod_iofile
  use mod_outfor, only : outfor
  implicit none
  integer(ip)          :: ivari,ivarp,ipois,iboib,ieles,ipoib,idime
  integer(ip)          :: inoib,iimbo,lunit
  integer(ip), pointer :: ntrad(:)
  integer(ip), save    :: ipass = 0
  real(rp),    pointer :: F(:,:),Fp(:),Fv(:),a(:,:),v(:,:),x(:,:)   ! Linear motion
  real(rp),    pointer :: T(:,:),Tp(:),Tv(:),z(:,:),w(:,:),s(:,:)   ! Angular motion
  real(rp),    pointer :: Mass,Rho,Vol,Momi(:),COG(:)
  integer(4)           :: istat
  character(7)         :: statu
  character(11)        :: forma
  character(6)         :: posit

  !---------------------------------------------------------------
  !
  ! Open IB files and output IB mesh
  !
  !---------------------------------------------------------------

  if( ittyp == ITASK_INITIA ) then

     call ibm_openfi(1_ip)
goto 10
     if( kfl_rstar == 0 .and. kfl_outfo /= 0 ) then

        if( kfl_outfo == 1 ) then

           if( 1 == 0 ) then
              statu='old'
              forma='formatted'
              posit='append'               
              call iofile(zero,lun_outpu_dom,fil_outpu_dom,'MESH',statu,forma,posit)
              ipois = 0
              ieles = 0
              ipois = npoin
              ieles = nelem
              do iimbo = 1,nimbo
                 call memgen(1_ip,imbou(iimbo)%nboib,0_ip)
                 do iboib = 1,imbou(iimbo)%nboib
                    gisca(iboib) = iimbo
                 end do
                 !gisca => imbou(iimbo)%ltyib
                 lunit = lun_outpu_dom
!                 call gidele(&
!                      ibsta_dom,ibsto_dom,mnoib,imbou(iimbo)%npoib,imbou(iimbo)%nboib,lunit,lexib,&
!                      imbou(iimbo)%ltyib,imbou(iimbo)%lnoib,imbou(iimbo)%cooin,gisca,ieles,ipois,'IBM',&
!                      kfl_markm)
                 ipois = ipois + imbou(iimbo)%npoib 
                 ieles = ieles + imbou(iimbo)%nboib
                 call memgen(3_ip,imbou(iimbo)%nboib,0_ip)
              end do
              call iofile(two,lun_outpu_dom,fil_outpu_dom,'MESH')

           else

              ipois = 0
              ieles = 0
              ipois = npoin
              ieles = nelem
              do iimbo = 1,nimbo
                 call memgen(1_ip,imbou(iimbo)%nboib,0_ip)
                 do iboib = 1,imbou(iimbo)%nboib
                    gisca(iboib) = iimbo
                 end do
                 !gisca => imbou(iimbo)%ltyib
                 lunit = lun_mshib
                
!                 call gidele(&
!                      ibsta_dom,ibsto_dom,mnoib,imbou(iimbo)%npoib,imbou(iimbo)%nboib,lunit,lexib,&
!                      imbou(iimbo)%ltyib,imbou(iimbo)%lnoib,imbou(iimbo)%cooin,gisca,ieles,ipois,'IBM',&
!                      kfl_markm)
                 ipois = ipois + imbou(iimbo)%npoib 
                 ieles = ieles + imbou(iimbo)%nboib
                 call memgen(3_ip,imbou(iimbo)%nboib,0_ip)
              end do
           end if

        else if( kfl_outfo == 300000 ) then

           allocate(ntrad(nimbo),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'NTRAD','ibm_output',ntrad)
           nullify(ntrad)
           ntrad(1) = 0
           do iimbo = 2,nimbo
              ntrad(iimbo) = ntrad(iimbo-1) + imbou(iimbo) % npoib
           end do
           write(lun_mshib) ndime,npoib,nboib,0_ip
           write(lun_mshib) (( ( imbou(iimbo) % cooib(idime,ipoib), idime = 1,ndime),ipoib = 1,imbou(iimbo)%npoib),iimbo = 1,nimbo)
           write(lun_mshib) (( ( imbou(iimbo) % lnoib(inoib,iboib) + ntrad(iimbo) , inoib = 1,mnoib),   &
                            iboib = 1,imbou(iimbo)%nboib),iimbo = 1,nimbo)
           deallocate(ntrad,stat=istat)

        else if( kfl_outfo == 30 ) then

           call runend('NO LONGER SUPPORTED')

        end if

     end if
10 continue

     call ibm_openfi(3_ip)

  end if


  !---------------------------------------------------------------
  !
  ! Initial solution, end of a time step and and of run
  !
  !---------------------------------------------------------------
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)   
     call ibm_outvar(ivari)
  end do

  !---------------------------------------------------------------
  !
  ! Particle properties
  !
  !---------------------------------------------------------------

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then

     if( INOTSLAVE ) then
        !
        ! Output results
        !
        if( ipass == 0 ) write(lun_outpu_ibm,99) nimbo
        ipass    = ipass + 1
        ioutp(1) = ittim
        routp(1) = cutim
        call outfor(46_ip,lun_outpu_ibm,' ')
        if( ipass == 1 ) then
           do iimbo = 1,nimbo
              Mass  => imbou(iimbo) % massa
              Rho => imbou(iimbo) % densi
              Vol => imbou(iimbo) % volum
              Momi  => imbou(iimbo) % momin
              COG  => imbou(iimbo) % posgr
              write(momod(modul) % lun_outpu,102) iimbo,Vol,Rho,Mass
              write(momod(modul) % lun_outpu,103) iimbo,COG
              write(momod(modul) % lun_outpu,104) iimbo,Momi
           end do
           flush(momod(modul) % lun_outpu)
        end if
        
        do iimbo = 1,nimbo
           F  => imbou(iimbo) % force
           Fp => imbou(iimbo) % pforce
           Fv => imbou(iimbo) % vforce
           a  => imbou(iimbo) % accel
           v  => imbou(iimbo) % velol
           x  => imbou(iimbo) % posil
           T  => imbou(iimbo) % torqu
           Tp => imbou(iimbo) % ptorqu
           Tv => imbou(iimbo) % vtorqu
           z  => imbou(iimbo) % accea
           w  => imbou(iimbo) % veloa
           s  => imbou(iimbo) % posia
           write(lun_outpu_ibm,101) iimbo,&
                x(1,1),x(2,1),x(3,1),v(1,1),v(2,1),v(3,1),a(1,1),a(2,1),a(3,1),   &
                s(1,1),s(2,1),s(3,1),w(1,1),w(2,1),w(3,1),z(1,1),z(2,1),z(3,1),   &
                F(1,1),F(2,1),F(3,1),Fv(1),Fv(2),Fv(3),Fp(1),Fp(2),Fp(3),T(1,1),T(2,1),T(3,1),&
                Tv(1),Tv(2),Tv(3),Tp(1),Tp(2),Tp(3)
        end do
        flush(lun_outpu_ibm)
     end if

  end if
  !
  ! Formats
  !
99 format('# ALYA Boundary set results',/,&
       & '#',/,&
       & '# HEADER',/,&
       & '# ISET, Column :   1',/,&
       & '# X,    Column :   2',/,&
       & '# Y,    Column :   3',/,&
       & '# Z,    Column :   4',/,&
       & '# VX,   Column :   5',/,&
       & '# VY,   Column :   6',/,&
       & '# VZ,   Column :   7',/,&
       & '# AX,   Column :   8',/,&
       & '# AY,   Column :   9',/,&
       & '# AZ,   Column :  10',/,&
       & '# S1,   Column :  11',/,&
       & '# S2,   Column :  12',/,&
       & '# S3,   Column :  13',/,&
       & '# W1,   Column :  14',/,&
       & '# W2,   Column :  15',/,&
       & '# W3,   Column :  16',/,&
       & '# Z1,   Column :  17',/,&
       & '# Z2,   Column :  18',/,&
       & '# Z3,   Column :  19',/,&
       & '# F1,   Column :  20',/,&
       & '# F2,   Column :  21',/,&
       & '# F3,   Column :  22',/,&
       & '# Fv1,  Column :  23',/,&
       & '# Fv2,  Column :  24',/,&
       & '# Fv3,  Column :  25',/,&
       & '# Fp1,  Column :  26',/,&
       & '# Fp2,  Column :  27',/,&
       & '# Fp3,  Column :  28',/,&
       & '# T1,   Column :  29',/,&
       & '# T2,   Column :  30',/,&
       & '# T3,   Column :  31',/,&
       & '# Tv1,  Column :  32',/,&
       & '# Tv2,  Column :  33',/,&
       & '# Tv3,  Column :  34',/,&
       & '# Tp1,  Column :  35',/,&
       & '# Tp2,  Column :  36',/,&
       & '# Tp3,  Column :  37',/,&
       & '# NUMVARIABLES :  37',/,&
       & '# NUMSETS      :  ',i4,/,&
       & '# START')

101 format(4x,i9,50(2x,e12.6))

102 format('# ISET, volum,densi,massa: ',i6,3(1x,e12.5))
103 format('# ISET, posgr: ',i6,3(1x,e12.5))
104 format('# ISET, : MOM. INERC. ',i6,6(1x,e12.5))

end subroutine ibm_output
