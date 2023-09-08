subroutine qua_reaphy()
  !------------------------------------------------------------------------
  !****f* Quanty/qua_reaphy
  ! NAME 
  !    qua_reaphy
  ! DESCRIPTION
  !    This routine reads the physical problem definition
  ! USES
  ! USED BY
  !    qua_turnon
  !------------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  use def_parame
  use def_inpout
  use def_master
  use def_quanty
  use def_domain
  use mod_memchk
  use def_solver
  use mod_ecoute, only :  ecoute
  use mod_outfor, only : outfor

  implicit none
  integer(ip)            :: nat,kk,jj,nespe,nklleno,nestallena,nesta,nk,nm,nl,mm,ndata,n1,nconta,nkk,nll,nmm
  integer(4)             :: istat
  character*20           :: fname
  character*5            :: Atch


  if( INOTSLAVE ) then
     !
     ! Read and write files
     !
     lispa = 0
     lisda = lun_pdata_qua                                ! Reading file
     lisre = lun_outpu_qua                                ! Writing file
     !
     ! Initializations (defaults)
     !
     kfl_timei_qua   = 0                    !  temporal problem d()/dt
     kfl_dftgs_qua   = 1                    ! DFT ground state
     kfl_potxc_qua   = 0                    ! tipo de XC  0::ninguno
     !             1::LDA
     !             2::GDA
     kfl_alele_qua   = 0                    ! all electron problem
     kfl_nolin_qua   = 0                    !  nonlinar problem
     kfl_perio_qua   = 0                    !  periodic problem
     kfl_coulo_qua   = 0                    ! Coulomb potential
     kfl_bfieldx_qua = 0                    ! B field
     kfl_bfieldy_qua = 0                    ! B field
     kfl_bfieldz_qua = 0                    ! B field
     klf_btemp_law   = 0                    ! B temp. law
     btemp_law_qua   = 0                    ! B temp. law
     kfl_efieldx_qua = 0                    ! E field
     kfl_efieldy_qua = 0                    ! E field
     kfl_efieldz_qua = 0                    ! E field
     klf_etemp_law   = 0                    ! E temp. law
     kfl_spinb_qua   = 0                    ! Spin B field
     kfl_spiorb_qua  = 0                    ! Spin orbit
     kfl_relat_qua   = 0                    ! relativistic
     kfl_vother_qua  = 0                    ! Other potential
     lawma_qua       = 0                    ! Law mass
     ncuanpal_qua    = 1                    ! Quantum main number N
     lcuanorb_qua    = 0                    ! Orbital Quantum number L (=0,..,N)
     nspin_qua       = 1                    ! Spin Quantum number S  (S=1 up, S=0 down)

     etemp_law_qua  = 0                     ! E temp. law
     eig_evol_qua   = 1.0_rp                ! autovalor que evoluciona
     coulo_qua      = 1.0_rp                ! Value of Coulomb potential
     bfieldx_qua    = 0.0_rp                ! Value of B field en x
     bfieldy_qua    = 0.0_rp                ! Value of B field en y
     bfieldz_qua    = 0.0_rp                ! Value of B field en z
     efieldx_qua    = 0.0_rp                ! Value of E field en x
     efieldy_qua    = 0.0_rp                ! Value of E field en y
     efieldz_qua    = 0.0_rp                ! Value of E field en z
     frecuencie     = 0_rp                  ! frecuencia temporal de E o B
     spinb_qua      = 0.0_rp                ! Value of Spin B field
     spiorb_qua     = 0.0_rp                ! Value of Spin orbit
     relat_qua      = 0.0_rp                ! Value of relativistic
     vother_qua     = 0.0_rp                ! Value of Other potential
     law_vother_qua = 0.0                   ! Law Other potential
     w_vother_qua   = 0.0_rp
     x0_qua         = 0.0_rp                ! X_0 position of the Nucleous 
     y0_qua         = 0.0_rp                ! Y_0 position of the Nucleous 
     z0_qua         = 0.0_rp                ! Z_0 position of the Nucleous 
     massa_qua      = 1.0_rp                ! Mass


     !
     ! Reach the section
     !
     call ecoute('qua_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('qua_reaphy')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDPH')

        if(words(1)=='NESPH') then                    ! time depedent problem swicht
           if(words(2)=='ON   ') then
              kfl_spher = 1
           else
              kfl_spher = 0
           end if

        elseif(words(1)=='TIMED') then                    ! time depedent problem swicht
           if(words(2)=='ON   ') then
              kfl_timei_qua = 1
           else
              kfl_timei_qua = 0
           end if
           eig_evol_qua=getrea('VALUE',1.0_rp,'#autovalor evolucionado')

        elseif(words(1)=='DFTGS') then                    ! DFT GROUND STATE
           if(words(2)=='ON   ') then
              kfl_dftgs_qua = 1
           else
              kfl_dftgs_qua = 0
           end if
           kfl_potxc_qua =getint('VALUE',0_ip,'#XC potential')

        elseif(words(1)=='ALLEE') then                    ! DFT all electron
           if(words(2)=='ON   ') then
              kfl_alele_qua = 1
           else
              kfl_alele_qua = 0
           end if
           !kfl_potxc_qua =getint('VALUE',0_ip,'#XC potential')

        elseif(words(1)=='COULO') then                    ! Coulomb potential
           if(words(2)=='ON   ') then
              kfl_coulo_qua = 1
           else
              kfl_coulo_qua = 0
           end if
           coulo_qua=getrea('VALUE',1.0_rp,'#Coulomb potential')

        else if(words(1)=='BFLDX') then               ! B FIELD x
           if(words(2)=='ON   ') then
              kfl_bfieldx_qua= 1
           else
              kfl_bfieldx_qua= 0
           end if
           bfieldx_qua=getrea('VALUE',0.0_rp,'#B field x')

        else if(words(1)=='BFLDY') then               ! B FIELD y
           if(words(2)=='ON   ') then
              kfl_bfieldy_qua= 1
           else
              kfl_bfieldy_qua= 0
           end if
           bfieldy_qua=getrea('VALUE',0.0_rp,'#B field y')

        else if(words(1)=='BFLDZ') then               ! B FIELD z
           if(words(2)=='ON   ') then
              kfl_bfieldz_qua= 1
           else
              kfl_bfieldz_qua= 0
           end if
           bfieldz_qua=getrea('VALUE',0.0_rp,'#B field z')

        else if(words(1)=='BTEMP') then               ! TEMPORAL LAW OF B
           if(words(2)=='ON   ') then
              klf_btemp_law= 1
           else
              klf_btemp_law= 0
           end if
           btemp_law_qua=getrea('VALUE',1.0_rp,'#B temp. law')

        else if(words(1)=='EFLDX') then                ! E FIELD x
           if(words(2)=='ON   ') then
              kfl_efieldx_qua= 1
           else
              kfl_efieldx_qua= 0
           end if
           efieldx_qua=getrea('VALUE',0.0_rp,'#E field x') 

        else if(words(1)=='EFLDY') then                ! E FIELD y
           if(words(2)=='ON   ') then
              kfl_efieldy_qua= 1
           else
              kfl_efieldy_qua= 0
           end if
           efieldy_qua=getrea('VALUE',0.0_rp,'#E field y') 

        else if(words(1)=='EFLDZ') then                ! E FIELD z
           if(words(2)=='ON   ') then
              kfl_efieldz_qua= 1
           else
              kfl_efieldz_qua= 0
           end if
           efieldz_qua=getrea('VALUE',0.0_rp,'#E field z') 

        else if(words(1)=='ETEMP') then               ! TEMPORAL LAW OF E
           if(words(2)=='ON   ') then
              klf_etemp_law= 1
           else
              klf_etemp_law= 0
           end if

           etemp_law_qua=getrea('VALUE',1.0_rp,'#E temp. law')

        else if(words(1)=='SPINB') then               ! Spin - B interaction
           if(words(2)=='ON   ') then
              kfl_spinb_qua= 1
           else
              kfl_spinb_qua= 0
           end if
           spinb_qua=getrea('VALUE',0.0_rp,'#SPIN B int.')

        else if(words(1)=='SPINO') then               ! Spin - Orbita interaction
           if(words(2)=='ON   ') then
              kfl_spiorb_qua= 1
           else
              kfl_spiorb_qua= 0
           end if
           spiorb_qua=getrea('VALUE',0.0_rp,'#SPIN orbita int.')

        else if(words(1)=='RELAT') then               ! RELATIVISTIC
           if(words(2)=='ON   ') then
              kfl_relat_qua= 1
           else
              kfl_relat_qua= 0
           end if
           relat_qua=getrea('VALUE',0.0_rp,'#RELATIVISTIC')

        else if(words(1)=='OTHER') then               ! Other potential
           if(words(2)=='ON   ') then
              kfl_vother_qua= 1
           else
              kfl_vother_qua= 0
           end if
           vother_qua=getrea('VALUE',0.0_rp,'#OTHER POTENTIAL')

        else if(words(1)=='LAWOT') then               ! law Other potential
           if(words(2)=='ON   ') then
              law_vother_qua=1
           else
              law_vother_qua=0
           end if
           w_vother_qua=getrea('VALUE',0.0_rp,'#OTHER POTENTIAL')

           !  else if(words(1)=='ATOMS') then               ! ATOMOS ?
           !     nat = getint('VALUE',0_ip,'#atomos')
           !    atomo_qua(1)%natoms = nat
           !  else if(words(1)=='ESPEC') then               ! ESPECIES ?
           !     nespe = getint('VALUE',0_ip,'#especies')
           !    atomo_qua(1)%nespec = nespe
        else if(words(1)=='FILEX') then               ! file ATOMOS 

           atomo_qua(1)%wprob=getcha('FILE ','atoms.xyz','#file name')
           atomo_qua(1)%wprob=atomo_qua(1)%wprob(1:5)//'.xyz'
           !'atoms.xyz'


           ! leo aca! y san seacabo
           ! abro file 
           fname=adjustl(atomo_qua(1)%wprob)
           open(unit=111,file= fname,status = 'unknown')

           read(111,'(i3)') atomo_qua(1)%natoms  ! numero total de atomos
           read(111,'(i3)') atomo_qua(1)%nespec  ! numero total de especies

           call qua_memphy(1_ip)
           call qua_memphy(2_ip)
 
           ! leo atoms.xyz
           do kk=1,atomo_qua(1)%natoms
              read(111,*) &
                   atomo_qua(1)%espe(kk),atomo_qua(1)%Tiespe(kk),atomo_qua(1)%coorx(kk), &
                   atomo_qua(1)%coory(kk),atomo_qua(1)%coorz(kk),atomo_qua(1)%spin(kk),atomo_qua(1)%tipoPP(kk)        
           enddo

           !     read(111,*) ndata  ! numero de especies

           nestallena=0
           do kk=1,atomo_qua(1)%nespec 
              read(111,*) n1,Atch,nesta 
              if(Atch=='empty') then
                 nestallena=nestallena+0
                 nklleno=0
              elseif(Atch=='He') then
                 nestallena=nestallena+1
                 nklleno=1
              elseif(Atch=='Ne') then
                 nestallena=nestallena+5
                 nklleno=2
              elseif(Atch=='Ar') then
                 nestallena=nestallena+9
                 nklleno=3
              elseif(Atch=='Kr') then
                 nestallena=nestallena+18 
                 nklleno=4
              elseif(Atch=='Xe') then
                 nestallena=nestallena+27 
                 nklleno=5
              elseif(Atch=='Rn') then
                 nestallena=nestallena+43 
                 nklleno=6
              endif

 
              nconta=0
              do nk=1,nklleno
                 do nl=0,nk-1
                    do nm=1,2*nl+1
                       nconta=nconta+1
                       especie_qua(n1)%nocupa(nconta)=2.0_rp
                    enddo
                 enddo
              enddo

              do jj=nestallena+1,nestallena+nesta
                 read(111,*) nk,nl,mm,especie_qua(n1)%nocupa(jj)
              enddo

           enddo
           ! cierro atoms.xyz

         close(111)

           ! leyó la fila 'atoms.xyz'  


        else if(words(1)=='X0POS') then               ! X0 POSITION
           x0_qua=getrea('VALUE',0.0_rp,'#X0')
        else if(words(1)=='Y0POS') then               ! Y0 POSITION
           y0_qua=getrea('VALUE',0.0_rp,'#Y0')
        else if(words(1)=='Z0POS') then               ! Z0 POSITION
           z0_qua=getrea('VALUE',0.0_rp,'#Z0')

        else if(words(1)=='MASSA') then               ! massa         
           massa_qua=getrea('MASSA',0.5_rp,'#Mass')

        else if(words(1)=='LAWMA') then               ! Arrays interpolation      
           if(words(2)=='CONST') then
              lawma_qua=0
           end if

        else if(words(1)=='NCPPA') then               ! Quantum main number N        
           ncuanpal_qua=getint('NCPPA',1_ip,'#NCUANP')

        else if(words(1)=='NORBI') then               ! Orbital Quantum number l        
           lcuanorb_qua=getint('NORBI',0_ip,'#NORBI')

        else if(words(1)=='NSPIN') then               ! Spin Quantum number l        
           nspin_qua=getint('NSPIN',1_ip,'#NSPIN')

        end if

        call ecoute('qua_reaphy')

     end do

     call armo_estructuras()

  end if

end subroutine qua_reaphy

subroutine  armo_estructuras()

  use def_kintyp, only    :  ip,rp
  use def_parame
  use def_inpout
  use def_master
  use def_quanty
  use def_domain
  use mod_memchk
  use def_solver

  implicit none
  integer(ip)     :: kato,kk,kespe,cont_espe,nespec,istat,kocup,nspin_tot

  cont_espe=0
  nespec   =0

  if(kfl_dftgs_qua==1 .or. kfl_alele_qua ==1) then

     do kato=1,atomo_qua(1)%natoms

        do kespe=1,kato-1  
           if(atomo_qua(1)%espe(kato) /= atomo_qua(1)%espe(kespe) ) then
              cont_espe=cont_espe+1
           endif
        enddo

        if(cont_espe==kato-1) then
	   nespec=nespec + 1 
	   especie_qua(nespec)%atnumber = atomo_qua(1)%espe(kato)  
	   especie_qua(nespec)%nspin    = atomo_qua(1)%spin(kato)          
           call qua_especiename(especie_qua(nespec)%atnumber,atomo_qua(1)%tipoPP(kato),especie_qua(nespec)%wprob1,especie_qua(nespec)%wprob2)
           cont_espe=0
        endif

     enddo


     ! para chequear
     if(nespec/=atomo_qua(1)%nespec) then
        call outfor(4_ip,lun_outpu_qua,'CANNOT SURE THE DATA FOR ALL THE SPECIES REQUIRE')
     endif


     ! open PP files, upload it and determination of species, locas, etc
     lmaximo=0
     do kato=1,atomo_qua(1)%nespec

        if(kfl_alele_qua/=1) then

           ! si es NO all-electron

           call leeespecie(kato) 

           if(lmaximo.lt.especie_qua(kato)%nlmax) lmaximo = especie_qua(kato)%nlmax
           !	  call leeespecie(especie_qua(kato)%wprob2,especie_qua(kato)%atnumber,especie_qua(kato)%ncp, &
           !	  especie_qua(kato)%ncl,especie_qua(kato)%ncm,especie_qua(1)%nlmax,especie_qua(1)%nlocal, &
           !	  especie_qua(kato)%valencia,especie_qua(kato)%nrad,especie_qua(kato)%radio,especie_qua(kato)%ppseu, &
           !	  especie_qua(kato)%ppphi,especie_qua(kato)%nspin) 

        else
           ! perform all electron calculation then
           call qua_idespe(especie_qua(kato)%atnumber,especie_qua(kato)%nspin,especie_qua(kato)%ncp, &
                especie_qua(kato)%ncl,especie_qua(kato)%ncm)

	   especie_qua(atomo_qua(1)%Tiespe(kato))%valencia = especie_qua(kato)%atnumber

        endif

     enddo

  endif

  return
end subroutine armo_estructuras
