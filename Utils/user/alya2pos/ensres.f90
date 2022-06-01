subroutine ensres(&
     ittim,npoin,nelem,parr1,pari1,lexis,&
     lbxis,ltype,leinv,ndime,namda,wopos,scavec,scatyp,intrea,&
     wbytes,kfl_filt,rttim,kfl_markm,kfl_multi,kfl_field)

  use def_kintyp, only          :  ip,rp,varna_pos,varnu_pos,nunam_pos
  use def_kintyp, only          :  tipoe_ens,nppti_ens,ncoun_pos
  use def_elmtyp
  use def_kintyp, only          :  nmax_ensi

  implicit none
  integer(ip),    intent(in)    :: ittim,npoin,nelem,ndime
  integer(ip),    intent(in)    :: lexis(*)
  integer(ip),    intent(in)    :: lbxis(*)
  integer(ip),    intent(in)    :: ltype(*)
  integer(ip),    intent(in)    :: leinv(*)
  real(rp),       intent(in)    :: rttim
  real(rp),       intent(in)    :: parr1(*)
  integer(ip),    intent(in)    :: pari1(*)
  character(150), intent(in)    :: namda
  character(5),   intent(in)    :: wopos(*)
  character(5),   intent(in)    :: scavec
  character(5),   intent(in)    :: scatyp
  character(5),   intent(in)    :: kfl_filt
  character(5),   intent(in)    :: intrea
  character(5),   intent(in)    :: wbytes
  integer(ip),    intent(in)    :: kfl_markm,kfl_field
  integer(ip),    intent(in)    :: kfl_multi
  character(8)                  :: chtim
  character(20)                 :: wopo2(3)
  integer(ip)                   :: ipoin,ielem,ielty,idime,i,kpoin,ipos
  integer(ip)                   :: mpoin,nbyte,jelty,istpp,itise
  integer(ip)                   :: iesta,iesto
  real(rp)                      :: dummr
  character(150)                :: filva

  if( ittim >= 0 ) then

     istpp=1
     do while(istpp<=nmax_ensi)
        if(wopos(1)==varna_pos(1,istpp)) then
           istpp=nmax_ensi
        else if(varna_pos(1,istpp)=='NULL') then
           varna_pos(1,istpp) = wopos(1)        ! OUTVAR NAME
           varna_pos(2,istpp) = wopos(2)        ! SCALA/VECTOR
           varna_pos(3,istpp) = kfl_filt        ! FILTE/NOFIL
           varna_pos(4,istpp) = scatyp          ! NELEM/NPOIN/NBOUN
           varnu_pos          = varnu_pos+1
           istpp=nmax_ensi
        end if
        istpp=istpp+1
     end do

     istpp=1
     do while(istpp<=10000)
        if(rttim==tipoe_ens(istpp)) then
           istpp=10000
        else if(tipoe_ens(istpp)==-1.0_rp) then
           tipoe_ens(istpp)= rttim
           nppti_ens = nppti_ens + 1
           istpp     = 10000
        end if
        istpp=istpp+1
     end do

     write(nunam_pos,'(i6)') nppti_ens      !   <<-- to write an integer to a character
     if(nppti_ens<10) then
        write(nunam_pos,'(a,i1)') '00000',nppti_ens
     else if(nppti_ens<100) then
        write(nunam_pos,'(a,i2)') '0000',nppti_ens
     else if(nppti_ens<1000) then
        write(nunam_pos,'(a,i3)') '000',nppti_ens
     else if(nppti_ens<10000) then
        write(nunam_pos,'(a,i4)') '00',nppti_ens
     else if(nppti_ens<100000) then
        write(nunam_pos,'(a,i5)') '0',nppti_ens
     end if
     !
     ! Write result
     !
     ncoun_pos = ncoun_pos + 1
     if(ittim<10) then
        write(chtim,'(a,i1)') '0000000',ittim
     else if(ittim<100) then
        write(chtim,'(a,i2)') '000000',ittim
     else if(ittim<1000) then
        write(chtim,'(a,i3)') '00000',ittim
     else if(ittim<10000) then
        write(chtim,'(a,i4)') '0000',ittim
     else if(ittim<100000) then
        write(chtim,'(a,i5)') '000',ittim
     else if(ittim<1000000) then
        write(chtim,'(a,i6)') '00',ittim
     else if(ittim<10000000) then
        write(chtim,'(a,i7)') '0',ittim
     end if

     !
     ! Element types
     !
     if (ndime == 2_ip) then
        iesta=10
        iesto=29
     else if(ndime == 3_ip) then
        iesta=30
        iesto=50
     end if

     !
     ! Open file
     !
     filva = trim(namda)//'.ensi.'//trim(wopos(1))//'-'//trim(nunam_pos)
     open(unit=102,file=trim(filva),status='unknown',form='formatted')

     if( scavec == 'SCALA' ) then
        !
        ! SCALAR
        !
        if (scatyp == 'NELEM') then
           !
           ! Scalar per-element
           !
           write(102,100) 'Alya Ensight Gold --- Scalar per-element variables file'
           write(102,100) 'part'
           write(102,110) 1_ip
           write(102,100) 'block'
           do ielty=iesta,iesto    ! Loop element types
              do ielem = 1,nelem
                 if (abs(ltype(ielem)) == ielty) then
                    write(102,120) parr1(ielem)
                 end if
              end do
           end do

        else
           !
           ! Scalar per-node
           !
           write(102,100) 'Alya Ensight Gold --- Scalar per-node variables file'
           write(102,100) 'part'
           write(102,110) 1_ip
           write(102,100) 'coordinates'
           do ipoin = 1,npoin
              write(102,120) parr1(ipoin)
           end do
        end if

     else if( scavec == 'VECTO' ) then
        !
        ! VECTOR
        !
        if (scatyp == 'NELEM') then
           !
           ! Vector per-element
           !
           write(102,100) 'Alya Ensight Gold --- Vector per-element variables file'
           write(102,100) 'part'
           write(102,110) 1_ip
           write(102,100) 'block'
           do ielty=iesta,iesto    ! Loop element types
              do idime = 1,ndime
                 do ielem = 1,nelem
                    ipos = ndime*(ielem - 1) + idime
                    if (abs(ltype(ielem)) == ielty) then
                       write(102,120) parr1(ipos) ! Ordered as vx_e1, vx_e2, ... vx_en, vy_e1, vy_e2 ...
                    end if
                 end do
              end do
           end do

        else
           !
           ! Vector per-node
           !
           if (kfl_field==0) then
              write(102,100) 'Alya Ensight Gold --- Vector per-node variables file'
              write(102,100) 'part'
              write(102,110) 1_ip
              write(102,100) 'coordinates'
              do idime = 1,ndime
                 kpoin = idime-ndime
                 do ipoin = 1,npoin
                    kpoin = kpoin + ndime
                    write(102,120) parr1(kpoin)
                 end do
              end do
           else if (kfl_field == 1) then
              ! esto es muy util para crear un field con desplazamientos...
              do ipoin = 1,npoin
                 kpoin = (ipoin-1)*ndime + 1
                 write(102,9000) ipoin,parr1(kpoin),parr1(kpoin+1),parr1(kpoin+2)
              end do
9000          format(i9,2x,3(2x,e16.8E3))
           end if
        end if

        if (ndime == 2) then
           dummr = 0.0_rp
           do ipoin = 1,npoin
              write(102,120) dummr
           end do
        end if

     end if
     !
     ! Close file
     !
     close(unit=102)

  else if( ittim == -1 ) then
     !
     ! Rewrite geometry file
     !
     rewind(100)
     itise = 1
     !tipoe_ens(nppti_ens) = cutim
     rewind(100)
     write(100,'(a)') '#'
     write(100,'(a)') '# Alya generated post-process files'
     write(100,'(a)') '# Ensight Gold Format'
     write(100,'(a)') '#'
     write(100,  50 ) '# Problem name:   ', adjustl(trim(namda))
     write(100,'(a)') '#'
     write(100,'(a)') 'FORMAT'
     write(100,'(a)') 'type:    ensight gold'
     write(100,'(a)') 'GEOMETRY'
     write(100,  60 ) 'model:   ', itise,  adjustl(trim(namda))//'.ensi.geo'
     write(100,'(a)') 'VARIABLE'
     flush(100)
     !
     ! Write case end of file
     !
     itise = 1
     if( ncoun_pos /= 0 ) then
        do istpp=1,varnu_pos
           !
           ! write the generic file name
           !
           filva=adjustl(trim(namda))//'.ensi.'//trim(varna_pos(1,istpp))//'-******'

           if(varna_pos(2,istpp)=='VECTO'.and.varna_pos(3,istpp)/='FILTE') then
              !
              ! SCALA
              !
              if (varna_pos(4,istpp) == 'NELEM') then
                 ! NELEM
                 write(100,70) 'vector per element:',itise,varna_pos(1,istpp),adjustl(trim(filva))
              else
                 ! NPOIN
                 write(100,70) 'vector per node:',itise,varna_pos(1,istpp),adjustl(trim(filva))
              end if
           else if(varna_pos(2,istpp)=='SCALA'.and.varna_pos(3,istpp)/='FILTE') then
              !
              ! SCALA
              !
              if (varna_pos(4,istpp) == 'NELEM') then
                 ! NELEM
                 write(100,70) 'scalar per element:',itise,varna_pos(1,istpp),adjustl(trim(filva))
              else
                 ! NPOIN
                 write(100,70) 'scalar per node:',itise,varna_pos(1,istpp),adjustl(trim(filva))
              end if
           end if

        end do

        write(100,'(a)') 'TIME'
        write(100, 80  ) 'time set:               ',itise
        write(100, 80  ) 'number of steps:        ',nppti_ens
        write(100,'(a)') 'filename start number:        1 '
        write(100,'(a)') 'filename increment:           1 '
        write(100,'(a)') 'time values: '
        write(100,'(10(1x,e12.5))') (tipoe_ens(i),i=1,nppti_ens)
        flush(100)

        !nppti_ens     = nppti_ens + 1
        !nppva_ens     = 0
        !kfl_statu_ens = 2
        !ncoun_pos     = 0
     end if
  end if

10 format(a)
15 format(2a)
20 format(i10)
25 format(20i10)
30 format(e12.5)
50 format(2a)
60 format(a,3x,i4,3x,a)
70 format(a,4x,i4,4x,a,4x,a)
80 format(a,3x,i8)
100 format(a)
110 format(i10)
120 format(e16.8E3)
122 format(f20.8)

end subroutine ensres
