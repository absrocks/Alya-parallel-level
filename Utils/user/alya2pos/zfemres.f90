!-----------------------------------------------------------------------
!> @addtogroup Alya2pos
!> @{
!> @file    zfemres.f90
!> @author  Mariano Vazquez
!> @brief   Postprocess tool
!> @date    16/11/1966
!> @details Postprocess tool
!> @} 
!-----------------------------------------------------------------------
subroutine zfemres(&
     ittim,npoin,nelem,parr1,pari1,lexis,&
     lbxis,ltype,pdime,namda,wopos,scavec,intrea,wbytes,kfl_filt,&
     rttim,kfl_markm,kfl_multi) 


  use def_kintyp, only          :  ip,rp,varna_pos,varnu_pos,nunam_pos
  use def_kintyp, only          :  tipoe_ens,nppti_ens,ncoun_pos
  use def_elmtyp
  use def_kintyp, only          :  nmax_ensi

  implicit none
  integer(ip),    intent(in)    :: ittim  !< 
  integer(ip),    intent(in)    :: npoin  !< 
  integer(ip),    intent(in)    :: nelem  !< 
  integer(ip),    intent(in)    :: pdime  !< 
  integer(ip),    intent(in)    :: lexis(*)!< 
  integer(ip),    intent(in)    :: lbxis(*)!< 
  integer(ip),    intent(in)    :: ltype(*)!< 
  real(rp),       intent(in)    :: rttim!< 
  real(rp),       intent(in)    :: parr1(*)!< 
  integer(ip),    intent(in)    :: pari1(*)!< 
  character(150), intent(in)    :: namda!< 
  character(5),   intent(in)    :: wopos(*)!< 
  character(5),   intent(in)    :: scavec!< 
  character(5),   intent(in)    :: kfl_filt!< 
  character(5),   intent(in)    :: intrea!< 
  character(5),   intent(in)    :: wbytes!< 
  integer(ip),    intent(in)    :: kfl_markm!< 
  integer(ip),    intent(in)    :: kfl_multi!< 
  character(8)                  :: chtim
  character(20)                 :: wopo2(3)
  integer(ip)                   :: ipoin,ielem,ielty,idime,i,kpoin
  integer(ip)                   :: mpoin,nbyte,jelty,istpp,itise,ititi
  real(rp)                      :: dummr
  character(150)                :: filva


  if( ittim >= 0 ) then
    
     istpp=1
     do while(istpp<=nmax_ensi)        
        if(wopos(1)==varna_pos(1,istpp)) then
           istpp=nmax_ensi
        else if(varna_pos(1,istpp)=='NULL') then
           varna_pos(1,istpp) = wopos(1)
           varna_pos(2,istpp) = wopos(2)
           varna_pos(3,istpp) = kfl_filt
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
     filva = trim(namda)//'.zfem.'//trim(wopos(1))//'-'//trim(nunam_pos)
     open(unit=102,file=trim(filva),status='unknown',form='formatted')
     write(102,'(a)') 'REAL'
     filva= '<'//trim(wopos(1))//'>'
     write(102,'(a)') trim(filva)

     if( scavec == 'SCALA' ) then
        ititi=1
        write(102,*) ititi
        do ipoin = 1,npoin
           write(102,*) parr1(ipoin)
        end do
     else if( scavec == 'VECTO' ) then
        ititi=3
        write(102,*) ititi
        do ipoin = 1,npoin
           kpoin= (ipoin-1)*3
           write(102,*) parr1(kpoin+1),parr1(kpoin+2),parr1(kpoin+3)
        end do
     end if   

     close(unit=102)

  else if( ittim == -1 ) then
     !
     ! Rewrite geometry file
     !
     rewind(100)
     itise                = 1
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
!    if (kfl_markm == 4) then
!       write(100,  50 ) 'scalar per element:   ', adjustl(trim(title))//'.ensi.LELCH'     
!    end if
     flush(100_4)
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
              write(100,70) 'vector per node:',itise,varna_pos(1,istpp),adjustl(trim(filva))          
           else if(varna_pos(2,istpp)=='SCALA'.and.varna_pos(3,istpp)/='FILTE') then
              write(100,70) 'scalar per node:',itise,varna_pos(1,istpp),adjustl(trim(filva))       
           end if
        end do
        write(100,'(a)') 'TIME' 
        write(100, 80  ) 'time set:               ',itise
        write(100, 80  ) 'number of steps:        ',nppti_ens
        write(100,'(a)') 'filename start number:        1 '
        write(100,'(a)') 'filename increment:           1 '
        write(100,'(a)') 'time values: '
        write(100,'(10(1x,e12.5))') (tipoe_ens(i),i=1,nppti_ens)           
        flush(100_4)

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

end subroutine zfemres
