subroutine vu_res(&
     ittim,npoin,nelem,parr1,pari1,lun_asc,lun_bin,lexis,&
     lbxis,ltype,pdime,namda,wopos,scavec,intrea,wbytes,&
     rttim,kfl_markm,kfl_multi,kfl_bound) 


  use def_kintyp, only          :  ip,rp,cenam,nnode,cetop,cepos,nelty
  use def_kintyp, only          :  ltyp2,lllll,nllll,llll2,jttim
  use def_kintyp, only          :  ncoun_pos,intost
  use def_elmtyp

  implicit none
  integer(ip),    intent(in)    :: ittim,npoin,nelem,lun_asc,lun_bin,pdime
  integer(ip),    intent(in)    :: lexis(*)
  integer(ip),    intent(in)    :: lbxis(*)
  integer(ip),    intent(in)    :: ltype(*)
  real(rp),       intent(in)    :: rttim
  real(rp),       intent(in)    :: parr1(*)
  integer(ip),    intent(in)    :: pari1(*)
  character(150), intent(in)    :: namda
  character(5),   intent(in)    :: wopos
  character(5),   intent(in)    :: scavec
  character(5),   intent(in)    :: intrea
  character(5),   intent(in)    :: wbytes
  integer(ip),    intent(in)    :: kfl_markm
  integer(ip),    intent(in)    :: kfl_multi
  integer(ip),    intent(in)    :: kfl_bound
  character(8)                  :: chtim
  character(20)                 :: wopo2(3)
  integer(ip)                   :: ipoin,ielem,ielty,idime
  integer(ip)                   :: mpoin,nbyte,jelty
  character(150)                :: fil_bin,fil_asc
  character(150)                :: celem
  !
  ! File name
  !
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
  if( kfl_multi == 0 ) then
     fil_asc = trim(namda)//'-'//trim(chtim)//'.res.vu'
     fil_bin = trim(namda)//'-'//trim(chtim)//'.res.vubin'
     if( ittim /= jttim ) then
        ncoun_pos = 0
        if( jttim /= -1 ) then
           close(lun_asc)
           close(lun_bin)
        end if
        open(unit=lun_asc,file=trim(fil_asc),form='formatted') 
        open(unit=lun_bin,file=trim(fil_bin),form='unformatted') 
     end if
  else
     ncoun_pos = 0
     fil_asc = trim(namda)//'-'//trim(wopos)//'-'//trim(chtim)//'.res.vu'
     fil_bin = trim(namda)//'-'//trim(wopos)//'-'//trim(chtim)//'.res.vubin'
     open(unit=lun_asc,file=trim(fil_asc),form='formatted') 
     open(unit=lun_bin,file=trim(fil_bin),form='unformatted') 
  end if

  !
  ! Initialization
  !
  jttim = ittim
  mpoin = npoin
  if( wbytes == '2BYTE' ) then
     nbyte = 2
  else if( wbytes == '4BYTE' ) then
     nbyte = 4
  else if( wbytes == '8BYTE' ) then
     nbyte = 8
  end if

  write(lun_asc,205) rttim

  if( scavec == 'SCALA' ) then
     !
     ! SCALAR
     !
     ncoun_pos=ncoun_pos+4_ip 
     write(lun_asc,200) trim(wopos),trim(fil_bin),mpoin,ncoun_pos
     ncoun_pos=ncoun_pos+4_ip+nbyte*mpoin
     write(lun_asc,210) 
     if( kfl_markm >= 3 ) then
        do ielty=1,nllll
           if(lllll(ielty)/=0) then
              celem=intost(ielty)
              jelty=llll2(ielty)
              write(lun_asc,220) trim(wopos),&
                   trim(cepos(jelty)),trim(wopos),'Connec'//trim(celem),&
                   'Zone'//trim(cepos(jelty))//'_'//trim(celem)
           end if
        end do
     else
        do ielty=1,nelty
           if(lexis(ielty)/=0) then
              celem=intost(ielty)
              write(lun_asc,220) trim(wopos),&
                   trim(cepos(ielty)),trim(wopos),'Connec'//trim(celem),&
                   'Zone'//trim(cepos(ielty))//'_'//trim(celem)
           end if
        end do
     end if
     if (kfl_bound == 1) then
        do ielty=1,nelty
           if(lbxis(ielty)/=0) then
              celem=intost(ielty)
              write(lun_asc,220) trim(wopos),&
                trim(cepos(ielty)),trim(wopos),'Connecb'//trim(celem),&
                'Zone'//trim(cepos(ielty))//'_'//trim(celem)
           end if
        end do
     end if
     write(lun_asc,230)
     if( intrea == 'INTEG' ) then
        write(lun_bin) (pari1(ipoin),ipoin=1,mpoin)
     else
        write(lun_bin) (parr1(ipoin),ipoin=1,mpoin)
     end if

  else if( scavec == 'VECTO' ) then
     !
     ! VECTOR
     !
     wopo2(1) = trim(wopos) // '_X'
     wopo2(2) = trim(wopos) // '_Y'
     wopo2(3) = trim(wopos) // '_Z'
     do idime = 1,pdime
        ncoun_pos=ncoun_pos+4_ip 
        write(lun_asc,200) trim(wopo2(idime)),trim(fil_bin),mpoin,ncoun_pos
        ncoun_pos=ncoun_pos+4_ip+nbyte*mpoin
        write(lun_asc,210) 
        if( kfl_markm >= 3 ) then
           do ielty=1,nllll
              if(lllll(ielty)/=0) then
                 celem=intost(ielty)
                 jelty=llll2(ielty)
                 write(lun_asc,220) trim(wopo2(idime)),&
                      trim(cepos(jelty)),trim(wopo2(idime)),'Connec'//trim(celem),&
                      'Zone'//trim(cepos(jelty))//'_'//trim(celem)
              end if
           end do
        else
           do ielty=1,nelty
              if(lexis(ielty)/=0) then
                 celem=intost(ielty)
                 write(lun_asc,220) trim(wopo2(idime)),&
                      trim(cepos(ielty)),trim(wopo2(idime)),'Connec'//trim(celem),&
                      'Zone'//trim(cepos(ielty))//'_'//trim(celem)
              end if
           end do
        end if
        if (kfl_bound == 1) then
           do ielty=1,nelty
              if(lbxis(ielty)/=0) then
                 celem=intost(ielty)
                 write(lun_asc,220) trim(wopo2(idime)),&
                   trim(cepos(ielty)),trim(wopo2(idime)),'Connecb'//trim(celem),&
                   'Zone'//trim(cepos(ielty))//'_'//trim(celem)
              end if
           end do
        end if
        write(lun_asc,230)
        if( intrea == 'INTEG' ) then
           write(lun_bin) (pari1((ipoin-1)*pdime+idime),ipoin=1,mpoin)
        else
           write(lun_bin) (parr1((ipoin-1)*pdime+idime),ipoin=1,mpoin)
        end if
     end do

  end if

  if( kfl_multi == 1 ) then
     close(unit=lun_asc)
     close(unit=lun_bin)
  end if
  !
  ! Format
  !
200 format('FIELD<double> ',a,'("',a,'",',i12,',',i12,');')
201 format('FIELD<float>  ',a,'("',a,'",',i12,',',i12,');')
205 format('TEXTE Time(" ',e13.6,'");')
210 format('SOLUTION Solution( ) =',/,'{')
  !220 format('   VARIABLE ',a,'( ',a,',',a,',','Connec',i1,',Zone',i1,');')
220 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,');')
221 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,',',a,');')

230 format('};')

end subroutine vu_res


