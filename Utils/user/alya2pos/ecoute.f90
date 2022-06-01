subroutine ecoute(subna)
  !-----------------------------------------------------------------------
  !
  ! Reads a string and interprets it as words and parameters.
  !
  !   - Maximum number of words and parameters = maxwp.
  !   - Only the first five characters of each word are decoded.
  !   - The underline characters '_' are discarted.
  !   - Lower case letters are converted to upper case.
  !   - Each word or parameter must be separated by ' ', '=', ':' or ','
  !   - A "comment" begins with '$', '!', '/'  or '/' in any place.
  !   - A line that has a comment beginning with '/' or '/'
  !     continuates in the next line.
  !   - A line beginning with title is not decoded. title directive.
  !   - A line beginning with include is an include directive.
  !   - A line beginning with echo turns 'on' or 'off' the echo.
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip,rp,fil_pdata
  use def_inpout
  implicit none
  !
  ! Var
  !
  character(*)     :: subna
  real(rp)          :: digit
  integer(ip)       :: first,firsp,i,last,lastp,ptrwo,npptr,nwptr
  integer(ip)       :: leng,flag,resum
  logical(lg)       :: newline=.false.
  logical(lg), save :: echo=.false. ! default echo off. to change it use: echo on
  !
  ! Data
  !      
  if(lispa==0) then
     nunit = lisda                                  !  initial data file
     lispa = 1
  end if
  !
  ! Begin
  !
  ccard=' '
  nnwor=0                                           ! initialize.
  nnpar=0
  nwptr=0
  npptr=0
  resum=0
  do i=1,maxwp
     words(i)=' '
     param(i)=0.0_rp
  end do
100 continue
  do while(((nnwor==0).and.(nnpar==0)&              ! don't return without answer
       .or.newline).or.resum==1)                    ! continue reading if / or \

     if (resum==0) then
        newline=.false.                             ! initialize.
        last=0
        lastp=0
     end if
     firsp=1
     resum=0
     read(nunit,10,end=101,err=1) ccard             ! read a card
     !     leng=lnbln1(ccard)                             ! calculate the length.
     leng=len_trim(ccard)                             ! calculate the length.

     decode_card: do while(last<leng)               ! decode all the card.
        first=last+1
        loop_first: do while(            &
             ccard(first:first)=='_'.or. &          ! jump null character (_)
             ccard(first:first)==' '.or. &          ! jump separators ( =:,)
             ccard(first:first)=='='.or. &
             ccard(first:first)==':'.or. &
             ccard(first:first)==','     )
           first=first+1
           if(first>leng) exit loop_first
        end do loop_first
        if(last==0) firsp=first                     ! save first to print card
        last=first
        loop_last: do while(             &
             ccard(last:last)/=' ' .and. &          ! look for separator ( =:,).
             ccard(last:last)/='=' .and. &
             ccard(last:last)/=':' .and. &
             ccard(last:last)/=',' .and. &
             ccard(last:last)/='$' .and. &          ! look for coment ($!/\).
             ccard(last:last)/='!' .and. &
             ccard(last:last)/='\' .and. &
             ccard(last:last)/='\'      )
           last=last+1
           if(last>leng) exit loop_last
        end do loop_last
        if(last<=251            .and.(&
             ccard(last:last)=='$'.or.&
             ccard(last:last)=='!')) leng=last-1
        if(last<=251            .and.(&
             ccard(last:last)=='\'.or.&
             ccard(last:last)=='\')) then
           leng=last-1                              ! deal with continuation
           newline=.true.                           ! set new line flag.
        end if
        last=last-1
        if(last>=first) then                        ! is it a word or a parameter
           lastp=last                               ! save last to print card
           call decod1(last-first+1_ip,&
                ccard(first:last),flag,digit)
           wname= ccard(first:last)                 ! keep the useful ccard part in wname
           if(flag==0) then                         ! it is a parameter.
              nnpar=nnpar+1                         ! # of parameters
              npptr=npptr+1                         ! integer :: to next parameter
              if(npptr>maxwp) go to 4               ! error.
              if(nwptr>npptr) npptr=nwptr
              param(npptr)=digit
           else                                     ! it is a word.
              nnwor=nnwor+1                         ! # of words
              nwptr=nwptr+1                         ! integer :: to next word
              if(nwptr>maxwp) go to 5               ! error.
              if(npptr>=nwptr) nwptr=npptr+1
              ptrwo=1
              do while ((first<=last).and.(ptrwo<=5))
                 words(nwptr)(ptrwo:ptrwo)=ccard(first:first)
                 ptrwo=ptrwo+1
                 first=first+1
                 do while (ccard(first:first)=='_') ! jump null character
                    first=first+1
                 end do
              end do
              call upcase(words(nwptr))             ! convert to upper case.
           end if
        end if ! (last>=first)
        if((nnwor==1).and.(nnpar==0).and.&          ! deal with title or include
             ((words(1)=='TITLE').or.((words(1)=='INCLU').and.&
             (subna/='NOREAD')))) then
           if(echo.and.(subna/='noecho'))&
                write(lisre,20) 'ecoute',ccard(firsp:leng)
           last=last+2
           do while(ccard(last:last)==' ')          ! remove blank spaces
              last=last+1
           end do
           if(leng<last) go to 6                    ! error
           ccard=ccard(last:leng)                   ! remove words(1) from ccard
           leng=leng-last+1
           if(words(1)=='TITLE') then               ! deal with titles
              firsp=1
              lastp=leng
           else                                     ! deal with include directive
              !if(nunit==lisin) go to 3              ! error
              last=1                                ! remove tail comments
              do while((last<=leng).and.&           ! look for end (last=leng)
                   ccard(last:last)/=' ')           ! look for separator ( )
                 last=last+1
              end do
              ccard=ccard(1:last-1)                 ! remove tail comments
              if(nunit==lisin) then
                 nunit=lisi1
              else
                 nunit=lisin
              end if
              call opincl()
              lastp=0                               ! to ignore the echo
              nnwor=0                               ! forget all
              nwptr=0
              words(1)=' '
           end if
           last=leng                                ! to break the do while
        else if(subna=='NOREAD') then
           newline=.false.
           nnwor=1
        end if
     end do decode_card

     if((words(1)=='ECHO') .and.&
          (nnpar==0).and.(nnwor==2)) then           ! deal with echo
        if(words(2)=='OFF') then
           echo=.false.
           if(subna/='noecho') write(lisre,20) 'ECOUTE','ECHO OFF'
        else
           echo=.true.
           if(subna/='noecho') write(lisre,20) 'ECOUTE','ECHO ON'
        endif
        nnwor=0                                     ! forget all
        nwptr=0
        do i=1,maxwp
           words(i)=' '
        end do
     else                                           ! print card
        if((echo).and.(firsp<=lastp).and.(subna/='noecho')) then
           if(newline) then
              lastp=lastp+2
              ccard(lastp-1:lastp)=' /'
           end if
           write(lisre,20) subna,ccard(firsp:lastp)
        end if
     end if

  end do ! while ((nnwor==0).and.(nnpar==0).or.newline)

  nwopa=max(npptr,nwptr)
  return
  !
  ! End of include
  !
101 continue 
  if(nunit/=lisin.and.nunit/=lisi1) then
     if(endst==1) then
        goto 2 ! error
     else
        words(1)='ENDFI'
        return
     end if
  end if
  close(unit=nunit,status='keep')
  if(nunit==lisi1) then
     nunit=lisin
  else
     nunit=lisda
  end if
  if(echo.and.(subna/='noecho'))&
       write(lisre,20) 'ECOUTE','END OF INCLUDE FILE'
  resum=1
  go to 100 ! for resume the error return to the same place.
  !
  ! Errors:
  !
1 call runend('ECOUTE: ERROR DETECTED WHEN READING')
2 call runend('ECOUTE: END OF FILE DETECTED IN SUBROUTINE '//trim(subna))
3 call runend('ECOUTE: ERROR: INCLUDE FROM INCLUDE')
4 call runend('ECOUTE: TOO MANY PARAM IN COMMAND  ')
5 call runend('ECOUTE: TOO MANY WORDS IN COMMAND  ')
6 call runend('ECOUTE: BLANK IS ILEGAL HERE       ')
7 call runend('ECOUTE: COULD NOT OPEN FILE: '//adjustl(trim(ccard)))
  !
  ! Format
  !
10 format(a250)
20 format(1x,a6,' <-- ',a)

end subroutine ecoute

subroutine opincl()
  !
  ! Include file cannot be opened: look the directory of the data file
  !
  use def_kintyp, only :  ip,fil_pdata,intost
  use def_inpout
  implicit none
  integer(ip)   :: i,istat
  integer(4)    :: istat4,nunit4
  character(20) :: wstat

  !inquire(unit=nunit,opened=lopen,exist=lexis,form=wform) 
  !print*,'popo=',trim(ccard)
  !print*,'caca=',nunit,lopen,lexis
  !print*,'pipi=',trim(wform)
  nunit4 = nunit
  open(unit=nunit4,file=adjustl(trim(ccard)),err=7,form='FORMATTED',status='OLD',IOSTAT=istat4)
  return
  
7 i=len(trim(fil_pdata))
  do while(i>1)
     i=i-1
     if(fil_pdata(i:i)=='/'.or.fil_pdata(i:i)==achar(92)) i=-i
  end do
  if(i<0) then
     i=-i
     open(unit=nunit,file=trim(fil_pdata(1:i))//trim(ccard),err=8,& 
          form='formatted',status='old')
     return
  else
     istat = int(istat4,ip)
     wstat = intost(istat)
     call runend('ECOUTE: COULD NOT OPEN FILE: '//trim(ccard)//'. ERROR '//trim(wstat))    
  end if
  
8 continue
  istat = int(istat4,ip)
  wstat = intost(istat)
  call runend('ECOUTE: COULD NOT OPEN FILE: '//adjustl(trim(fil_pdata(1:i))&
       //trim(ccard))//'. ERROR '//trim(wstat))

end subroutine opincl
