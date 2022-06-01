!***************************************************************
!*
!*   Module for input/output operations
!*
!***************************************************************
 MODULE InpOut
   use KindType
   IMPLICIT NONE
   SAVE
!
!***  File logical units
!
      integer(ip), parameter  :: lulog = 10     ! log   file
!
!***  Input file names
! 
     character(len=s_file) :: path,name,fname 
     character(len=s_file) :: finp,fpatch,fdom_bkg,fgeo_bkg,fbcs_bkg
     character(len=s_file) :: flog,fdat,fdom,fgeo,fgid,fdod,fdof,fker,fnsi,ftur,fgid_patch  
!
!***  Output files
!
      logical :: out_screen = .true.
      logical :: out_gid    = .true.
!
!***  List of Warnings 
!
      integer(ip), parameter :: maxwarn = 100
      integer(ip)            :: nwarn = 0
      character(len=s_mess)  :: warning(maxwarn) 
!
 CONTAINS
!
!
!
      subroutine get_input_path(fname,path)
!**************************************************************************
!*
!*    Gets the path of the input file
!*
!**************************************************************************
      implicit none
      character(len=*) :: fname,path
!
      integer(ip) :: ipos
      logical     :: go_on
!
      path(:) = ' ' 
      ipos  = LEN_TRIM(fname) + 1
      go_on = .true.
      do while(go_on)
         ipos = ipos - 1
         if(fname(ipos:ipos) == '/') then 
            go_on = .false.
            path(1:ipos) = fname(1:ipos)
         else if(ipos == 1) then
            go_on = .false. 
            path = './'
         end if
      end do
!
      return
      end subroutine get_input_path
!
!
!
      subroutine get_input_name(fname,name)
!**************************************************************************
!*
!*    Gets the path of the input file
!*
!**************************************************************************
      implicit none
      character(len=*) :: fname,name
!
      integer(ip) :: ipos,ilon
      logical     :: go_on
!
      name(:) = ' ' 
!
!     First remove the / (if any)
!
      ipos  = LEN_TRIM(fname) + 1
      ilon  = 0
      go_on = .true.
      do while(go_on)
         ipos = ipos - 1
         ilon = ilon + 1
         if(fname(ipos:ipos) == '/') then 
            go_on = .false.
            name(1:ilon) = fname(ipos+1:ipos+1+ilon)
         else if(ipos == 1) then
            go_on = .false. 
            name = fname
         end if
      end do
!
!     Then remove the point
!
      ipos  = LEN_TRIM(name) + 1  
      ilon  = LEN_TRIM(name)  
      go_on = .true.
      do while(go_on)
         ipos = ipos - 1
         if(name(ipos:ipos) == '.') then 
            go_on = .false.
            name(ipos:ilon) = ' '
         else if(ipos == 1) then
            go_on = .false. 
         end if
      end do
!
      return
      end subroutine get_input_name
!
!
!
      subroutine get_input_int(fname,block,line,value,nval,istat,message)
!**************************************************************************
!*
!*    Gets nval integer inputs from the file fname
!*
!*    INPUT:
!*    character*(*)   fname    Name of the file
!*    character*(*)   block    Block to search
!*    character*(*)   line     Line block to search
!*    integer         nval     Number of integers to read                     
!*
!*    OUTPUT:
!*    integer          istat    -1 ERROR  0 OK  1 WARNING
!*    integer          value    Values of the nval integers read
!*    character        message  Output message with the error description       
!*
!**************************************************************************
      implicit none
      character(len=*)      :: message,fname,block,line
      integer(ip)           :: nval,istat
      integer(ip)           :: value(nval)
!
      character(len=s_mess) :: mymessage
      character(len=s_long) :: card
      character(len=s_long) :: words(nwormax)
      logical               :: linefound,blockfound
      integer(ip)           :: ilen,nword,npar,ival,ipar
      real(rp)              :: param(nparmax),x0,xf,dx
!
!***  Initializations
!
      ilen = LEN(message)  
      message(:)  = ' '
      words(:)(:) = ' '
      istat = 0
!
!***  Opens the file
!
      open(90,FILE=TRIM(fname),STATUS='old',ERR=101)
!
!***  Search the line
!
      blockfound = .false.
      linefound  = .false.
      do while(.not.linefound)
         do while(.not.blockfound) 
             read(90,'(a256)',END=102) card
             call sdecode(card,words,param,nword,npar)
             if(words(1)(1:LEN_TRIM(block)).eq.block(1:LEN_TRIM(block))) blockfound=.true.
         end do
         read(90,'(a256)',END=103) card
         call sdecode(card,words,param,nword,npar)
         if(words(1)(1:LEN_TRIM(line)).eq.line(1:LEN_TRIM(line))) linefound = .true.
      end do
!
!***  Line format  FROM x0 TO xf INCREMENT dx
!***  Calculate npar
!
      if(TRIM(words(2)).eq.'FROM'.or.TRIM(words(2)).eq.'from') then
         x0 = param(1)
         xf = param(2)
         if(x0.gt.xf) goto 104
         dx = min(xf-x0,param(3))
         npar = INT((xf-x0)/dx)+1
         if(npar.gt.nparmax) then
            npar = nparmax
            istat = 1
            mymessage = 'get_input_int: warning. Too big number of parameters' 
            message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
         end if
         do ipar = 1,npar
           param(ipar) = x0 + (ipar-1)*dx
        end do
      end if
!                 
      if(npar.lt.nval) goto 105
!
      do ival = 1,nval
       value(ival) = INT(param(ival))
      end do
!
!***  Successful end
! 
      close(90)
      return
!
!***  List of errors
!
  101 istat = -1
      close(90)
      mymessage = 'get_input_int: error opening the input file '//TRIM(fname)
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  102 istat = -1
      close(90)
      mymessage = 'get_input_int: block '//TRIM(block)//' not found in the input file' 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  103 istat = -1
      close(90)
      mymessage = 'get_input_int: line '//TRIM(line)//' not found in the input file' 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  104 istat = -1
      close(90)
      mymessage = 'get_input_int: error in line '//line(1:LEN_TRIM(line)) 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  105 istat = -1
      close(90)
      mymessage = 'get_input_int: too few parameters in line '//TRIM(line) 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
!
      end subroutine get_input_int
!
!
!
      subroutine get_input_npar(fname,block,line,npar,istat,message)
!**************************************************************************
!*
!*    Gets the number of parameters associated to a line
!*
!*    INPUT:
!*    character*(*)   fname    Name of the file
!*    character*(*)   block    Block to search
!*    character*(*)   line     Line block to search
!*
!*    OUTPUT:
!*    integer         npar     Number of parameters                     
!*    integer         istat    -1 ERROR  0 OK  1 WARNING
!*    character       message  Output message with the error description       
!*
!**************************************************************************
      implicit none
!
      character(len=*) :: message,fname,block,line
      integer(ip)      :: npar,istat
!
      character(len=s_mess) :: mymessage
      character(len=s_name) :: card
      character(len=s_name) :: words(nwormax)      
      logical               :: linefound,blockfound
      integer(ip)           :: ilen,nword
      real(rp)              :: param(nparmax),x0,xf,dx
!
!***  Initializations
!
      ilen = LEN(message)  
      message(:)  = ' '
      words(:)(:) = ' '
      istat = 0
!
!***  Opens the file
!
      open(90,FILE=TRIM(fname),STATUS='old',ERR=101)
!
!***  Search the line
!
      blockfound = .false.
      linefound  = .false.
      do while(.not.linefound)
         do while(.not.blockfound) 
             read(90,'(a256)',END=102) card
             call sdecode(card,words,param,nword,npar)
             if(words(1)(1:LEN_TRIM(block)).eq.block(1:LEN_TRIM(block))) blockfound=.true.
         end do
         read(90,'(a256)',END=103) card
         call sdecode(card,words,param,nword,npar)
         if(words(1)(1:LEN_TRIM(line)).eq.line(1:LEN_TRIM(line))) linefound = .true.
      end do
!
!***  Line format  FROM x0 TO xf INCREMENT dx
!***  Calculate npar
!
      if(TRIM(words(2)).eq.'FROM'.or.TRIM(words(2)).eq.'from') then
         x0 = param(1)
             xf = param(2)
             if(x0.gt.xf) goto 104
             dx = min(xf-x0,param(3))
         npar = INT((xf-x0)/dx)+1
         if(npar.gt.nparmax) then
            npar = nparmax
            istat = 1
            mymessage = 'get_input_int4: warning. Too big number of parameters' 
            message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
         end if
      end if
!
!***  Successful end
!
      close(90)
      return
!
!***  List of errors
!
  101 istat = -1
      close(90)
      mymessage = 'get_input_npar : error opening the input file '//TRIM(fname)
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  102 istat = -1
      close(90)
      mymessage = 'get_input_npar : block '//TRIM(block)//' not found in the input file' 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  103 istat = -1
      close(90)
      mymessage = 'get_input_npar : line '//TRIM(line)//' not found in the input file' 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  104 istat = -1
      close(90)
      mymessage = 'get_input_npar: error in line '//line(1:LEN_TRIM(line)) 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
!
      end subroutine get_input_npar
!
!
!
      subroutine get_input_rea(fname,block,line,value,nval,istat,message)
!**************************************************************************
!*
!*    Gets nval real inputs from the file fname
!*
!*    INPUT:
!*    character*(*)   fname    Name of the file
!*    character*(*)   block    Block to search
!*    character*(*)   line     Line block to search
!*    integer         nval     Number of integers to read                     
!*
!*    OUTPUT:
!*    integer          istat    -1 ERROR  0 OK  1 WARNING
!*    real            value    Values of the nval integers read
!*    character       message  Output message with the error description       
!*
!**************************************************************************
      implicit none
      character(len=*)      :: message,fname,block,line
      integer(ip)           :: nval,istat
      real(rp)              :: value(nval)
!
      character(len=s_mess) :: mymessage
      character(len=s_name) :: card
      character(len=s_name) :: words(nwormax)
      logical               :: linefound,blockfound
      integer(ip)           :: ilen,nword,npar,ival,ipar
      real(rp)              :: param(nparmax),x0,xf,dx
!
!***  Initializations
!
      ilen = LEN(message)  
      message(:)  = ' '
      words(:)(:) = ' '
      istat = 0
!
!***  Opens the file
!
      open(90,FILE=TRIM(fname),STATUS='old',ERR=101)
!
!***  Search the line
!
      blockfound = .false.
      linefound  = .false.
      do while(.not.linefound)
         do while(.not.blockfound) 
             read(90,'(a256)',END=102) card
             call sdecode(card,words,param,nword,npar)
             if(words(1)(1:LEN_TRIM(block)).eq.block(1:LEN_TRIM(block))) blockfound=.true.
         end do
         read(90,'(a256)',END=103) card
         call sdecode(card,words,param,nword,npar)
         if(words(1)(1:LEN_TRIM(line)).eq.line(1:LEN_TRIM(line))) linefound = .true.
      end do
!
!***  Line format  FROM x0 TO xf INCREMENT dx
!***  Calculate npar
!
      if(TRIM(words(2)).eq.'FROM'.or.TRIM(words(2)).eq.'from') then
         x0 = param(1)
             xf = param(2)
             if(x0.gt.xf) goto 104
             dx = min(xf-x0,param(3))
         npar = INT((xf-x0)/dx)+1
         if(npar.gt.nparmax) then
            npar = nparmax
            istat = 1
            mymessage = 'get_input_rea: warning. Too big number of parameters' 
            message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
         end if
             do ipar = 1,npar
                param(ipar) = x0 + (ipar-1)*dx
             end do
      end if
!                
      if(npar.lt.nval) goto 105
!
      do ival = 1,nval
           value(ival) = param(ival)
        end do
!
!***  Successful end
!        
      close(90)
      return
!
  101 istat = -1
      close(90)
      mymessage = 'get_input_rea: error opening the input file '//TRIM(fname)
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  102 istat = -1
      close(90)
      mymessage = 'get_input_rea: block '//TRIM(block)//' not found in the input file' 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  103 istat = -1
      close(90)  
      mymessage = 'get_input_rea: line '//TRIM(line)//' not found in the input file' 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  104 istat = -1
      close(90)
      mymessage = 'get_input_rea: error in line '//line(1:LEN_TRIM(line)) 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  105 istat = -1
      close(90)
      mymessage = 'get_input_rea: too few parameters in line '//TRIM(line) 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
      end subroutine get_input_rea
!
!
!
      subroutine get_input_cha(fname,block,line,value,nval,istat,message)
!**************************************************************************
!*
!*    Gets nval character inputs from the file fname
!*
!*    NOTE: words are converted to UPPER case
!*
!*    INPUT:
!*    character*(*)   fname    Name of the file
!*    character*(*)   block    Block to search
!*    character*(*)   line     Line block to search
!*    integer         nval     Number of integers to read                     
!*
!*    OUTPUT:
!*    integer          istat    -1 ERROR  0 OK  1 WARNING
!*    character        value    Values of the nval integers read
!*    character        message  Output message with the error description       
!*
!*    CALLED BY: user
!*
!**************************************************************************
      implicit none
      integer(ip)           :: nval,istat
      character(len=*)      :: message,fname,block,line
      character(len=*)      :: value(nval)
!
      character(len=s_mess) :: mymessage
      character(len=s_name) :: card
      character(len=s_name) :: words(nwormax)
      logical               :: linefound,blockfound
      integer(ip)           :: ilen,nword,npar,ival,ipar,j
      real(rp)              :: param(nparmax)
!
!***  Initializations
!
      ilen = LEN(message)  
      message(:)  = ' '
      istat = 0
!
!***  Opens the file
!
      open(90,FILE=TRIM(fname),STATUS='old',ERR=101)
!
!***  Search the line
!
      blockfound = .false.
      linefound  = .false.
      do while(.not.linefound)
         do while(.not.blockfound) 
             read(90,'(a256)',END=102) card
             call sdecode(card,words,param,nword,npar)
             if(words(1)(1:LEN_TRIM(block)).eq.block(1:LEN_TRIM(block))) blockfound=.true.
         end do
         read(90,'(a256)',END=103) card
         call sdecode(card,words,param,nword,npar)
         if(words(1)(1:LEN_TRIM(line)).eq.line(1:LEN_TRIM(line))) linefound = .true.
      end do
!  
      if((nword-1).lt.nval) goto 104
!
      do ival = 1,nval
           value(ival)(1:LEN_TRIM(words(ival+1))) = words(ival+1)(1:LEN_TRIM(words(ival+1)))
!
!***     Fill the rest with ' '
!
           do j=LEN_TRIM(words(ival+1))+1,LEN(value(ival))
              value(ival)(j:j)=' '
           end do
        end do
!
!***  Convert to upper case
!
!      do ival = 1,nval
!           call upcase(value(ival))
!       end do
!
!***  Successful end
!        
      close(90)
      return
!
!***  List of errors
!
  101 istat = -1
      close(90)
      mymessage = 'get_input_cha: error opening the input file '//TRIM(fname)
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  102 istat = -1
      close(90)
      mymessage = 'get_input_cha: block '//TRIM(block)//' not found in the input file' 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  103 istat = -1
      close(90)
      mymessage = 'get_input_cha: line '//TRIM(line)//' not found in the input file' 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
  104 istat = -1
      close(90)
      mymessage = 'get_input_cha: too few parameters in line '//TRIM(line) 
      message(1:MIN(ilen,s_mess)) = mymessage(1:MIN(ilen,s_mess))
      return
!
      end subroutine get_input_cha
!
!
!
      subroutine sdecode(card,words,param,nword,npar)
!********************************************************************
!*
!*    This routine decodes a string card(s_name) into words and parameters 
!*
!********************************************************************
      implicit none
      integer(ip)  ::    nword,npar
      character(len=s_name) ::  card
      character(len=s_name) ::  words(nwormax)
      character(len=1)      ::  sstring(s_name)
      real(rp)              ::  param(nparmax)
!
      integer(ip) ::  ipos,first,last,nstr,lflag,i,ii,iword
      real(rp)    ::  digit
!
!***  Initializations
!
      nword = 0
      npar  = 0     
      ipos  = 0
      do while(1.ne.0)
!                                    ! First position
         ipos = ipos + 1
         if(ipos.gt.s_name) return 
  10     if(card(ipos:ipos).eq.' '.or.card(ipos:ipos).eq.'=') then
                ipos = ipos + 1 
            if(ipos.gt.s_name) return 
            go to 10
         end if
         first = ipos
!
         ipos = ipos + 1
         if(ipos.gt.s_name) return 
  20     if(card(ipos:ipos).ne.' '.and.card(ipos:ipos).ne.'=') then
                ipos = ipos + 1 
            if(ipos.gt.s_name) return 
            go to 20
         end if
         last = ipos-1
!
         nstr = last-first+1
         
         ii = 0
         do i=first,last
            ii = ii + 1
            sstring(ii) = card(i:i)
         end do
         call decod1(sstring,nstr,lflag,digit)
             if(lflag.eq.0) then
            npar = npar + 1
                    param(npar)= digit
             else if(lflag.eq.1) then
                nword = nword + 1
            words(nword)(:) = ' '
                words(nword)(1:nstr) = card(first:last)
         end if
!
      end do
      return
      end subroutine sdecode
!
!
!
      subroutine upcase(word)
!*********************************************************************** 
!*
!*    This routine converts word to upper case 
!*
!*********************************************************************** 
      implicit none
      character(len=*) :: word
      integer(ip)      :: iposi,ioctv,item1,item2,item3
      integer(ip)      :: ilen
!
      item1='141'o
      item2='172'o
      item3='40'o 
!
      ilen=LEN_TRIM(word)
!
      do iposi=1,ilen                                     ! process all positions
        ioctv=ichar(word(iposi:iposi))                    ! octal value
        if(item1.le.ioctv.and.item2.ge.ioctv) then        ! it is a lower case
             ioctv=ioctv-item3                            ! equivalent upper case
             word(iposi:iposi)=char(ioctv)                ! convert it to upcase
        end if
      end do ! iposi=1,ilen
!
      end subroutine upcase      
!
!
!
      subroutine decod1(string,nstr,lflag,digit)
!*******************************************************************
!*
!*    This subroutine decodes a single string(1:nstr) 
!*    
!*    If string(1:nstr) is a string returns lflag = 1
!*    If string(1:nstr) is a number returns lflag = 0 and the digit
!*
!*******************************************************************
      implicit none
      integer(ip) :: nstr,lflag
      integer(ip) :: istr,decim
      real(rp)    :: digit
      character(len=1) :: string(s_name)
!         
      lflag = 0                                             ! Number by default
      istr  = 1
!      
      do while(istr.le.nstr)
        decim = ichar(string(istr))                         ! decimal value
        if(decim.lt.48.or.decim.gt.57) then                 ! It is not a num.
           if(decim.ne.43 .and.decim.ne.45 .and.  &          ! discard + -
              decim.ne.68 .and.decim.ne.69 .and.  &          ! discard D E
              decim.ne.100.and.decim.ne.101.and.  &          ! discard d e
              decim.ne.46) then                             ! discard .
              lflag = 1
              istr  = nstr
          end if
        end if
        istr = istr+1
      end do
!
      if(lflag.eq.0) digit = stof(string,nstr)              ! It's a number
!      
      return
      end subroutine decod1
!
!
!
      real(rp) function stof(string,nstr)
!**************************************************************
!*
!*    This routine converts a real/integer number stored in a
!*    string(1:nstr) into a real(rp)  digit format
!*
!**************************************************************
      implicit none
      integer(ip) ::   nstr
      character(len=1) :: string(*)
!
      integer(ip) :: i,ipos,nsign,esign,nvalu
      integer(ip) :: expo,valu(s_name)  
      logical     :: next
!
      stof = 0.0_rp 
!
!***  Sing decoding
!
      ipos = 1
      if(ichar(string(ipos)).eq.43) then         !  + sign
        nsign = 1
            ipos  = ipos + 1
      else if(ichar(string(ipos)).eq.45) then    !  - sign
        nsign = -1
            ipos  = ipos + 1
          else                                       !  no sing (+)
            nsign = 1
            ipos  = ipos 
          end if
!
!***  Base decoding
!
      nvalu = 0
      next  = .true.
      do while(next)
        if((ichar(string(ipos)).eq.68 ).or. &       ! D
           (ichar(string(ipos)).eq.69 ).or. &       ! E
           (ichar(string(ipos)).eq.100).or. &       ! d
           (ichar(string(ipos)).eq.101).or. &       ! e
           (ichar(string(ipos)).eq.46 )) then       ! .
           next = .false.
            else
              nvalu = nvalu + 1
              valu(nvalu) = stof1(string(ipos))
              ipos = ipos + 1
           if(ipos.eq.(nstr+1)) then
                next = .false.
            ipos = ipos - 1
               end if
            end if
          end do
      do i = 1,nvalu
           stof = stof + valu(i)*1d1**(nvalu-i)
      end do
!
!***  Decimal decoding
!
      if((ichar(string(ipos)).eq.46   ).and.  &
                       ipos  .ne.nstr) then
        ipos = ipos + 1
        nvalu = 0
        next  = .true.
        do while(next)
          if((ichar(string(ipos)).eq.68 ).or. &        ! D
             (ichar(string(ipos)).eq.69 ).or. &       ! E
             (ichar(string(ipos)).eq.100).or. &       ! d
             (ichar(string(ipos)).eq.101)) then      ! e
             next = .false.
              else
                 nvalu = nvalu + 1
                 valu(nvalu) = stof1(string(ipos))
                 ipos = ipos + 1
             if(ipos.eq.(nstr+1)) then
                  next = .false.
              ipos = ipos - 1
                 end if
              end if
            end do
            do i = 1,nvalu
               stof = stof + valu(i)*1d1**(-i)
        end do
          end if
!
!***  Exponent
!
      if(((ichar(string(ipos)).eq.68 ).or. &        ! D
          (ichar(string(ipos)).eq.69 ).or. &        ! E
          (ichar(string(ipos)).eq.100).or. &        ! d
          (ichar(string(ipos)).eq.101)).and. &      ! e
                        ipos  .ne.nstr) then
        ipos = ipos + 1
        if(ichar(string(ipos)).eq.43) then         !  + sign
           esign = 1
               ipos  = ipos + 1
        else if(ichar(string(ipos)).eq.45) then    !  - sign
           esign = -1
               ipos  = ipos + 1
            else                                       !  no sing (+)
               esign = 1
               ipos  = ipos 
            end if     
!        
        nvalu = 0
        next  = .true.
        do while(next)
           nvalu = nvalu + 1
           valu(nvalu) = stof1(string(ipos))
               ipos = ipos + 1
           if(ipos.eq.(nstr+1)) then
                  next = .false.
              ipos = ipos - 1
               end if
            end do
            expo = 0
            do i = 1,nvalu
             expo = expo + valu(i)*10**(nvalu-i)
        end do
!
        if(esign.eq.1) then
           stof = stof*(10.0_rp**expo)
        else if(esign.eq.-1) then
           stof = stof/(10.0_rp**expo)
        end if
!
         end if
!
       stof = nsign*stof
       end function stof
!
!
!
     integer(ip) function stof1(string1)
!**************************************************************
!*
!*    Decodes a character*1 string
!*
!**************************************************************
     implicit none
     character(len=1) :: string1
!      
        if(string1.eq.'0') then
          stof1 = 0
        else if(string1.eq.'1') then
          stof1 = 1      
        else if(string1.eq.'2') then
          stof1 = 2      
        else if(string1.eq.'3') then
          stof1 = 3      
        else if(string1.eq.'4') then
          stof1 = 4      
        else if(string1.eq.'5') then
          stof1 = 5      
        else if(string1.eq.'6') then
          stof1 = 6      
        else if(string1.eq.'7') then
          stof1 = 7      
        else if(string1.eq.'8') then
          stof1 = 8      
        else if(string1.eq.'9') then
          stof1 = 9
       else
          stof1= 0
!          print *, 'error stof1, inpout.f90, string1 :', string1
        end if
      return
      end function stof1
!
 END MODULE InpOut 
