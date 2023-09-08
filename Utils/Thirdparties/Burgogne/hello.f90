!Error handling
subroutine raiseError(err)
    INTEGER(4) :: err
    IF (err == 1001) STOP "Error: Your license is not valid!"
    IF (err == 1002) STOP "Error: Your license has expired!"
    IF (err == 1003) STOP "Error: Are you trying to cheat me?"
    STOP "Unkown ERROR!"
end subroutine


!Convert str to int
subroutine str2int(str,int,stat)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,intent(out)         :: stat

    read(str,*,iostat=stat)  int
end subroutine str2int

subroutine strlen(st, l)
      INTEGER(4) :: i
      character  :: st*(*)
      i = len(st)
      do while (st(i:i) < '' .OR. st(i:i) == ' ')
        !print *, i,": '",st(i:i),"'"
        i = i - 1
      enddo
      l = i
end subroutine strlen

!Check the license
!DEC$ ATTRIBUTES FORCEINLINE :: check
subroutine check(name)

    !External vars
    CHARACTER(LEN=500)   :: name

    !Local vars
    INTEGER(4) :: l = 0           ! File Name lenght.
    INTEGER(4) :: initDate = 0    ! Initial Date of valid license.
    INTEGER(4) :: endDate  = 0    ! Finall Date of valid license.
    INTEGER(4) :: actDate  = 0    ! Activation Date limit.
    INTEGER(4) :: lastUse  = 0    ! Last use of the license.
    INTEGER(4) :: maxHour  = 0    ! Max avaiable hours to execute Alya.
    INTEGER(4) :: hourCounter = 0 ! Actual counter of used hours.
    INTEGER(4) :: today           ! Today Date.
    INTEGER(4) :: stat            ! Aux var for converting str2int.
    CHARACTER(8) :: date          ! Today Data in string format.


    !Calling C function to decrypt the license.
    external licensetopass

    call strlen(name, l)
    call licensetopass(name, l, initDate, endDate, actDate, lastUse, maxHour, hourCounter)

    !Get the date of today
    call date_and_time(DATE=date)
    call str2int(date, today, stat)

    !Debug.
    !print *, "initDate is: ", initDate
    !print *, "endDate is:  ", endDate
    !print *, "actDate is:  ", actDate
    !print *, "lastUse is:  ", lastUse
    !print *, "today is:    ", today
    !print *, "maxHour is:  ", maxHour
    !print *, "hourCounter: ", hourCounter

    !Raising error if license is not valid.
    IF (today < initDate)THEN 
        call raiseError(1001)
    ELSE IF (today > endDate)THEN 
        call raiseError(1002)
    ELSE IF (actDate == -1 .AND. lastUse == -1)THEN
        call raiseError(1003)
    ELSE IF (today > actDate .AND. actDate /= -1) THEN
        call raiseError(1002)
    ELSE IF (today < lastUse)THEN
        call raiseError(1003)
    ELSE IF (maxHour <= hourCounter)THEN
        call raiseError(1002)
    END IF
end subroutine check

subroutine actualize(licName, newHourCounter)
    INTEGER(4)    :: l = 0           ! File Name lenght.
    REAL(8)       :: newHourCounter  ! Actualized Counter
    INTEGER(4)    :: today           ! Today Date.
    CHARACTER(8)  :: date            ! Today Data in string format.
    CHARACTER(500):: licName         ! License File
    INTEGER(4)    :: stat            ! Aux var for converting str2int.

    external cactualize

    call strlen(licName, l)

    call date_and_time(DATE=date)
    call str2int(date, today, stat)

    call cActualize(licName, l, newHourCounter, today)

end subroutine actualize


!program hello
!    CHARACTER(50)   :: s = "Alya-License.lic"
!    CHARACTER(50)   :: licenseName
!    INTEGER(4)      :: hourCounter = 50
!
!    CALL getarg(1, licenseName)
!
!    print *,""
!    print *,"Starting the best HelloWorld in the World!!!"
!    print *,""
!
!    !Check License
!    call check(licenseName)
!
!    print *, "Hello World ¬¬"
!    print *,""
!
!   !Actualize License
!    call actualize(s, hourCounter)
!
!    print *, "Bye"
!    print *,""
!
!end program Hello
