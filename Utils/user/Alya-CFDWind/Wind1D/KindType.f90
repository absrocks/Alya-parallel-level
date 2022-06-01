!***************************************************************
!*
!*    Module for variable kind definition
!*
!***************************************************************
     MODULE KindType
     implicit none
     save
!
     integer    , parameter :: ip = 4
     integer(ip), parameter :: rp = 8
!
     integer(ip), parameter :: s_name = 256    ! generic string
     integer(ip), parameter :: s_file = 256    ! file
     integer(ip), parameter :: s_mess = 256    ! message
     integer(ip), parameter :: s_long = 256

!
!*** Number of words and parameters
! 
     integer(ip), parameter :: nwormax = 128 
     integer(ip), parameter :: nparmax = 128

     END MODULE KindType
