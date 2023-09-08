program test01
    use unitt
    implicit none
    call assert_true(.TRUE.)
    call assert_false(.FALSE.)
    !call assert_true(.FALSE.)
    call assert_false(.TRUE.)
end program test01
