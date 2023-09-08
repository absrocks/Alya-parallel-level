module mod_optim
  use def_kintyp, only       :  ip,rp,lg
  !
  !     Data structure for optimization by swapping
  !
  !
  !     Indeces of the triangles in the plane
  !
  integer(ip)     :: indx(3,4,14,6)=RESHAPE((/&
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
                                ! 2 sides
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
                                ! 3 Sides, 1 triangle,  1 possibility
       1,2,3,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
                                ! 4 sides, 2 triangles, 2 posibilities 
       1,2,3,  1,3,4,  0,0,0, 0,0,0, &
       2,3,4,  2,4,1,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
                                ! 5 sides, 3 triangles, 5 possibilities
       1,2,3,  1,4,5,  1,3,4, 0,0,0, &
       2,3,4,  1,4,5,  1,2,4, 0,0,0, &
       1,2,5,  2,3,4,  2,4,5, 0,0,0, &
       3,4,5,  1,2,3,  1,3,5, 0,0,0, &
       1,2,5,  3,4,5,  2,3,5, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
       0,0,0,  0,0,0,  0,0,0, 0,0,0, &
                                ! 6 sides, 4 triangles, 14 possibilities
       1,2,3,  1,5,6,  1,3,4,  1,4,5,&                          
       2,3,4,  1,2,6,  2,4,5,  2,5,6,&
       3,4,5,  1,2,3,  1,3,6,  6,3,5,&
       4,5,6,  2,3,4,  1,2,4,  1,4,6,&
       1,5,6,  3,4,5,  1,2,5,  2,3,5,&
       1,2,6,  6,4,5,  6,2,3,  6,3,4,&
       3,4,5,  1,5,6,  1,2,3,  1,3,5,&
       1,2,3,  4,5,6,  1,3,4,  1,4,6,&
       2,3,4,  1,5,6,  1,2,4,  1,4,5,&
       1,2,6,  4,5,6,  2,3,4,  2,4,6,& 
       1,2,6,  3,4,5,  2,3,5,  2,5,6,&
       2,3,4,  1,5,6,  1,2,5,  2,4,5,&  
       1,2,3,  4,5,6,  4,6,3,  1,3,6,&  
       1,2,6,  3,4,5,  6,2,3,  6,3,5 /),(/3,4,14,6/))
  !
  !     Position of the repetition
  ! 
  integer(ip)           :: pos(56,6)=RESHAPE((/&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,&
                                !2 sides 
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,&
                                !3 sides
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,&
                                !4 sides
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,&
                                !5 sides
       11,5,0,8,0,0,13,0,0,14,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0,&
                                !6 sides
       10,17,31,36,14,21,48,44,18,27,52,56,&
       22,33,35,32,26,25,47,43,37,30,55,51,&
       42,34,29,0,49,38,0,0,39,46,0,0,&
       41,50,45,0,53,54,0,0,0,0,0,0,&
       0,0,0,0,0,0,0,0 /),(/56,6/))
  !
  !     External neighboring
  !
  integer(ip)       ::  exte(2,6,14,6)=RESHAPE((/&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
                                !2 sides 
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
                                !3 sides 
       1,3,1,1,1,  2,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
                                !4 sides 
       1,3,1,1,2,  1,2,2,0,0,0,0,&
       2,2,1,3,1,  1,2,1,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
                                !5 sides
       1,3,1,1,3,  1,2,1,2,2,0,0,& 
       3,3,1,3,1,  1,2,1,2,2,0,0,&  
       1,3,2,3,2,  1,3,1,1,2,0,0,&
       2,3,2,1,1,  3,1,1,3,2,0,0,&   
       1,3,3,3,2,  3,2,1,1,2,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
       0,0,0,0,0,  0,0,0,0,0,0,0,&
                                !6 sides
       1,3,1,1,3,1,  4,1,2,1,2,2,&
       2,3,1,3,1,1,  3,1,4,1,2,2,&
       2,3,2,1,1,3,  1,1,4,2,3,2,& 
       3,3,2,3,2,1,  1,3,1,1,4,2,&
       3,3,4,3,2,3,  2,1,1,1,1,2,&
       1,3,3,1,4,1,  2,1,2,2,1,2,&
       3,3,3,1,1,3,  1,1,2,1,2,2,&
       1,3,1,1,3,1,  2,3,2,1,4,2,&
       3,3,1,3,1,1,  4,1,2,1,2,2,&
       1,3,3,3,3,1,  2,3,2,1,1,2,&
       1,3,3,3,2,3,  2,1,4,1,1,2,&
       3,3,1,3,1,1,  4,1,2,1,2,2,&
       1,3,1,1,3,2,  2,3,2,1,4,2,&
       1,3,3,1,2,3,  2,1,4,2,1,2   /),(/2,6,14,6/))
  !
  !     Internal neighboring 1
  !
  integer(ip)       ::  inte1(2,3,14,6)=RESHAPE((/&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 2 sides 
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 3 sides 
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 4 sides 
       1,2,0,0,0,0,&
       1,2,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 5 sides 
       1,2,3,2,0,0,&
       1,2,3,2,0,0,&
       1,1,3,3,0,0,&
       1,2,3,3,0,0,&
       1,1,3,1,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 6 sides 
       1,2,3,2,4,2,&
       1,2,3,2,4,2,&
       1,2,4,3,3,3,&
       1,2,4,3,3,1,&
       1,3,3,1,4,1,&
       1,1,3,2,4,2,&
       1,2,4,2,4,3,&
       1,2,3,2,4,1,&
       1,2,3,2,4,2,&
       1,1,4,1,4,3,&
       1,1,4,3,3,1,&
       1,2,4,2,3,2,&
       1,2,4,1,3,3,&
       1,1,3,2,4,1/),(/2,3,14,6/))
  !
  !     Internal neighboring 2
  !
  integer(ip)       ::  inte2(2,3,14,6)=RESHAPE((/&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 2 sides 
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 3 sides 
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 4 sides 
       2,3,0,0,0,0,&
       2,3,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 5 sides 
       3,3,2,3,0,0,&
       3,1,2,3,0,0,&
       3,2,2,2,0,0,&
       3,1,2,2,0,0,&
       3,2,2,2,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 6 sides 
       3,3,4,3,2,3,&
       3,3,4,3,2,1,&
       4,1,3,1,2,2,&
       4,1,3,2,2,2,&
       3,2,4,2,2,2,&
       3,3,4,3,2,3,&
       4,1,2,3,3,2,&
       3,3,4,3,2,2,&
       3,1,4,3,2,3,&
       4,2,2,2,3,2,&
       4,2,3,2,2,2,&
       4,3,3,1,2,3,&
       4,3,3,1,2,2,&
       3,3,4,3,2,2/),(/2,3,14,6/))
  !
  !     Internal edges 
  !
  integer(ip)       ::  indedg(2,3,14,6)=RESHAPE((/&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 2 sides 
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 3 sides 
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 4 sides 
       1,3,0,0,0,0,&
       2,4,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 5 sides 
       1,3,1,4,0,0,&
       1,4,2,4,0,0,&
       2,5,2,4,0,0,&
       1,3,5,3,0,0,&
       2,5,3,5,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
       0,0,0,0,0,0,&
                                ! 6 sides 
       1,3,1,4,1,5,&
       2,4,2,5,2,6,&
       3,5,3,6,3,1,&
       4,6,4,1,4,2,&
       5,1,5,2,5,3,&
       6,2,6,3,6,4,&
       1,5,1,3,5,3,&
       1,3,1,4,6,4,&
       2,4,1,4,1,5,&
       2,4,4,6,6,2,&
       2,5,2,6,3,5,&
       2,4,2,5,1,5,&
       1,3,6,3,6,4,&
       2,6,3,6,5,3/),(/2,3,14,6/))
  !
  !     Internal element 
  !
  integer(ip)       ::  indedglm(3,14,6)=RESHAPE((/&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
                                ! 2 sides 
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
                                ! 3 sides 
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
                                ! 4 sides 
       1,0,0,&
       1,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
                                ! 5 sides 
       1,3,0,&
       2,3,0,&
       1,3,0,&
       2,3,0,&
       1,3,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
       0,0,0,&
                                ! 6 sides 
       1,3,4,&
       3,4,2,&
       4,3,2,&
       4,3,2,&
       3,4,2,&
       3,4,2,&
       4,4,4,&
       3,4,2,&
       3,4,2,&
       4,4,4,&
       4,4,3,&
       4,4,3,&
       4,4,3,&
       3,3,4/),(/3,14,6/))
  !
  !     Number of combinations
  !
  integer(ip)           :: lcatn(6)=(/0,0,1,2,5,14/)
  !
  !     Number of triangles
  !
  integer(ip)           :: ltri(6)=(/0,0,1,2,3,4/)

contains

  subroutine optimsh()
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only    : memor_msh
    use mod_mshtol, only    : ptoelmb,tetoteb
    implicit none
    integer(ip)           :: ndim,npoin,nnode,nnosi,nnofa
    integer(ip)           :: nelem,ipoin,nedg,isto,inode,jpoin,ielem
    integer(ip)           :: ip1,ip2,ip3,ineigh,nhole
    integer(ip),pointer   :: elem(:,:),eltoel(:,:),ptoel1(:),ptoel2(:)
    integer(ip),pointer   :: lmark(:),lhole(:)
    real(rp),pointer      :: rsize(:),coor(:,:)
    real(rp)              :: rlen,c00,rx,ry,rz,rl
    integer(ip)           :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))  
    character*80          :: filename
    integer(4)                      :: istat
    !
    !     This sub optimizes a mesh by calling optim3d
    !  
    ndim=3_ip
    nnosi=2_ip
    nnofa=3_ip
    nnode=4_ip
    c00=0.0d+00
    filename="naca"
    !
    !     Read the mesh
    !
    call readoptim(nelem,npoin,nnode,ndim,elem,coor,filename)
    !
    !    Allocate rsize 
    !
    allocate(rsize(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'RSIZE','optimsh',rsize)
    allocate(lhole(max(1_ip,npoin/100_ip)),stat=istat)
    call memchk(zero,istat,memor_msh,'LHOLE','optimsh',lhole)
    !
    !     Get the elements surrounding the points
    !
    call ptoelmb(elem,nelem,npoin,nnode,ptoel1,ptoel2)
    !
    !    Allocate lmark 
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','optimsh',lmark)
    !
    !     Fill rsize
    ! 
    do ipoin=1,npoin
       lmark(ipoin)=ipoin
       nedg=0_ip
       rlen=c00
       do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(isto) 
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             if(lmark(jpoin)/=ipoin)then
                lmark(jpoin)=ipoin
                nedg=nedg+1_ip
                rx=coor(1,ipoin)-coor(1,jpoin) 
                ry=coor(2,ipoin)-coor(2,jpoin) 
                rz=coor(3,ipoin)-coor(3,jpoin) 
                rl=sqrt(rx*rx+ry*ry+rz*rz)
                rlen=rlen+rl
             endif
          enddo
       enddo
       rsize(ipoin)=rlen/real(nedg)
    enddo
    !
    !     Get the elements surrounding elements 
    !
    call tetoteb(elem,nnode,nelem,ptoel1,ptoel2,npoin,eltoel)
    !
    !    Deallocate ptoel1 and ptoel2 
    !
    call memchk(2_ip,istat,memor_msh,'PTOEL1','optimsh',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','optimsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','optimsh',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','optimsh',0_ip)
    !
    !     Write the mesh
    !
    !call writeoptim(nelem,npoin,nnode,ndim,elem,coor)
    !
    !     Clean up lmark
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Mark the boundary points
    !
    do ielem=1,nelem
       if(elem(1,ielem)==0)cycle
       do inode=1,nnode
          ineigh=eltoel(inode,ielem)
          if(ineigh==0)then
             ip1=elem(ltab(1,inode),ielem)   
             ip2=elem(ltab(2,inode),ielem)   
             ip3=elem(ltab(3,inode),ielem)   
             lmark(ip1)=1_ip 
             lmark(ip2)=1_ip 
             lmark(ip3)=1_ip 
          endif
       enddo
    enddo
    !
    !     Optimize the mesh
    !
    call optim3d(nelem,ndim,npoin,nnode,elem,coor,eltoel,nnosi,rsize,nnofa,lmark,lhole,nhole)
    !
    !     Write the mesh
    !
    !call writeoptim(nelem,npoin,nnode,ndim,elem,coor)

    call memchk(2_ip,istat,memor_msh,'LHOLE','optimsh',lhole)
    deallocate(lhole,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHOLE','optimsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'RSIZE','optimsh',rsize)
    deallocate(rsize,stat=istat)
    if(istat/=0) call memerr(2_ip,'RSIZE','optimsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELEM','optimsh',elem)
    deallocate(elem,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELEM','optimsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'ELTOEL','optimsh',eltoel)
    deallocate(eltoel,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELTOEL','optimsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'COOR','optimsh',coor)
    deallocate(coor,stat=istat)
    if(istat/=0) call memerr(2_ip,'COOR','optimsh',0_ip)

  end subroutine optimsh

  subroutine optim3d(nelem,ndim,npoin,nnode,elem,coor,eltoel,nnosi,rsize,&
       nnofa,lboup,lhole,nhole)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only    : memor_msh
    implicit none
    integer(ip),intent(in)          :: ndim,npoin,nnode,nnosi,nnofa
    integer(ip),intent(inout)       :: nelem,nhole
    integer(ip),pointer             :: elem(:,:),eltoel(:,:)
    real(rp),intent(in)             :: rsize(npoin)
    integer(ip),intent(inout)       :: lboup(npoin)
    real(rp),intent(inout)          :: coor(ndim,npoin)
    real(rp),pointer                :: rqual(:),rqualed(:)
    real(rp)                        :: c00,qmax,rvol,qtol,rq
    integer(ip)                     :: ielem,iter,niter,nedge,ihole,nelem0,inode 
    integer(ip)                     :: nedge0,iedge,ip1,ip2,ip3,ip4,iemax,ipoin,nbedg
    integer(ip),pointer             :: lrenu(:),ledge(:,:),ledglm(:),lptet(:)
    integer(ip),pointer             :: lhole(:),lhash1(:,:),lhash2(:),lbhash1(:),lbhash2(:)
    logical(lg),pointer             :: ledgact(:)
    integer(4)                      :: istat
    integer(ip)                     :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))  

    c00=0.0d+00
    niter=3_ip
    !
    !     This subroutine is responsible for improving the quality of a mesh 
    !
    !     On input: 
    !               - the mesh connectivity elem 
    !               - the coordinates array coor 
    !               - the face neighbors eltoel
    !
    !               - lptet carries the tet number for each point.
    !                 If the points is not optimized by a move, it becomes 0
    !                 It is set initially to zero for boundary points but
    !                 may be set to one due to inner edges touching boundary points   
    !               - ledglm is the edge to element pointer
    !
    !
    !
    !
    !
    !
    allocate(lptet(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPTET','optim3d',lptet)
    !
    !     Mark lptet from the elements to identify possible isolated points
    !
    do ielem=1,nelem
       if(elem(1,ielem)/=0)then
          lptet(elem(1,ielem))=1 
          lptet(elem(2,ielem))=1 
          lptet(elem(3,ielem))=1 
          lptet(elem(4,ielem))=1 
       endif
    enddo
    !
    !    Correct in lptet the boundary points and isolated points
    !
    do ipoin=1,npoin
       if(lboup(ipoin)==1)then
          lptet(ipoin)=0_ip
       else
          if(lptet(ipoin)==0)then
             lboup(ipoin)=1
           endif
       endif
    enddo
    !
    !     Check that the points belonging to the boundary faces are in lboup
    !
    do ielem=1,nelem
       if(elem(1,ielem)==0)cycle
       do inode=1,nnode
          if(eltoel(inode,ielem)==0)then
             ip1=elem(ltab(1,inode),ielem)
             if(lboup(ip1)==0)then
                write(*,*)'Warning, boundary point not marked'  
                lptet(ip1)=0_ip
                lboup(ip1)=1_ip
             endif
             ip1=elem(ltab(2,inode),ielem)
             if(lboup(ip1)==0)then
                write(*,*)'Warning, boundary point not marked'  
                lptet(ip1)=0_ip
                lboup(ip1)=1_ip
             endif
             ip1=elem(ltab(3,inode),ielem)
             if(lboup(ip1)==0)then
                write(*,*)'Warning, boundary point not marked'  
                lptet(ip1)=0_ip
                lboup(ip1)=1_ip
             endif
          endif
       enddo
    enddo
    !
    !     Allocate quality array
    !
    allocate(rqual(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'RQUAL','optim3d',rqual)
    !
    !     Set tolerance below which nothing will be done
    !
    qtol=2.0d+00
    !
    !     First compute the volume and quality of the old mesh
    !
    qmax=c00
    do ielem=1,nelem
       if(elem(1,ielem)==0)cycle
       call vol3d(elem,nelem,nnode,coor,ndim,npoin,ielem,rvol)
       call qual3d(elem,nelem,nnode,coor,ndim,npoin,ielem,rvol,rqual(ielem))
       if(rqual(ielem)>qmax)qmax=rqual(ielem) 
    enddo
    !
    !     Give diagnostics
    !
    write(*,*)'Initial quality:',qmax
    !
    !     Get the inner edges for the whole mesh 
    !  
    call gtedge(nelem,ndim,npoin,nnode,nedge,elem,coor,ledge,&
         ledglm,eltoel,nnosi,nnofa,lbhash1,lbhash2,lboup)
    !
    !     Initialize nhole
    !
    !nhole=0_ip 
    !
    !     Allocate edge quality, renumbering and active array
    !
    !allocate(lhole(nedge/100),stat=istat)
    !call memchk(zero,istat,memor_msh,'LHOLE','optim3d',lhole)
    allocate(lrenu(nedge),stat=istat)
    call memchk(zero,istat,memor_msh,'LRENU','optim3d',lrenu)
    allocate(rqualed(nedge),stat=istat)
    call memchk(zero,istat,memor_msh,'RQUALED','optim3d',rqualed)
    allocate(ledgact(nedge),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDGACT','optim3d',ledgact)
    do iedge=1,nedge
       ledgact(iedge)=.true.
    enddo
    !
    !     Then loop on iterations
    !
    do iter=1,niter 
       !
       !     Get the quality of each edge
       !
       call gtqualedg(nedge,ledge,ledglm,rqualed,nnode,nelem,&
            elem,rqual,eltoel)
       !
       !     Order the edges depending on their quality
       !
       call orderedg(npoin,rqualed,lrenu,nedge)
       !
       !     DBG
       !      
       !do iedge=1,nedge
       !   write(*,*)iedge,lrenu(iedge),rqualed(lrenu(iedge))
       !enddo
       !
       !     Build the hash table for the inner edges
       !
       call gtedghash2(lhash1,lhash2,npoin,ledge,nedge)
       !
       !     Swap the elements
       !     On output, the new edges are added in ledge but they 
       !     do not have a correct pointer to volume
       !
       call swap3dmsh(nelem,nnode,elem,ndim,npoin,coor,lrenu,ledge,ledglm,&
            eltoel,rqual,nhole,lhole,lhash1,lhash2,ledgact,nedge,&
            rqualed,qtol)
       !
       !     DBG
       !
       do ielem=1,nelem
          if(elem(1,ielem)/=0)then
             call vol3d(elem,nelem,nnode,coor,ndim,npoin,ielem,rvol)
             if(rvol<=c00)then
                write(*,*)'Error optim, negative volume 1'
                stop
             endif
             call qual3d(elem,nelem,nnode,coor,ndim,npoin,ielem,rvol,rq)
             if(abs(rq-rqual(ielem))>1.0d-02)then
                write(*,*)'Error optim different quality 1'
                write(*,*)ielem,rq,rqual(ielem)
                stop
             endif
          endif
       enddo
       !
       !     DBG:  verify mesh consistency
       !  
       call chkmsh(elem,nnode,nelem,eltoel,npoin)
       !
       !     Improve point placement
       !  
       call move3dmsh(nelem,nnode,elem,ndim,npoin,coor,lptet,eltoel,&
            rsize,rqual,lboup,nhole,qtol) 
       !
       !     DBG
       !
       do ielem=1,nelem
          if(elem(1,ielem)/=0)then
             call vol3d(elem,nelem,nnode,coor,ndim,npoin,ielem,rvol)
             if(rvol<=c00)then
                write(*,*)'Error optim, negative volume 2'
                stop
             endif
             call qual3d(elem,nelem,nnode,coor,ndim,npoin,ielem,rvol,rq)
             if(abs(rq-rqual(ielem))>1.0d-02)then
                write(*,*)'Error optim different quality 2'
                write(*,*)ielem,rq,rqual(ielem)
                stop
             endif
          endif
       enddo
       !
       !     Prepare data structure for next round
       !
       call updatnewedg(ledge,nedge,nelem,nnode,elem,npoin,ledglm,lptet,&
            ledgact,rqualed,lrenu,nhole,lbhash1,lbhash2,lboup)
       !
       !     Get current max quality
       !
       qmax=c00
       do ielem=1,nelem
          if(elem(1,ielem)/=0)then
             if(rqual(ielem)>qmax)then
                qmax=rqual(ielem)
                iemax=ielem
             endif
             !
             !     DBG
             !
             call vol3d(elem,nelem,nnode,coor,ndim,npoin,ielem,rvol)
             if(rvol<=c00)then
                write(*,*)'Error optim, negative volume'
                stop
             endif
             call qual3d(elem,nelem,nnode,coor,ndim,npoin,ielem,rvol,rq)
             if(abs(rq-rqual(ielem))>1.0d-02)then
                write(*,*)'Error optim different quality'
                write(*,*)ielem,rq,rqual(ielem)
                stop
             endif
          endif
       enddo
       write(*,*)'For this round,qmax=',qmax,iemax
       call writeoptim(nelem,npoin,nnode,ndim,elem,coor)

    enddo
    !
    !     Compact elem
    ! 
    !nelem0=nelem
    !nelem=0_ip    
    !do ielem=1,nelem0
    !   if(elem(1,ielem)/=0)then
    !      nelem=nelem+1
    !      elem(1,nelem)=elem(1,ielem) 
    !      elem(2,nelem)=elem(2,ielem) 
    !      elem(3,nelem)=elem(3,ielem) 
    !      elem(4,nelem)=elem(4,ielem) 
    !      rqual(nelem)=rqual(ielem) 
    !   endif
    !enddo
    !
    !     Give diagnostic
    !
    qmax=c00   
    do ielem=1,nelem
       if(elem(1,ielem)==0)cycle
       if(rqual(ielem)>qmax)qmax=rqual(ielem)
       !
       !     DBG compare rqual with qual
       ! 


    enddo
    write(*,*)'Final quality:',qmax

    call memchk(2_ip,istat,memor_msh,'LHASH1','optim3d',lhash1)
    deallocate(lhash1,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHASH1','optim3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHASH2','optim3d',lhash2)
    deallocate(lhash2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHASH2','optim3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'RQUALED','optim3d',rqualed)
    deallocate(rqualed,stat=istat)
    if(istat/=0) call memerr(2_ip,'RQUALED','optim3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LRENU','optim3d',lrenu)
    deallocate(lrenu,stat=istat)
    if(istat/=0) call memerr(2_ip,'LRENU','optim3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'RQUAL','optim3d',rqual)
    deallocate(rqual,stat=istat)
    if(istat/=0) call memerr(2_ip,'RQUAL','optim3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LPTET','optim3d',lptet)
    deallocate(lptet,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPTET','optim3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGACT','optim3d',ledgact)
    deallocate(ledgact,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGACT','optim3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGE','optim3d',ledge)
    deallocate(ledge,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGE','optim3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LEDGLM','optim3d',ledglm)
    deallocate(ledglm,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDGLM','optim3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LBHASH1','optim3d',lbhash1)
    deallocate(lbhash1,stat=istat)
    if(istat/=0) call memerr(2_ip,'LBHASH1','optim3d',0_ip)
    call memchk(2_ip,istat,memor_msh,'LBHASH2','optim3d',lbhash2)
    deallocate(lbhash2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LBHASH2','optim3d',0_ip)

  end subroutine optim3d

  subroutine gtedge(nelem,ndim,npoin,nnode,nedge,elem,coor,ledge,&
       ledglm,eltoel,nnosi,nnofa,lbhash1,lbhash2,lboup)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only    : memor_msh
    use mod_mshtol, only    : ptoelm,ptoelmb
    use mod_memchk
    implicit none
    integer(ip),intent(in)          :: nelem,ndim,npoin,nnode,nnosi,nnofa
    integer(ip),intent(inout)       :: nedge
    integer(ip),intent(in)          :: elem(nnode,nelem),eltoel(nnode,nelem)
    integer(ip),intent(in)          :: lboup(npoin)
    real(rp),intent(in)             :: coor(ndim,npoin)
    integer(ip),pointer             :: ledge(:,:),ledglm(:)
    integer(ip),pointer             :: lface(:,:),ptoel1(:),ptoel2(:)
    integer(ip),pointer             :: lbhash1(:),lbhash2(:),lmark(:)
    integer(ip)                     :: ielem,inode,npoin2,ipoin,isto,iface,inofa
    integer(ip)                     :: jpoin,nbedge1,sumt,iplace,jsto,kpoin,nface
    integer(ip)                     :: ip1,ip2
    logical(lg)                     :: found 
    integer(ip)                     :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))  
    integer(4)              :: istat
    !
    !     This subroutine gets the boundary and inner edges for the whole mesh
    !
    !
    !     Allocate lmark
    !  
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','gtedge',lmark)
    !
    !     First get the boundary edges
    !
    !
    !     Count the external faces
    !
    nface=0_ip 
    do ielem=1,nelem
       if(elem(1,ielem)==0)cycle
       do inode=1,nnode
          if(eltoel(inode,ielem)==0)then 
             nface=nface+1
          endif
       enddo
    enddo
    !
    !    Allocate lface
    !
    allocate(lface(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LFACE','gtedge',lface)
    !
    !     Store the external faces
    !
    nface=0_ip 
    do ielem=1,nelem
       if(elem(1,ielem)==0)cycle
       do inode=1,nnode
          if(eltoel(inode,ielem)==0)then 
             nface=nface+1
             lface(1,nface)=elem(ltab(1,inode),ielem)
             lface(2,nface)=elem(ltab(2,inode),ielem)
             lface(3,nface)=elem(ltab(3,inode),ielem)
          endif
       enddo
    enddo
    !
    !     Get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Get the hash table of the boundary edges
    !
    call gtedghash(nface,nnofa,lbhash1,lbhash2,npoin,lface,ptoel1,ptoel2)
    !
    !     Get the elements surrounding the points
    !
    call ptoelmb(elem,nelem,npoin,nnode,ptoel1,ptoel2)
    !
    !     Clean up lmark
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Count the internal edges
    ! 
    nedge=0_ip
    do ipoin=1,npoin
       lmark(ipoin)=ipoin
       do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1_ip
          ielem=ptoel1(isto)
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             if(lmark(jpoin)/=ipoin .and. jpoin>ipoin)then
                lmark(jpoin)=ipoin
                !
                !     Check if this is a boundary edge
                !
                !
                !     If ipoin or jpoin are not boundary points, then
                !     this is an inner edge
                !
                if(lboup(ipoin)/=1 .or. lboup(jpoin)/=1)then

                   nedge=nedge+1

                else
                   !
                   !     Should check boundary edge hash table
                   !
                   sumt=jpoin+ipoin
                   found=.false.
                   do jsto=lbhash2(sumt),lbhash2(sumt+1)-1
                      kpoin=lbhash1(jsto)
                      if(kpoin==ipoin)then
                         found=.true.
                         exit 
                      endif
                   enddo

                   if(found .eqv. .false.)then
                      nedge=nedge+1
                   endif
                endif
             endif
          enddo
       enddo
    enddo
    !
    !     Allocate ledge and ledglm
    !
    allocate(ledge(nnosi,nedge),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDGE','gtedge',ledge)
    allocate(ledglm(nedge),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDLEM','gtedge',ledglm)
    !
    !     Clean up lmark
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Store the internal edges and the edge to element pointer
    !
    nedge=0_ip
    do ipoin=1,npoin
       lmark(ipoin)=ipoin
       do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1_ip
          ielem=ptoel1(isto)
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             if(lmark(jpoin)/=ipoin .and. jpoin>ipoin)then
                lmark(jpoin)=ipoin
                sumt=jpoin+ipoin
                !
                !     Check if this is a boundary edge
                !
                if(lboup(ipoin)==0 .or. lboup(jpoin)==0)then

                   nedge=nedge+1
                   ledge(1,nedge)=ipoin
                   ledge(2,nedge)=jpoin
                   ledglm(nedge)=ielem

                else
                   !
                   !     Should check boundary edge hash table
                   !
                   found=.false.
                   do jsto=lbhash2(sumt),lbhash2(sumt+1)-1
                      kpoin=lbhash1(jsto)
                      if(kpoin==ipoin)then
                         found=.true.
                         exit 
                      endif
                   enddo

                   if(found .eqv. .false.)then
                      nedge=nedge+1
                      ledge(1,nedge)=ipoin
                      ledge(2,nedge)=jpoin
                      ledglm(nedge)=ielem
                   endif
                endif
             endif
          enddo
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'PTOEL1','gtedge',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','gtedge',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','gtedge',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','gtedge',0_ip)
    call memchk(2_ip,istat,memor_msh,'LFACE','gtedge',lface)
    deallocate(lface,stat=istat)
    if(istat/=0) call memerr(2_ip,'LFACE','gtedge',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','gtedge',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','gtedge',0_ip)

  end subroutine gtedge

  subroutine gtedghash(nelem,nnode,lhash1,lhash2,npoin,elem,ptoel1,ptoel2)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only    : memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)          :: nelem,npoin,nnode
    integer(ip),intent(in)          :: elem(nnode,nelem)
    integer(ip),intent(in)          :: ptoel2(npoin+1),ptoel1(*)
    integer(ip),pointer             :: lhash1(:),lhash2(:)
    integer(ip),pointer             :: lmark(:)
    integer(ip)                     :: npoin2,ipoin,isto,ielem,inode,jpoin
    integer(ip)                     :: sumt,nedge,iplace
    integer(4)              :: istat
    !
    !     This sub get the hash table (lhash2,lhash1) for the external edges 
    !     lhash1 carries the point with the lowest number   
    !
    !     Allocate lmark
    !  
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','gtedghash',lmark)
    !
    !     Allocate the hash table of the edges
    !
    npoin2=2*npoin
    allocate(lhash2(npoin2+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LHASH2','gtedghash',lhash2)
    !
    !     Count the edges
    ! 
    do ipoin=1,npoin
       lmark(ipoin)=ipoin
       do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1_ip
          ielem=ptoel1(isto)
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             if(lmark(jpoin)/=ipoin .and. jpoin>ipoin)then
                sumt=jpoin+ipoin+1_ip
                lhash2(sumt)=lhash2(sumt)+1_ip
                lmark(jpoin)=ipoin
             endif
          enddo
       enddo
    enddo
    !
    !     Add in lhash2
    !
    lhash2(1)=1_ip
    do ipoin=2,npoin2+1
       lhash2(ipoin)=lhash2(ipoin)+lhash2(ipoin-1)
    enddo
    !
    !     Allocate lhash1
    !
    nedge=lhash2(npoin2+1)-1  
    allocate(lhash1(nedge),stat=istat)
    call memchk(zero,istat,memor_msh,'LHASH1','gtedghash',lhash1)
    !
    !     Clean up lmark
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Store in lhash1
    !   
    do ipoin=1,npoin
       lmark(ipoin)=ipoin
       do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1_ip
          ielem=ptoel1(isto)
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             if(lmark(jpoin)/=ipoin .and. jpoin>ipoin)then
                sumt=jpoin+ipoin
                iplace=lhash2(sumt)
                lhash1(iplace)=ipoin
                lhash2(sumt)=iplace+1
                lmark(jpoin)=ipoin
             endif
          enddo
       enddo
    enddo
    !
    !     Reshuffle
    !
    do ipoin=npoin2+1,2,-1 
       lhash2(ipoin)=lhash2(ipoin-1)
    enddo
    lhash2(1)=1_ip

    call memchk(2_ip,istat,memor_msh,'LMARK','gtedghash',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','gtedghash',0_ip)

  end subroutine gtedghash

  subroutine gtedghash2(lhash1,lhash2,npoin,ledge,nedge)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only    : memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)          :: npoin,nedge
    integer(ip),intent(in)          :: ledge(2,nedge)
    integer(ip),pointer             :: lhash1(:,:),lhash2(:)
    integer(ip)                     :: npoin2,ipoin,isto,ielem,inode,jpoin
    integer(ip)                     :: sumt,iplace,iedge,ip1,ip2
    integer(4)              :: istat
    !
    !     This sub get the hash table (lhash2,lhash1) for the edges in ledge 
    !     lhash1(1,isto) is the point with lowest index 
    !     lhash1(2,isto) is the edge number 
    !   
    !
    !     Allocate the hash table of the edges
    !
    npoin2=2*npoin
    if(.not.associated(lhash2))then
       allocate(lhash2(npoin2+1),stat=istat)
       call memchk(zero,istat,memor_msh,'LHASH2','gtedghash2',lhash2)
    else
       call memrea(npoin2+1_ip,memor_msh,'LHASH2','gthash2',lhash2)
       lhash2=0_ip
    endif
    !
    !     Count the edges
    !
    do iedge=1,nedge
       ip1=ledge(1,iedge)
       ip2=ledge(2,iedge)
       sumt=ip1+ip2+1_ip
       lhash2(sumt)=lhash2(sumt)+1_ip
    enddo
    !
    !     Add in lhash2
    !
    lhash2(1)=1_ip
    do ipoin=2,npoin2+1
       lhash2(ipoin)=lhash2(ipoin)+lhash2(ipoin-1)
    enddo
    !
    !     Allocate lhash1
    !
    if(.not.associated(lhash1))then
       allocate(lhash1(2,nedge),stat=istat)
       call memchk(zero,istat,memor_msh,'LHASH1','gtedghash2',lhash1)
    else
       call memrea(nedge,memor_msh,'LHASH1','gthash1',lhash1)
       lhash1=0_ip
    endif
    !
    !     Store in lhash1
    !   
    do iedge=1,nedge
       ip1=ledge(1,iedge)
       ip2=ledge(2,iedge)
       sumt=ip1+ip2
       iplace=lhash2(sumt)
       lhash1(1,iplace)=ip1
       lhash1(2,iplace)=iedge
       lhash2(sumt)=lhash2(sumt)+1
    enddo
    !
    !     Reshuffle
    !
    do ipoin=npoin2+1,2,-1 
       lhash2(ipoin)=lhash2(ipoin-1)
    enddo
    lhash2(1)=1_ip

  end subroutine gtedghash2

  subroutine gtqualedg(nedge,ledge,ledglm,rqualed,nnode,nelem,elem,rqual,eltoel)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only    : memor_msh
    use mod_memchk
    implicit none
    integer(ip)             :: nedge,nelem,nnode 
    integer(ip)             :: ledge(2,nedge),ledglm(nedge),elem(nnode,nelem)
    integer(ip)             :: eltoel(nnode,nelem)
    real(rp)                :: rqualed(nedge),rqual(nelem)
    integer(ip)             :: ielem0,ielem,ienext,ipoin
    integer(ip)             :: ipos1,ipos2,iedge,ip1,ip2,idir
    real(rp)                :: qual
    integer(ip)                     :: ltab(2,4,4)=RESHAPE((/&
         0,0, 3,4, 4,2, 2,3,&
         4,3, 0,0, 1,4, 3,1,&
         2,4, 4,1, 0,0, 1,2,&
         3,2, 1,3, 2,1, 0,0 /),(/2,4,4/))
    integer(4)              :: istat
    !
    !     This sub gets the quality for each inner edge
    !     Is is assumed that ledglm is valid at this point
    ! 
    !
    do iedge=1,nedge 
       !
       !     Get the element associated with this edge
       !
       ielem0=ledglm(iedge) 
       ip1=ledge(1,iedge)
       ip2=ledge(2,iedge)

       if(ielem0<=0)then
          write(*,*)'Error in gtqualedg, ielem0=',ielem0
          stop
       endif 

       if(elem(1,ielem0)==ip1)then
          ipos1=1_ip
       else if(elem(2,ielem0)==ip1)then
          ipos1=2_ip
       else if(elem(3,ielem0)==ip1)then
          ipos1=3_ip
       else if(elem(4,ielem0)==ip1)then
          ipos1=4_ip
       else
          write(*,*)'Error in qualed 1, point not found'
          write(*,*)'iedge=',iedge,'ielem0=',ielem0
          write(*,*)'ip1=',ip1,'ip2',ip2
          write(*,*)elem(1,ielem0),elem(2,ielem0),elem(3,ielem0),elem(4,ielem0)
          stop
       endif

       if(elem(1,ielem0)==ip2)then
          ipos2=1_ip
       else if(elem(2,ielem0)==ip2)then
          ipos2=2_ip
       else if(elem(3,ielem0)==ip2)then
          ipos2=3_ip
       else if(elem(4,ielem0)==ip2)then
          ipos2=4_ip
       else
          write(*,*)'Error in qualed 2, point not found'
          write(*,*)'iedge=',iedge,'ielem0=',ielem0
          write(*,*)'ip1=',ip1,'ip2',ip2
          write(*,*)elem(1,ielem0),elem(2,ielem0),elem(3,ielem0),elem(4,ielem0)
          stop
       endif
       !
       !     Initialize quality
       !
       idir=ltab(1,ipos2,ipos1)
       ipoin=elem(ltab(2,ipos2,ipos1),ielem0)
       ielem=ielem0
       qual=rqual(ielem)
       !
       !     Loop on the elements surrounding iedge
       !
       do 
          ienext=eltoel(idir,ielem)
          !
          !     Any error
          !
          if(ienext<=0)then
             write(*,*)'Error in gtqualedg for edge:',ip1,ip2
             write(*,*)'ienext=',ienext
             stop
          endif
          !
          !     Check exit condition
          !
          if(ienext==ielem0)exit
          !
          !     Update quality
          !
          if(rqual(ienext)>qual)qual=rqual(ienext) 

          if(elem(1,ienext)==ipoin)then
             idir=1_ip
          else if(elem(2,ienext)==ipoin)then
             idir=2_ip
          else if(elem(3,ienext)==ipoin)then
             idir=3_ip
          else
             idir=4_ip
          endif

          if(eltoel(1,ienext)==ielem)then
             ipoin=1_ip
          else if(eltoel(2,ienext)==ielem)then
             ipoin=2_ip
          else if(eltoel(3,ienext)==ielem)then
             ipoin=3_ip
          else
             ipoin=4_ip
          endif

          ipoin=elem(ipoin,ienext)
          ielem=ienext

       enddo
       !
       !     Remember quality
       !
       rqualed(iedge)=qual

    enddo

  end subroutine gtqualedg

  subroutine orderedg(npoin,rqualed,lrenu,nedge)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only    : memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)          :: npoin,nedge
    integer(ip),intent(inout)       :: lrenu(nedge)
    real(rp),intent(in)             :: rqualed(nedge)
    real(rp)                        :: c00,cbig,qmax,qmin,dq,c101,c05
    real(rp)                        :: qual,dqual,qc,epsil 
    integer(ip)                     :: iedge,iplace,ipos
    integer(ip),pointer             :: ledg2(:)
    integer(4)              :: istat
    !
    !     This subroutine orders the edges depending on their quality by hashing
    !
    c00=0.0d+00
    c05=0.5d+00
    c101=1.01d+00
    cbig=1.0d+12
    epsil=0.99999d+00 
    !
    !     Allocate hash table
    !
    allocate(ledg2(nedge+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LEDG2','orderedg',ledg2)
    !
    !     Get the min/max
    !
    qmax=c00
    qmin=cbig

    do iedge=1,nedge 
       qual=rqualed(iedge) 
       if(qual>qmax)qmax=qual
       if(qual<qmin)qmin=qual
    enddo
    !
    !     Compute center
    !
    qc=(qmax+qmin)*c05
    dqual=(qmax-qc)*c101
    qmin=qc-dqual
    qmax=qc+dqual 
    !
    !     Do we really have to sort?
    !
    qmin=min(qmin,qmax*epsil)
    !
    !     Compute resolution
    !
    dqual=(qmax-qmin)/real(nedge)
    !
    !     Hash in ledg2
    ! 
    do iedge=1,nedge
       iplace=floor((qmax-rqualed(iedge))/dqual)+2_ip
       ledg2(iplace)=ledg2(iplace)+1_ip
    enddo
    !
    !     Add
    !
    ledg2(1)=1_ip 
    do iedge=2,nedge+1
       ledg2(iedge)=ledg2(iedge)+ledg2(iedge-1)
    enddo
    !
    !     Get the renumbering
    !
    do iedge=1,nedge
       iplace=floor((qmax-rqualed(iedge))/dqual)+1_ip
       ipos=ledg2(iplace)
       lrenu(ipos)=iedge
       ledg2(iplace)=ipos+1_ip
    enddo

    call memchk(2_ip,istat,memor_msh,'LEDG2','orderedg',ledg2)
    deallocate(ledg2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LEDG2','orderedg',0_ip)

  end subroutine orderedg

  subroutine swap3dmsh(nelem,nnode,elem,ndim,npoin,coor,lrenu,ledge,ledglm,&
       eltoel,rqual,nhole,lhole,lhash1,lhash2,ledgact,nedge,&
       rqualed,qtol)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only    : memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)          :: ndim,npoin,nnode
    integer(ip),intent(inout)       :: nelem,nhole,nedge
    integer(ip),intent(in)          :: lhash2(2*npoin+1)
    integer(ip),intent(in)          :: lhash1(2,*)
    integer(ip),pointer             :: lrenu(:)
    real(rp),pointer                :: rqualed(:)
    integer(ip),pointer             :: lhole(:),eltoel(:,:),elem(:,:)
    real(rp),pointer                :: rqual(:)
    integer(ip),pointer             :: ledglm(:),ledge(:,:)
    logical(lg),pointer             :: ledgact(:)
    real(rp),intent(in)             :: coor(ndim,npoin),qtol
    integer(ip)                     :: ielem,iedge,jedge,nedgnew,iopt,ncont,ncont2 
    real(rp)                        :: rq,rvol,c00 
    integer(4)              :: istat
    !
    !     This subroutine is responsible for improving the quality of a mesh 
    !     by swapping the elements around the edges. The new edges are added to the list.
    !     The old one are marked as inactive and deleted during compression
    !
    !     On input: the list of edges ledge(2,nedge)
    !
    !     On output: the list of old and new active edges  
    !
    ! 
    c00=0.0d+00
    !
    !     Initialize counter of optimized edges
    !
    ncont=0_ip
    !
    !     Set nedgnew
    !  
    nedgnew=nedge
    !
    !     Loop on edges
    !
    do iedge=1,nedge
       jedge=lrenu(iedge)
       !write(*,*)iedge,jedge,nhole+nelem,nelem,nhole
       !
       !     DBG
       ! 
       !call vol3d(elem,nelem,nnode,coor,ndim,npoin,102,rvol)
       !if(rvol<=c00)then
       !   write(*,*)'Error optim, negative volume'
       !   stop
       !endif 
       !call qual3d(elem,nelem,nnode,coor,ndim,npoin,102,rvol,rq)
       !if(abs(rq-rqual(102))>1.0d-02)then
       !   write(*,*)'Error optim different quality'
       !   write(*,*)102,rq,rqual(102)
       !   write(*,*)iedge,elem(1,102),elem(2,102),elem(3,102),elem(4,102)
       !   stop
       !endif  

       iopt=0_ip
       call swap3d(nelem,nnode,elem,ndim,npoin,coor,jedge,ledge,ledglm,eltoel,&
            rqual,nhole,lhole,lhash1,lhash2,nedgnew,ledgact,iopt,qtol) 
       ncont=ncont+iopt

       !if(iedge>2000)then 
       !call chkmsh(elem,nnode,nelem,eltoel,npoin)
       !endif
    enddo
    write(*,*)ncont,'edges optimized out of',nedge
    !
    !     Check mesh consistency
    !
    call chkmsh(elem,nnode,nelem,eltoel,npoin)
    !
    !     Compact the old (and still active) and new (and active) edges
    !  
    nedge=0_ip
    do iedge=1,nedgnew
       if(ledgact(iedge) .eqv. .true.)then
          nedge=nedge+1_ip
          ledge(1,nedge)=ledge(1,iedge)  
          ledge(2,nedge)=ledge(2,iedge)
          ledglm(nedge)=ledglm(iedge)
       endif
    enddo
    !
    !     Resize rqualedg
    !
    call memrea(nedge,memor_msh,'RQUALED','swap3dmsh',rqualed)
    call memrea(nedge,memor_msh,'LRENU','swap3dmsh',lrenu)

  end subroutine swap3dmsh

  subroutine move3dmsh(nelem,nnode,elem,ndim,npoin,coor,lptet,eltoel,rsize,rqual,&
       lboup,nhole,qtol) 
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only    : memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)          :: nelem,ndim,npoin,nnode,nhole
    integer(ip),intent(in)          :: elem(nnode,nelem)
    integer(ip),intent(inout)       :: lptet(npoin)
    integer(ip),intent(in)          :: eltoel(nnode,nelem),lboup(npoin)
    real(rp),intent(in)             :: rsize(npoin),qtol
    real(rp),intent(inout)          :: coor(ndim,npoin),rqual(nelem)
    integer(ip)                     :: ipoin,ielem,ip1,ip2,ip3,ip4,ncont,ncont2,ichk
    integer(ip),pointer             :: lmark(:) 
    real(rp)                        :: rq,rvol,c00
    integer(4)              :: istat
    !
    !     This subroutine is responsible for improving the quality of a mesh 
    !     by optimizing the point placement
    !
    !     On input: lptet >0 if the point is candiate for  further optimization
    ! 
    !     On output: lptet = 0 if not worth the optimization or not done
    !
    !
    c00=0.0d+00 
    !
    !     Allocate lmark
    !
    allocate(lmark(nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','move3dmsh',lmark)
    !
    !     Initialize the counter of optimized points ncont
    !
    ncont=0_ip 
    !
    !     Initialize the counter of tested points ncont2
    !
    ncont2=0_ip 
    !
    !     Scatter lptet
    !
    do ielem=1,nelem
       if(elem(1,ielem)/=0)then
          ip1=elem(1,ielem)
          ip2=elem(2,ielem)
          ip3=elem(3,ielem)
          ip4=elem(4,ielem)
          !
          !     Scatter if >0
          !
          if(lptet(ip1)>0)lptet(ip1)=ielem
          if(lptet(ip2)>0)lptet(ip2)=ielem
          if(lptet(ip3)>0)lptet(ip3)=ielem
          if(lptet(ip4)>0)lptet(ip4)=ielem

       endif
    enddo
    !
    !     Loop on points
    ! 
    do ipoin=1,npoin
       !
       !    Do not move boundary points
       ! 
       if(lboup(ipoin)==0)then
          !
          !     Has the point been marked for further optimization?
          !
          if(lptet(ipoin)>0)then
             !
             !     DBG
             !
             !write(*,*)ipoin
             !call vol3d(elem,nelem,nnode,coor,ndim,npoin,77,rvol)
             !if(rvol<=c00)then
             !   write(*,*)'Error optim, negative volume'
             !   stop
             !endif 
             !call qual3d(elem,nelem,nnode,coor,ndim,npoin,77,rvol,rq)
             !if(abs(rq-rqual(77))>1.0d-02)then
             !   write(*,*)'Error optim different quality'
             !   write(*,*)77,rq,rqual(77)
             !   write(*,*)ipoin,elem(1,77),elem(2,77),elem(3,77),elem(4,77)
             !   stop
             !endif 
             !
             !     Move the point to optimize quality
             !
             ichk=0_ip
             call move3d(ipoin,coor,ndim,npoin,elem,nnode,nelem,eltoel,lptet,&
                  rsize,lmark,rqual,nhole,qtol,ichk)
             ncont=ncont+ichk
             ncont2=ncont2+1_ip           

          endif
       else
          lptet(ipoin)=0_ip
       endif

    enddo

    write(*,*)ncont,'points optimized out of',ncont2

    call memchk(2_ip,istat,memor_msh,'LMARK','orderedg',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','orderedg',0_ip)

  end subroutine move3dmsh

  subroutine vol3d(elem,nelem,nnode,coor,ndim,npoin,ielem,vol)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: nelem,ndim,npoin,nnode,ielem
    integer(ip),intent(in)          :: elem(nnode,nelem)
    real(rp),intent(in)             :: coor(ndim,npoin)
    real(rp),intent(inout)          :: vol
    integer(ip)                     :: ip1,ip2,ip3,ip4
    real(rp)                        :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4    
    real(rp)                        :: x21,y21,z21,x31,y31,z31,x41,y41,z41    
    real(rp)                        :: c16
    !
    !     This subroutine computes the volume of a tetrahedron
    !
    c16=1.0d+00/6.0d+00   

    ip1=elem(1,ielem)
    ip2=elem(2,ielem)
    ip3=elem(3,ielem)
    ip4=elem(4,ielem)

    x1=coor(1,ip1)
    y1=coor(2,ip1)
    z1=coor(3,ip1)
    x2=coor(1,ip2)
    y2=coor(2,ip2)
    z2=coor(3,ip2)
    x3=coor(1,ip3)
    y3=coor(2,ip3)
    z3=coor(3,ip3)
    x4=coor(1,ip4)
    y4=coor(2,ip4)
    z4=coor(3,ip4)

    x21=x2-x1
    y21=y2-y1
    z21=z2-z1

    x31=x3-x1
    y31=y3-y1
    z31=z3-z1

    x41=x4-x1
    y41=y4-y1
    z41=z4-z1

    vol=c16*(x21*(y31*z41-z31*y41) + x31*(z21*y41-y21*z41) + x41*(y21*z31-z21*y31))

  end subroutine vol3d

  subroutine qual3d(elem,nelem,nnode,coor,ndim,npoin,ielem,vol,qual)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: nelem,ndim,npoin,nnode,ielem
    integer(ip),intent(in)          :: elem(nnode,nelem)
    real(rp),intent(in)             :: coor(ndim,npoin),vol
    real(rp),intent(inout)          :: qual
    integer(ip)                     :: ip1,ip2,ip3,ip4
    real(rp)                        :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4    
    real(rp)                        :: p1(3),p2(3),p3(3),p4(3),p5(3),p6(3),ptmp1(3)    
    real(rp)                        :: l1,l2,l3,l4,l5,l6,s1,s2,s3,s4,s,lmax    
    real(rp)                        :: c12,c30,alpha

    c12=1.0d+00/2.0d+00
    c30=3.0d+00 
    alpha=sqrt(6.0d+00)/12.0d+00
    !
    !     This subroutine computes the quality of the tetrahedron 
    !
    ip1=elem(1,ielem)
    ip2=elem(2,ielem)
    ip3=elem(3,ielem)
    ip4=elem(4,ielem)
    !
    !     Get the sides
    !
    p1(1)=coor(1,ip2)-coor(1,ip1)
    p1(2)=coor(2,ip2)-coor(2,ip1)
    p1(3)=coor(3,ip2)-coor(3,ip1)
    p2(1)=coor(1,ip3)-coor(1,ip1)
    p2(2)=coor(2,ip3)-coor(2,ip1)
    p2(3)=coor(3,ip3)-coor(3,ip1)
    p3(1)=coor(1,ip4)-coor(1,ip1)
    p3(2)=coor(2,ip4)-coor(2,ip1)
    p3(3)=coor(3,ip4)-coor(3,ip1)
    p4(1)=coor(1,ip3)-coor(1,ip2)
    p4(2)=coor(2,ip3)-coor(2,ip2)
    p4(3)=coor(3,ip3)-coor(3,ip2)
    p5(1)=coor(1,ip4)-coor(1,ip2)
    p5(2)=coor(2,ip4)-coor(2,ip2)
    p5(3)=coor(3,ip4)-coor(3,ip2)
    p6(1)=coor(1,ip4)-coor(1,ip3)
    p6(2)=coor(2,ip4)-coor(2,ip3)
    p6(3)=coor(3,ip4)-coor(3,ip3)
    !
    !     Get the length of the sides
    !
    l1=sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
    l2=sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))
    l3=sqrt(p3(1)*p3(1)+p3(2)*p3(2)+p3(3)*p3(3))
    l4=sqrt(p4(1)*p4(1)+p4(2)*p4(2)+p4(3)*p4(3))
    l5=sqrt(p5(1)*p5(1)+p5(2)*p5(2)+p5(3)*p5(3))
    l6=sqrt(p6(1)*p6(1)+p6(2)*p6(2)+p6(3)*p6(3))
    !
    !     Get the max
    !
    if(l1>l2)then
       lmax=l1
    else
       lmax=l2
    endif

    if(l3>lmax)lmax=l3
    if(l4>lmax)lmax=l4
    if(l5>lmax)lmax=l5
    if(l6>lmax)lmax=l6
    !
    !     Get the area of the faces
    !
    ptmp1(1)= p6(2)*p4(3)-p6(3)*p4(2)
    ptmp1(2)=-p6(1)*p4(3)+p6(3)*p4(1)
    ptmp1(3)= p6(1)*p4(2)-p6(2)*p4(1)
    s1=sqrt(ptmp1(1)*ptmp1(1)+ptmp1(2)*ptmp1(2)+ptmp1(3)*ptmp1(3))

    ptmp1(1)= p3(2)*p2(3)-p3(3)*p2(2)
    ptmp1(2)=-p3(1)*p2(3)+p3(3)*p2(1)
    ptmp1(3)= p3(1)*p2(2)-p3(2)*p2(1)
    s2=sqrt(ptmp1(1)*ptmp1(1)+ptmp1(2)*ptmp1(2)+ptmp1(3)*ptmp1(3))

    ptmp1(1)= p1(2)*p3(3)-p1(3)*p3(2)
    ptmp1(2)=-p1(1)*p3(3)+p1(3)*p3(1)
    ptmp1(3)= p1(1)*p3(2)-p1(2)*p3(1)
    s3=sqrt(ptmp1(1)*ptmp1(1)+ptmp1(2)*ptmp1(2)+ptmp1(3)*ptmp1(3))

    ptmp1(1)= p2(2)*p1(3)-p2(3)*p1(2)
    ptmp1(2)=-p2(1)*p1(3)+p2(3)*p1(1)
    ptmp1(3)= p2(1)*p1(2)-p2(2)*p1(1)
    s4=sqrt(ptmp1(1)*ptmp1(1)+ptmp1(2)*ptmp1(2)+ptmp1(3)*ptmp1(3))

    s=(s1+s2+s3+s4)*c12

    qual=alpha*lmax*s/(c30*vol)

  end subroutine qual3d

  subroutine move3d(ipoin,coor,ndim,npoin,elem,nnode,nelem,eltoel,lptet,rsize,lmark,&
       rqual,nhole,qtol,ichk)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: nelem,ndim,npoin,nnode,ipoin,nhole
    integer(ip),intent(in)          :: elem(nnode,nelem)
    integer(ip),intent(in)          :: eltoel(nnode,nelem)
    integer(ip),intent(inout)       :: lmark(nelem),lptet(npoin),ichk
    real(rp),intent(inout)          :: coor(ndim,npoin),rqual(nelem)
    real(rp),intent(in)             :: rsize(npoin),qtol
    integer(ip)                     :: ip1,ip2,ip3,ienext
    integer(ip)                     :: ie,iter,maxiter,nball
    real(rp)                        :: pold(ndim),rnopoold(ndim),c10 
    real(rp)                        :: d1,d2,d3,rpoin(ndim),c00
    real(rp)                        :: nrmal1(ndim),nrmal2(ndim),optdir(3,6),lenloc,epsil,rl 
    real(rp)                        :: qini,qiniloc,pinit(ndim),coef1,coef2,TOLGAIN,c05
    real(rp)                        :: coefmin,lenabs,qini0,ploc(ndim),q1,qend
    real(rp)                        :: rpoint(ndim),rinit(ndim),rloc(ndim),qual,pnew(3)
    integer(ip)                     :: chkimp,idir,idiropt,ldir(6),check,ihostt 
    integer(ip)                     :: isloc,iimpr,tab(6),ielem,iball,ipos,jchk 
    integer(ip)                     :: ineigh1,ineigh,optim 
    integer(ip)                     :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))  
    integer(ip),parameter           :: mball=100
    integer(ip)                     :: lball(mball),elemt(4,mball)
    real(rp)                        :: rq(mball),rq2(mball),modul(mball)
    !
    !     This subroutine optimizes the point position by improving the quality of the
    !     patch around the point
    !
    !     On input: -ipoin, the point to be optimized
    !               -lptet contains an element containing ipoin
    !
    !     On output: lptet=0 for not optimized points
    !                lptet>0 for optimized points   
    !
    c10=1.0d+00
    c00=0.0d+00
    c05=0.5d+00
    epsil=1.0d-12
    coef1=0.1d+00
    coef2=2.0d+00
    coefmin=coef1*coef1!*coef1
    TOLGAIN=0.1d+00
    tab(1)=2_ip
    tab(2)=1_ip
    tab(3)=4_ip
    tab(4)=3_ip
    tab(5)=6_ip
    tab(6)=5_ip
    !
    !     First loop, get the ball of the point. Check if it is worth the trouble 
    !
    ie=lptet(ipoin)
    !
    !     DBG
    !
    jchk=0_ip
    if(elem(1,ie)==ipoin)then
       jchk=1_ip 
    else if(elem(2,ie)==ipoin)then
       jchk=1_ip 
    else if(elem(3,ie)==ipoin)then
       jchk=1_ip 
    else if(elem(4,ie)==ipoin)then
       jchk=1_ip
    endif 
    if(jchk==0)then
       write(*,*)'Error in move3d, ipoin=',ipoin
       stop
    endif 
      
    elemt(1,1)=elem(1,ie)
    elemt(2,1)=elem(2,ie)
    elemt(3,1)=elem(3,ie)
    elemt(4,1)=elem(4,ie)
    nball=1_ip
    iball=0_ip
    lball(1)=ie
    lmark(ie)=ipoin

    do

       if(iball==nball)exit
       iball=iball+1_ip
       ie=lball(iball)
       !
       !     Get the point in the element
       !
       if(elem(1,ie)==ipoin)then
          ipos=1_ip
       else if(elem(2,ie)==ipoin)then
          ipos=2_ip
       else if(elem(3,ie)==ipoin)then
          ipos=3_ip
       else
          ipos=4_ip
       endif
       !
       !     Get the surrounding elements
       !
       ineigh1=eltoel(ltab(1,ipos),ie)
       if(ineigh1==0)cycle
       if(lmark(ineigh1)/=ipoin)then
          lmark(ineigh1)=ipoin
          nball=nball+1   
          if(nball>mball)then
             write(*,*)'Error in move3d, nball==mball'
             stop
          endif        
          lball(nball)=ineigh1
          elemt(1,nball)=elem(1,ineigh1)
          elemt(2,nball)=elem(2,ineigh1)
          elemt(3,nball)=elem(3,ineigh1)
          elemt(4,nball)=elem(4,ineigh1)
       endif
       ineigh1=eltoel(ltab(2,ipos),ie)
       if(ineigh1==0)cycle
       if(lmark(ineigh1)/=ipoin)then
          lmark(ineigh1)=ipoin
          nball=nball+1          
          if(nball>mball)then
             write(*,*)'Error in move3d, nball==mball'
             stop
          endif        
          lball(nball)=ineigh1
          elemt(1,nball)=elem(1,ineigh1)
          elemt(2,nball)=elem(2,ineigh1)
          elemt(3,nball)=elem(3,ineigh1)
          elemt(4,nball)=elem(4,ineigh1)
       endif
       ineigh1=eltoel(ltab(3,ipos),ie)
       if(ineigh1==0)cycle
       if(lmark(ineigh1)/=ipoin)then
          lmark(ineigh1)=ipoin
          nball=nball+1          
          if(nball>mball)then
             write(*,*)'Error in move3d, nball==mball'
             stop
          endif        
          lball(nball)=ineigh1
          elemt(1,nball)=elem(1,ineigh1)
          elemt(2,nball)=elem(2,ineigh1)
          elemt(3,nball)=elem(3,ineigh1)
          elemt(4,nball)=elem(4,ineigh1)
       endif

    enddo
    !
    !     Compute initial quality
    !
    qini=c00

    do ie=1,nball
       ielem=lball(ie)
       if(rqual(ielem)>qini)qini=rqual(ielem)
    enddo
    !
    !     Is it worth the trouble?
    !
    if(qini<qtol)then
       lptet(ipoin)=0_ip   
       return
    endif
    !
    !     Remember the position
    !
    pold(1)=coor(1,ipoin)
    pold(2)=coor(2,ipoin)
    pold(3)=coor(3,ipoin)
    !
    !     Fill optimal direction
    !
    optdir(1,1)= 1
    optdir(2,1)= 0
    optdir(3,1)= 0
    optdir(1,2)=-1
    optdir(2,2)= 0
    optdir(3,2)= 0
    optdir(1,3)= 0
    optdir(2,3)= 1
    optdir(3,3)= 0
    optdir(1,4)= 0
    optdir(2,4)=-1
    optdir(3,4)= 0
    optdir(1,5)= 0
    optdir(2,5)= 0
    optdir(3,5)= 1
    optdir(1,6)= 0
    optdir(2,6)= 0
    optdir(3,6)=-1

    lenloc=rsize(ipoin)
    !
    !     Initialize optimization
    !
    lenloc=lenloc*coef1
    lenabs=coef1
    qini0=qini
    pinit(1)=pold(1)
    pinit(2)=pold(2)
    pinit(3)=pold(3)
    !
    !     Initialize flag
    ! 
    check=0_ip
    !
    !     Loop iteratively
    !   
    maxiter=10_ip

    do iter=1,maxiter
       !
       !     Reset optimization flag
       !
       chkimp=0_ip
       qiniloc=qini
       !
       !     Loop on direction
       !
       do idir=1,6
          !
          !     Compute new position
          !
          pnew(1)=pinit(1)+lenloc*optdir(1,idir)
          pnew(2)=pinit(2)+lenloc*optdir(2,idir)
          pnew(3)=pinit(3)+lenloc*optdir(3,idir)
          !
          !     Transfer data to ipoin
          !
          coor(1,ipoin)=pnew(1)
          coor(2,ipoin)=pnew(2)
          coor(3,ipoin)=pnew(3)
          !
          !     Verify correctness and quality 
          !
          q1=c00

          do ie=1,nball
             call vol3d(elemt,nball,nnode,coor,ndim,npoin,ie,modul(ie))
             if(modul(ie)<=c00)exit
             call qual3d(elemt,nball,nnode,coor,ndim,npoin,ie,modul(ie),rq(ie))
             if(rq(ie)>qini)exit
             if(rq(ie)>q1)q1=rq(ie)
          enddo

          if(ie>nball)then
             !
             !     Successful, remember the point
             !
             ploc(1)=pnew(1) 
             ploc(2)=pnew(2) 
             ploc(3)=pnew(3)
             qini=q1 
             chkimp=1_ip
             idiropt=idir
             do ie=1,nball
                rq2(ie)=rq(ie)
             enddo

          endif

       enddo
       !
       !     Did we optimize something
       !
       if(chkimp==0)then
          !
          !     Optimization unsuccessfull for this round 
          !

          !
          !     Initialize ldir
          !
          ldir(1)=0_ip
          ldir(2)=0_ip
          ldir(3)=0_ip
          ldir(4)=0_ip
          ldir(5)=0_ip
          ldir(6)=0_ip
          !
          !      Reduce the steps
          !
          lenloc=lenloc*coef1
          lenabs=lenabs*coef1
          !
          !     Did we reduce too much ?
          !
          if(lenabs<coefmin)exit

       else
          !
          !     Optimization successfull
          !
          !
          !     Is it worth the trouble?
          !
          !if((qiniloc-qini)<TOLGAIN)exit
          !
          !     Mark the oposite direction
          !
          ldir(tab(idiropt))=0_ip
          !
          !     Assign the valid new point
          !
          pinit(1)=ploc(1)
          pinit(2)=ploc(2)
          pinit(3)=ploc(3)
          check=1_ip

       endif
    enddo

    coor(1,ipoin)=pinit(1)
    coor(2,ipoin)=pinit(2)
    coor(3,ipoin)=pinit(3)
    !
    !     Did we optimize something?
    !
    if(check==0)then
       lptet(ipoin)=0_ip
       return
    endif
    !write(*,*)'Optimized:',qini0,qini
    !
    !     Update quality
    !
    ichk=1_ip   
    do ie=1,nball
       ielem=lball(ie)
       rqual(ielem)=rq2(ie)
    enddo

  end subroutine move3d

  subroutine swap3d(nelem,nnode,elem,ndim,npoin,coor,iedge,ledge,ledglm,eltoel,&
       rqual,nhole,lhole,lhash1,lhash2,nedgnew,ledgact,iopt,qtol)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only    : memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)          :: ndim,npoin,nnode,iedge
    integer(ip),intent(inout)       :: nelem,nhole,nedgnew,iopt
    integer(ip),intent(in)          :: lhash2(2*npoin+1)
    integer(ip),intent(in)          :: lhash1(2,*)
    integer(ip),pointer             :: ledglm(:),ledge(:,:)
    logical(lg),pointer             :: ledgact(:)
    integer(ip),pointer             :: lhole(:),eltoel(:,:),elem(:,:)
    real(rp),pointer                :: rqual(:)
    real(rp),intent(in)             :: coor(ndim,npoin),qtol
    real(rp)                        :: qini,c00,qend,modul
    real(rp)                        :: qnew,q1,c20
    real(rp)                        :: newvol(112),newqual(112),goodqual(112)
    integer(ip)                     :: nshell,ishell,lshell(100),ncat
    integer(ip)                     :: ielem,icat,badtry(56),ip1,ip2,ip3,itri,ntri,ntri2
    integer(ip)                     :: newelem(nnode,100),knew,k,iposa,iposb,ipoin
    integer(ip)                     :: goodelem(nnode,100),ielem0,ipa,ipb,idir,ipnt
    integer(ip)                     :: defpnt(100),externa(100),externb(100),ienext
    integer(ip)                     :: viewa(2,100),viewb(2,100),iview,inedg2
    integer(ip)                     :: iview1,iview2,ipos1,ipos2,inedg,ielem1a,ielem1b
    integer(ip)                     :: ielem2a,ielem2b,iedg,nnew,neleh,itet1,itet2
    integer(ip)                     :: cont2,nbtry,ipos,narray,ielem1,ielem2
    integer(4)                      :: istat
    integer(ip)                     :: ltab(2,4,4)=RESHAPE((/&
         0,0, 3,4, 4,2, 2,3,&
         4,3, 0,0, 1,4, 3,1,&
         2,4, 4,1, 0,0, 1,2,&
         3,2, 1,3, 2,1, 0,0 /),(/2,4,4/))

    badtry=0_ip
    nbtry=56_ip
    c00=0.0d+00
    c20=2.0d+00
    newvol=c00 
    newqual=c00
    narray=112_ip
    nnew=100_ip
    neleh=nelem
    !
    !     First get the shell
    !
    ielem0=ledglm(iedge) 
    ipa=ledge(1,iedge)
    ipb=ledge(2,iedge)

    if(elem(1,ielem0)==ipa)then
       iposa=1_ip
    else if(elem(2,ielem0)==ipa)then
       iposa=2_ip
    else if(elem(3,ielem0)==ipa)then
       iposa=3_ip
    else if(elem(4,ielem0)==ipa)then
       iposa=4_ip
    else
       write(*,*)'Error in swap3d 1, point not found'
       stop
    endif

    if(elem(1,ielem0)==ipb)then
       iposb=1_ip
    else if(elem(2,ielem0)==ipb)then
       iposb=2_ip
    else if(elem(3,ielem0)==ipb)then
       iposb=3_ip
    else if(elem(4,ielem0)==ipb)then
       iposb=4_ip
    else
       write(*,*)'Error in swap3d 2, point not found'
       stop
    endif
    !
    !     Initialize 
    !
    idir=ltab(1,iposb,iposa)
    ipnt=ltab(2,iposb,iposa)
    ipoin=elem(ipnt,ielem0)

    ielem=ielem0
    nshell=1_ip
    lshell(1)=ielem0
    defpnt(1)=elem(idir,ielem0)
    externa(1)=eltoel(iposa,ielem0)
    externb(1)=eltoel(iposb,ielem0)
    !
    !     Loop on the elements surrounding iedge
    !
    do 
       ienext=eltoel(idir,ielem)
       !
       !     Check exit condition
       !
       if(ienext==ielem0)exit

       if(elem(1,ienext)==ipoin)then
          idir=1_ip
       else if(elem(2,ienext)==ipoin)then
          idir=2_ip
       else if(elem(3,ienext)==ipoin)then
          idir=3_ip
       else if(elem(4,ienext)==ipoin)then
          idir=4_ip
       else
          write(*,*)'Error in swap3d 3, point not found'
          stop
       endif

       if(eltoel(1,ienext)==ielem)then
          ipoin=1_ip
       else if(eltoel(2,ienext)==ielem)then
          ipoin=2_ip
       else if(eltoel(3,ienext)==ielem)then
          ipoin=3_ip
       else if(eltoel(4,ienext)==ielem)then
          ipoin=4_ip
       else
          write(*,*)'Error in swap3d 4, point not found'  
          stop
       endif
       !
       !     Remember needed information
       !
       nshell=nshell+1_ip
       lshell(nshell)=ienext
       defpnt(nshell)=elem(idir,ienext)
       externa(nshell)=eltoel(ltab(1,ipoin,idir),ienext)
       externb(nshell)=eltoel(ltab(2,ipoin,idir),ienext)
       !
       !     Update for next round
       !
       ipoin=elem(ipoin,ienext)
       ielem=ienext

    enddo
    !
    !     For the moment until 6 elements
    !
    if(nshell>6)then
       ledgact(iedge)=.false.
       return
    endif
    !
    !     Get the quality of the initial shell
    ! 
    qini=c00
    do ishell=1,nshell
       ielem=lshell(ishell)
       if(rqual(ielem)>qini)qini=rqual(ielem)  
    enddo
    !
    !     Is it worth the trouble?
    !
    if(qini<qtol)then
       ledgact(iedge)=.false.
       return
    endif
    !
    !     Transfer qini to qnew, the best actual configuration 
    !
    qnew=qini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !     Construct the new tetrahedra to test
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ncat=lcatn(nshell)
    ntri=ltri(nshell) 
    ntri2=2*ntri

    do icat=1,ncat
       !
       !     Did we already test this possibility?
       !
       if(badtry(icat)==1)cycle
       !
       !     Update the quality for this round
       !
       qend=c00 
       cont2=0_ip  
       !
       !     Test the new tetrahedra
       !
       do itri=1,ntri
          !
          !     Position of this triangle in the complete list
          !   
          ipos=(icat-1)*ntri+itri
          !
          !     Position of the first tetra in the complete list
          !   
          itet1=2*(ipos-1)+1_ip
          !
          !     Upper tetrahedron
          !
          ip1=defpnt(indx(1,itri,icat,nshell))
          ip2=defpnt(indx(2,itri,icat,nshell))
          ip3=defpnt(indx(3,itri,icat,nshell))

          cont2=cont2+1_ip
          newelem(1,cont2)=ip1
          newelem(2,cont2)=ip2
          newelem(3,cont2)=ip3
          newelem(4,cont2)=ipb
          !
          !     Did we compute the volume already?
          !
          if(newvol(itet1)==c00)then
             !
             !     Compute volume
             ! 
             call vol3d(newelem,nnew,nnode,coor,ndim,npoin,cont2,newvol(itet1))
             !
             !     Mark the next combinations
             !
             call mrkarray(ipos,narray,newvol,1_ip,newvol(itet1),nshell)

          endif
          !
          !     Do we have a valid element?
          !
          if(newvol(itet1)<=c00)then
             !
             !     Mark all the combinations containing this element
             !
             call mrkbad(ipos,nbtry,nshell,badtry)
             !
             !     Goto the end of the round
             !
             qend=c20*qnew
             exit 

          endif
          !
          !     Did we compute the quality already?
          !
          if(newqual(itet1)==c00)then
             !
             !     Compute quality
             ! 
             call qual3d(newelem,nnew,nnode,coor,ndim,npoin,cont2,newvol(itet1),newqual(itet1))
             !
             !     Mark the next combinations
             !
             call mrkarray(ipos,narray,newqual,1_ip,newqual(itet1),nshell)

          endif
          !
          !     Is it worth following?
          !
          if(newqual(itet1)>qnew)then
             !
             !     Mark all the combinations containing this element
             !
             call mrkbad(ipos,nbtry,nshell,badtry)
             !
             !     Goto the end of the round
             !
             qend=c20*qnew
             exit

          endif
          !
          !     Update local quality 
          !
          if(newqual(itet1)>qend)qend=newqual(itet1)
          !
          !     Lower tetrahedron
          !
          cont2=cont2+1_ip
          newelem(1,cont2)=ipa
          newelem(2,cont2)=ip1
          newelem(3,cont2)=ip2
          newelem(4,cont2)=ip3
          itet2=itet1+1_ip
          !
          !     Did we compute the volume already?
          !
          if(newvol(itet2)==c00)then
             !
             !     Compute volume
             ! 
             call vol3d(newelem,nnew,nnode,coor,ndim,npoin,cont2,newvol(itet2))
             !
             !     Mark the next combinations
             !
             call mrkarray(ipos,narray,newvol,2_ip,newvol(itet2),nshell)

          endif
          !
          !     Do we have a valid element?
          !
          if(newvol(itet2)<=c00)then
             !
             !     Mark all the combinations containing this element
             !
             call mrkbad(ipos,nbtry,nshell,badtry)
             !
             !     Goto the end of the round
             !
             qend=c20*qnew
             exit

          endif
          !
          !     Did we compute the quality already?
          !
          if(newqual(itet2)==c00)then
             !
             !     Compute quality
             ! 
             call qual3d(newelem,nnew,nnode,coor,ndim,npoin,cont2,newvol(itet2),newqual(itet2))
             !
             !     Mark the next combinations
             !
             call mrkarray(ipos,narray,newqual,2_ip,newqual(itet2),nshell)

          endif
          !
          !     Is it worth following?
          !
          if(newqual(itet2)>qnew)then
             !
             !     Mark all the combinations containing this element
             !
             call mrkbad(ipos,nbtry,nshell,badtry)
             !
             !     Goto the end of the round
             !
             qend=c20*qnew
             exit

          endif
          !
          !     Update local quality 
          !
          if(newqual(itet2)>qend)qend=newqual(itet2)

       enddo
       !
       !     Compare the quality of this configuration
       !
       if(qend>qnew)cycle 
       !
       !     We have a better configuration, dump it
       !
       qnew=qend
       knew=icat

       cont2=0_ip
       do itri=1,ntri
          cont2=cont2+1_ip
          ipos=(icat-1)*ntri+itri
          itet1=2*(ipos-1)+1_ip
          goodqual(cont2)=newqual(itet1)
          goodelem(1,cont2)=newelem(1,cont2)
          goodelem(2,cont2)=newelem(2,cont2)
          goodelem(3,cont2)=newelem(3,cont2)
          goodelem(4,cont2)=newelem(4,cont2) 

          cont2=cont2+1_ip
          itet2=itet1+1_ip
          goodqual(cont2)=newqual(itet2)
          goodelem(1,cont2)=newelem(1,cont2)
          goodelem(2,cont2)=newelem(2,cont2)
          goodelem(3,cont2)=newelem(3,cont2)
          goodelem(4,cont2)=newelem(4,cont2) 

       enddo

    enddo
    !
    !     Did we optimize something?
    !
    if(qnew==qini)then
       ledgact(iedge)=.false.
       return
    endif
    iopt=1
    !
    !     As the optimization was successfull, write the new mesh
    !
    !write(*,*)'Edge optimized nshell=',nshell,qini,qend

    !if(nshell==5)then
    !   call writeshell(nelem,npoin,nnode,ndim,elem,coor,lshell,nshell)
    !endif
    !
    !     Special case 3-->2
    !
    if(nshell==3)then
       nhole=nhole+1_ip
       call memrea(nelem,memor_msh,'LHOLE','swap3d',lhole)
       lhole(nhole)=lshell(3)
       elem(1,lshell(3))=0_ip
    endif
    !
    !     Get the new places
    !
    do ishell=nshell+1,ntri2
       if(nhole>0)then
          lshell(ishell)=lhole(nhole)
          nhole=nhole-1_ip
       else
          nelem=nelem+1_ip
          lshell(ishell)=nelem
       endif
    enddo
    !
    !     Resize elem
    !
    call memrea(nelem,memor_msh,'ELEM','swap3d',elem)
    call memrea(nelem,memor_msh,'ELEM','swap3d',rqual)
    call memrea(nelem,memor_msh,'ELTOEL','swap3d',eltoel)
    !
    !     Copy the elements back in elem
    !
    do itri=1,ntri2
       ielem=lshell(itri)
       elem(1,ielem) = goodelem(1,itri)
       elem(2,ielem) = goodelem(2,itri)
       elem(3,ielem) = goodelem(3,itri)
       elem(4,ielem) = goodelem(4,itri)
       rqual(ielem)  = goodqual(itri)
    enddo
    !
    !     Update external neighborings
    !
    do ishell=1,nshell
       !
       !     Get the triangle facing the ishell side
       !
       itri=2*(exte(1,ishell,knew,nshell)-1)+1
       !
       !    From this deduce the new tetrahedra
       !
       ielem1=lshell(itri)
       ielem2=lshell(itri+1)
       !
       !     Get the position of the viewer in the triangle
       !
       ipos=exte(2,ishell,knew,nshell)
       !
       !     From this update eltoel for the elements in the ring
       !
       eltoel(ipos,ielem1)=externa(ishell)
       eltoel(ipos+1,ielem2)=externb(ishell)
       !
       !     Get the viewer for the outer elements
       !  
       call viewer(neleh,eltoel,nnode,iview,lshell(ishell),externa(ishell))
       viewa(1,ishell)=ielem1
       viewa(2,ishell)=iview 
       call viewer(neleh,eltoel,nnode,iview,lshell(ishell),externb(ishell))
       viewb(1,ishell)=ielem2
       viewb(2,ishell)=iview
       !
       !     Update edge to element pointer
       !
       call updatedg(nedgnew,ledge,ledglm,ipb,defpnt(ishell),ielem1,lhash1,lhash2,npoin,ledgact)
       call updatedg(nedgnew,ledge,ledglm,ipa,defpnt(ishell),ielem2,lhash1,lhash2,npoin,ledgact)
       if(ishell==nshell)then
          ipoin=defpnt(1)
       else   
          ipoin=defpnt(ishell+1)
       endif
       call updatedg(nedgnew,ledge,ledglm,defpnt(ishell),ipoin,ielem1,lhash1,lhash2,npoin,ledgact)

    enddo
    !
    !     Second loop to avoid side effects
    !
    do ishell=1,nshell
       ielem=externa(ishell)
       if(ielem/=0)then
          eltoel(viewa(2,ishell),ielem)=viewa(1,ishell)
       endif
       ielem=externb(ishell)
       if(ielem/=0)then
          eltoel(viewb(2,ishell),ielem)=viewb(1,ishell)
       endif
    enddo
    !
    !     Tetrahedra facing each other
    !
    do itri=1,ntri
       ipos=2*(itri-1)+1
       ielem1=lshell(ipos)
       ielem2=lshell(ipos+1)
       eltoel(4,ielem1)=ielem2
       eltoel(1,ielem2)=ielem1
    enddo
    !
    !     Update internal neighborings
    !
    inedg2=ntri-1
    do inedg=1,inedg2
       !
       !     Get the first element with inner relations
       !
       ipos1=2*(inte1(1,inedg,knew,nshell)-1)+1_ip
       ielem1a=lshell(ipos1)
       ielem1b=lshell(ipos1+1)
       iview1=inte1(2,inedg,knew,nshell)
       !
       !     Get the second element with inner relations
       ! 
       ipos2=2*(inte2(1,inedg,knew,nshell)-1)+1_ip
       ielem2a=lshell(ipos2)
       ielem2b=lshell(ipos2+1)
       iview2=inte2(2,inedg,knew,nshell)
       !
       !     Update eltoel
       ! 
       eltoel(iview1,ielem1a)=ielem2a
       eltoel(iview1+1,ielem1b)=ielem2b
       eltoel(iview2,ielem2a)=ielem1a
       eltoel(iview2+1,ielem2b)=ielem1b

    enddo
    !
    !     Add the new edges 
    !
    do iedg=1,inedg2
       ip1=defpnt(indedg(1,iedg,knew,nshell))
       ip2=defpnt(indedg(2,iedg,knew,nshell))
       ielem=lshell(2*(indedglm(iedg,knew,nshell)-1)+1)
       call addedg(nedgnew,ledge,ledglm,ip1,ip2,ielem,ledgact,lhash1,lhash2,npoin)
    enddo
    !
    !     Mark the edge as deleted
    !
    ledgact(iedge)=.false.     

  end subroutine swap3d

  subroutine mrkbad(ipos,nbtry,nshell,badtry)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: nbtry,nshell
    integer(ip),intent(in)          :: ipos
    integer(ip),intent(inout)       :: badtry(nbtry)
    integer(ip)                     :: icomb,ipost 

    if(ipos==0)then
       write(*,*)'Error in mrkbad'
       stop
    endif

    ipost=pos(ipos,nshell)

    do while(ipost/=0)

       icomb=(ipost-1)/nshell+1
       badtry(icomb)=1_ip
       ipost=pos(ipost,nshell)

    enddo

  end subroutine mrkbad

  subroutine mrkarray(ipos,narray,array,ishift,rval,nshell)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: narray,nshell,ishift
    integer(ip),intent(in)          :: ipos
    real(rp),intent(inout)          :: array(narray)
    real(rp),intent(in)             :: rval
    integer(ip)                     :: iplace,ipost 

    if(ipos==0)then
       write(*,*)'Error in mrkarray'
       stop
    endif

    ipost=pos(ipos,nshell)

    do while(ipost/=0)

       iplace=2*(ipost-1)+ishift
       array(iplace)=rval
       ipost=pos(ipost,nshell)

    enddo

  end subroutine mrkarray

  subroutine viewer(nelem,eltoel,nnode,ipos,ielem1,ielem2)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nnode,nelem,ielem1,ielem2
    integer(ip),intent(in)      :: eltoel(nnode,nelem)
    integer(ip),intent(out)     :: ipos

    if(ielem2==0)then
       ipos=0 
       return
    endif

    if(eltoel(1,ielem2)==ielem1)then
       ipos=1
    else if(eltoel(2,ielem2)==ielem1)then
       ipos=2
    else if(eltoel(3,ielem2)==ielem1)then
       ipos=3
    else
       ipos=4
    endif

  end subroutine viewer

  subroutine addedg(nedgnew,ledge,ledglm,pa,pb,ielem,ledgact,lhash1,lhash2,npoin)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only    : memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: pa,pb,ielem,npoin
    integer(ip),intent(inout)   :: nedgnew
    integer(ip),pointer         :: ledge(:,:),ledglm(:)
    integer(ip),intent(in)      :: lhash2(2*npoin+1)
    integer(ip),intent(in)      :: lhash1(2,*)
    logical(lg),pointer         :: ledgact(:)
    integer(4)                  :: istat
    integer(ip)                 :: minp,maxp,sump,iedge,jedge
    !
    !     This sub checks if an edge already exists and if not,
    !     add it to the database
    !
    sump=pa+pb
    minp=min(pa,pb) 
    !
    !     Look for the edge
    !
    do iedge=lhash2(sump),lhash2(sump+1)-1
       if(lhash1(1,iedge)==minp)then
          jedge=lhash1(2,iedge)
          ledgact(jedge)=.true.
          return
       endif
    enddo
    !
    !     As the edge has not been found, add it in the edge list
    !     not in the hash list
    !
    nedgnew=nedgnew+1 
    call memrea(nedgnew,memor_msh,'LEDGE','addedg',ledge)
    call memrea(nedgnew,memor_msh,'LEDGLM','addedg',ledglm)
    call memrea(nedgnew,memor_msh,'LEDGACT','addedg',ledgact)

    if(pa<pb)then
       minp=pa
       maxp=pb
    else
       minp=pb
       maxp=pa
    endif

    ledge(1,nedgnew)=minp 
    ledge(2,nedgnew)=maxp 
    ledglm(nedgnew)=ielem 
    ledgact(nedgnew)=.true. 

  end subroutine addedg

  subroutine updatedg(nedgnew,ledge,ledglm,pa,pb,ielem,lhash1,lhash2,npoin,ledgact)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only    : memor_msh
    use mod_memchk
    implicit none
    integer(ip),intent(in)      :: pa,pb,ielem,npoin,nedgnew
    integer(ip),intent(in)      :: ledge(2,nedgnew),lhash2(2*npoin+1)
    integer(ip),intent(in)      :: lhash1(2,*)
    integer(ip),intent(inout)   :: ledglm(nedgnew)
    logical(lg),intent(inout)   :: ledgact(nedgnew)
    integer(ip)                 :: sump,minp,maxp,iedge,jedge

    sump=pa+pb
    minp=min(pa,pb) 

    do iedge=lhash2(sump),lhash2(sump+1)-1
       if(lhash1(1,iedge)==minp)then
          jedge=lhash1(2,iedge)
          ledglm(jedge)=ielem
          ledgact(jedge)=.true.
          exit
       endif
    enddo

  end subroutine updatedg

  subroutine readoptim(nelem,npoin,nnode,ndim,elem,coor,filename)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(inout)  :: nelem,npoin
    integer(ip), intent(in)     :: nnode,ndim
    character*80,intent(in)     :: filename
    integer(ip),pointer         :: elem(:,:)
    real(rp),pointer            :: coor(:,:)
    integer(ip)                 :: ielem,jelem,ipoin,i
    integer(ip)                 :: iface,jface,isurf,imat,icond,ichk
    integer(4)                  :: istat
    character*80                :: ntext
    integer(ip), pointer        :: ptoel1(:),ptoel2(:) 
    !
    !     Open file
    !    
    open(unit=77,file=trim(filename)//'.dat',status='old')
    !
    !     Read points
    !
    read(77,*)ntext,npoin
    allocate(coor(ndim,npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'COOR','readgeo',coor)

    do i=1,npoin
       read(77,*) ipoin,coor(1,i),coor(2,i),coor(3,i)
       !read(77,*) coor(1,i),coor(2,i),coor(3,i)
    end do
    !
    !     Read elements
    !
    read(77,*)ntext,nelem
    allocate(elem(nnode,nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'ELEM','readgeo',elem) 
    do ielem=1,nelem
       read(77,*) jelem,elem(1,ielem),elem(2,ielem),elem(3,ielem),&
            elem(4,ielem),imat
    enddo
    !
    !     Close input file
    !
    close(77)

  end subroutine readoptim

  subroutine writeoptim(nelem,npoin,nnode,ndim,elem,coor)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)     :: nelem,npoin
    integer(ip), intent(in)     :: nnode,ndim
    integer(ip),intent(in)      :: elem(nnode,nelem)
    real(rp),intent(in)         :: coor(ndim,npoin)
    integer(ip)                 :: ipoin,ielem,ncont
    real(rp)                    :: rx,ry,rz

    open(unit=50,file='outoptim.msh',status='unknown')
    rewind 50

1   format('MESH dimension 3 ElemType Tetrahedra Nnode 4')
2   format('Coordinates')
3   format('#node number   coord_x   coord_y  coord_z')
100 format(i10,3e20.10)
200 format(5i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(50,1)
    write(50,2)
    write(50,3)
    do  ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(50,100)ipoin,rx,ry,rz
    enddo

    write(50,4)
    write(50,5)
    ncont=0_ip
    do  ielem=1,nelem
       if(elem(1,ielem)==0)cycle
       ncont=ncont+1
       write(50,200)ncont,elem(1,ielem),elem(2,ielem),elem(3,ielem),elem(4,ielem)
    enddo
    write(50,6)
    close(50)

  end subroutine writeoptim

  subroutine chkmsh(elem,nnode,nelem,eltoel,npoin)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)     :: nelem,nnode,npoin
    integer(ip),intent(in)      :: elem(nnode,nelem),eltoel(nnode,nelem)
    integer(ip),pointer         :: lmark(:) 
    integer(ip)                 :: ipa,ipb,ipc,ip1,ip2,ip3,ielem
    integer(ip)                 :: ineigh,inode,ipos,ierr
    integer(4)                      :: istat
    integer(ip)                     :: ltab(3,4)=RESHAPE((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))  

    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','readgeo',lmark) 

    do ielem=1,nelem
       if(elem(1,ielem)/=0)then
          do inode=1,nnode
             ineigh=eltoel(inode,ielem)
             if(ineigh/=0)then
                ipa=elem(ltab(1,inode),ielem) 
                ipb=elem(ltab(2,inode),ielem) 
                ipc=elem(ltab(3,inode),ielem) 

                lmark(ipa)=1
                lmark(ipb)=1
                lmark(ipc)=1

                ierr=0

                if(eltoel(1,ineigh)==ielem)then
                   ipos=1 
                else if(eltoel(2,ineigh)==ielem)then
                   ipos=2 
                else if(eltoel(3,ineigh)==ielem)then
                   ipos=3 
                else if(eltoel(4,ineigh)==ielem)then
                   ipos=4 
                else
                   write(*,*)'Error in chkmsh 1'
                   ierr=1
                endif

                ip1=elem(ltab(1,ipos),ineigh) 
                ip2=elem(ltab(2,ipos),ineigh) 
                ip3=elem(ltab(3,ipos),ineigh) 

                if(lmark(ip1)/=1)then
                   write(*,*)'Error in chkmsh 2'
                   ierr=1
                endif
                if(lmark(ip2)/=1)then
                   write(*,*)'Error in chkmsh 3'
                   ierr=1
                endif
                if(lmark(ip3)/=1)then
                   write(*,*)'Error in chkmsh 4'
                   ierr=1
                endif

                if(ierr==1)then
                   write(*,*)'ielem=',ielem,'ineigh=',ineigh
                   write(*,*)elem(1,ielem),elem(2,ielem),elem(3,ielem),elem(4,ielem)  
                   write(*,*)eltoel(1,ielem),eltoel(2,ielem),eltoel(3,ielem),eltoel(4,ielem)  
                   write(*,*)elem(1,ineigh),elem(2,ineigh),elem(3,ineigh),elem(4,ineigh)  
                   write(*,*)eltoel(1,ineigh),eltoel(2,ineigh),eltoel(3,ineigh),eltoel(4,ineigh) 
                   stop 
                endif

                lmark(ipa)=0
                lmark(ipb)=0
                lmark(ipc)=0

             endif
          enddo
       endif
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','chkmsh',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','chkmsh',0_ip)

  end subroutine chkmsh

  subroutine writeshell(nelem,npoin,nnode,ndim,elem,coor,lshell,nshell)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)     :: nelem,npoin
    integer(ip), intent(in)     :: nnode,ndim,nshell
    integer(ip),intent(in)      :: elem(nnode,nelem),lshell(nshell)
    real(rp),intent(in)         :: coor(ndim,npoin)
    integer(ip)                 :: ipoin,ielem,ishell
    real(rp)                    :: rx,ry,rz

    open(unit=50,file='outshell.msh',status='unknown')
    rewind 50

1   format('MESH dimension 3 ElemType Tetrahedra Nnode 4')
2   format('Coordinates')
3   format('#node number   coord_x   coord_y  coord_z')
100 format(i10,3e20.10)
200 format(5i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(50,1)
    write(50,2)
    write(50,3)
    do  ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(50,100)ipoin,rx,ry,rz
    enddo

    write(50,4)
    write(50,5)
    do  ishell=1,nshell
       ielem=lshell(ishell)   
       write(50,200)ielem,elem(1,ielem),elem(2,ielem),elem(3,ielem),elem(4,ielem)
    enddo
    write(50,6)
    close(50)

  end subroutine writeshell

  subroutine updatnewedg(ledge,nedge,nelem,nnode,elem,npoin,ledglm,lptet,&
       ledgact,rqualed,lrenu,nhole,lbhash1,lbhash2,lboup)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    use mod_mshtol, only    : ptoelmb
    implicit none
    !
    !     This sub is responsible for updating the edge to volume pointer
    !     for the new edges not updated during swap3d and the edges touching
    !     optimized points 
    !
    !     On input:  -ledge(1:nedge) contains the new edges without valid pointer to volume
    !                -lptet >0 for candidate points to optimization through movement
    ! 
    !     on output: the new edges in ledge with valid pointer to volume 
    ! 
    !
    integer(ip),intent(in)   :: nelem,nnode,npoin,nhole
    integer(ip),intent(inout):: nedge
    integer(ip),intent(in)   :: elem(nnode,nhole+nelem)
    integer(ip),intent(in)   :: lbhash2(2*npoin+1),lbhash1(*),lboup(npoin)
    integer(ip),intent(inout):: lptet(npoin)
    integer(ip),pointer      :: ledge(:,:),ledglm(:),lrenu(:)
    logical(lg),pointer      :: ledgact(:)
    real(rp),pointer         :: rqualed(:)
    integer(ip),pointer      :: ptoel1(:),ptoel2(:),lmark(:)
    integer(ip),pointer      :: lhash1(:,:),lhash2(:)
    integer(ip)              :: ip1,ip2,ipmin,nelem1,nelem2,isto      
    integer(ip)              :: inode,jpoin,ielem,iedge,ipoin      
    integer(ip)              :: isum,jsto
    logical(lg)              :: ifound     
    integer(4)               :: istat
    !
    !     Allocate lmark
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','updatnewedg',lmark) 
    !
    !     Get the elements surrounding the points 
    !     (special version as some elements may have been deleted)
    !
    call ptoelmb(elem,nelem,npoin,nnode,ptoel1,ptoel2)
    !
    !     Loop on edges
    !
    do iedge=1,nedge
       ip1=ledge(1,iedge) 
       ip2=ledge(2,iedge) 
       !
       !     Mark lmark
       !  
       lmark(ip1)=iedge
       !
       !     Loop on elements surrounding point
       ! 
       do isto=ptoel2(ip1),ptoel2(ip1+1)-1
          ielem=ptoel1(isto)
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             if(jpoin>ip1)then
                if(lmark(jpoin)/=iedge)then
                   lmark(jpoin)=iedge
                   if(jpoin==ip2)then 
                      !
                      !     Edge found, get the pointer
                      !
                      ledglm(iedge)=ielem 
                      goto 100

                   endif

                endif
             endif
          enddo
       enddo
       !
       !     Error, edge not found
       ! 
       write(*,*)'Error in updatnewedg, edge not found' 
       stop 

100    continue    

    enddo
    !
    !     Then get the hash table of the current edges
    !  
    call gtedghash2(lhash1,lhash2,npoin,ledge,nedge)
    !
    !     Clean up lmark
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Finally add in ledge all the edges of the points moved by optimization
    !     not already present
    !
    do ipoin=1,npoin
       if(lptet(ipoin)>0)then
          lmark(ipoin)=ipoin
          !
          !     Loop on elements surrouding points
          ! 
          do isto=ptoel2(ipoin),ptoel2(ipoin+1)-1
             ielem=ptoel1(isto)
             do inode=1,nnode 
                jpoin=elem(inode,ielem)
                if(jpoin>ipoin)then
                   if(lmark(jpoin)/=ipoin)then         
                      lmark(jpoin)=ipoin
                      !
                      !     Check if we already have this edge
                      ! 
                      isum=ipoin+jpoin
                      ifound=.false.
                      do jsto=lhash2(isum),lhash2(isum+1)-1
                         if(lhash1(1,jsto)==ipoin)then
                            ifound=.true.
                            exit
                         endif
                      enddo
                      !
                      !     Did we find this edge?
                      !
                      if(ifound .eqv. .false.)then
                         !
                         !     Is this edge a boundary edge?
                         !
                         if(lboup(jpoin)==0 .or. lboup(ipoin)==0)then 

                            nedge=nedge+1   
                            call memrea(nedge,memor_msh,'LEDGE','updatnewedg',ledge)
                            call memrea(nedge,memor_msh,'LEDGLM','updatnewedg',ledglm)
                            call memrea(nedge,memor_msh,'LEDGACT','updatnewedg',ledgact)
                            call memrea(nedge,memor_msh,'RQUALED','updatnewedg',rqualed)
                            call memrea(nedge,memor_msh,'LRENU','updatnewedg',lrenu)
                            ledge(1,nedge)=ipoin
                            ledge(2,nedge)=jpoin
                            ledglm(nedge)=ielem
                            ledgact(nedge)=.true. 

                         else
                            !
                            !    Should check boundary edge hash table
                            !
                            do jsto=lbhash2(isum),lbhash2(isum+1)-1
                               if(lbhash1(jsto)==ipoin)then
                                  ifound=.true.
                                  exit
                               endif
                            enddo
                            !
                            !     Did we find this edge?
                            !
                            if(ifound .eqv. .false.)then

                               nedge=nedge+1   
                               call memrea(nedge,memor_msh,'LEDGE','updatnewedg',ledge)
                               call memrea(nedge,memor_msh,'LEDGLM','updatnewedg',ledglm)
                               call memrea(nedge,memor_msh,'LEDGACT','updatnewedg',ledgact)
                               call memrea(nedge,memor_msh,'RQUALED','updatnewedg',rqualed)
                               call memrea(nedge,memor_msh,'LRENU','updatnewedg',lrenu)
                               ledge(1,nedge)=ipoin
                               ledge(2,nedge)=jpoin
                               ledglm(nedge)=ielem
                               ledgact(nedge)=.true.

                            endif
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
       endif
    enddo
    !
    !     And mark lptet as >0 for all this new edges
    !
    do iedge=1,nedge
       ip1=ledge(1,iedge) 
       ip2=ledge(2,iedge) 
       lptet(ip1)=1_ip
       lptet(ip2)=1_ip
    enddo

    call memchk(2_ip,istat,memor_msh,'LHASH1','updatedgnew',lhash1)
    deallocate(lhash1,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHASH1','updatedgnew',0_ip)
    call memchk(2_ip,istat,memor_msh,'LHASH2','updatedgnew',lhash2)
    deallocate(lhash2,stat=istat)
    if(istat/=0) call memerr(2_ip,'LHASH2','updatedgnew',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','updatnewedg',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','updatnewedg',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL1','updatnewedg',ptoel1)
    deallocate(ptoel1,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL1','updatnewedg',0_ip)
    call memchk(2_ip,istat,memor_msh,'PTOEL2','updatnewedg',ptoel2)
    deallocate(ptoel2,stat=istat)
    if(istat/=0) call memerr(2_ip,'PTOEL2','updatnewedg',0_ip)

  end subroutine updatnewedg

end module mod_optim
