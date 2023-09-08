subroutine qua_reabcs()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_reabcs
  ! NAME
  !    qua_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions for the Shrodinger
  !    equation.
  !
  !    * For conditions on boundaries, bvnat(iboun)
  !
  !         Neumann ... fixbo(iboun)=2 -> bvnat(1,iboun)
  !         Robin ..... fixbo(iboun)=3 -> bvnat(3,iboun)
  !
  !    * For conditions on nodes, Neumann and Robin conditions are stored
  !      temporarily in bvess_qua(ipoin) and then changed to conditions on 
  !      boundaries  in the routine qua_bcntoe (called at the end) 
  !
  !    * Conditions on nodes have priority over conditions on boundaries.
  !      This is done and explained in qua_bcntoe. CHANGE?
  ! OUTPUT 
  ! USES
  !    qua_membcs
  !    ecoute
  !    reacod
  !    runend
  !    reacod
  ! USED BY
  !    qua_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_quanty
  use mod_ecoute, only : ecoute
  implicit none
  integer(ip)  :: ipoin,ncodf,nbcod

  if( INOTSLAVE .and. .not. READ_AND_RUN() ) then
     !
     ! Initializations
     !
     kfl_conbc_qua = 1             ! Constant boundary conditions
     !
     ! Read flags
     !
     ncodf         =  1            ! Phi has 1 degree of freedom 
     nbcod         = -1

     call ecoute('qua_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('qua_reabcs')
     end do

     if(exists('NONCO')) then
        kfl_conbc_qua = 0
        call qua_membcs(2_ip)
     else
        kfl_conbc_qua = 1
     end if
     !
     ! Allocate memory
     !
     call qua_membcs(1_ip)
     do ipoin = 1,npoin
        kfl_fixno_qua(1,ipoin) = -1
     end do
     !
     ! Read data
     !
     call ecoute('qua_reabcs')
     do while( words(1) /= 'ENDBO' )

        if( words(1) == 'CODES' .and. exists('NODES') ) then 
           !
           ! User-defined codes on nodes ACA!!!!!!!!!! codigo dirichlet
           !
           if( kfl_conbc_qua == 0 ) then
              iffun     =  1
              kfl_funno => kfl_funno_qua
           else
              iffun      =  0
           end if
           ifloc     =  0
           ifbop     =  0
           kfl_fixno => kfl_fixno_qua ! numero
           bvess     => bvess_qua     ! valor de dirich
           call reacod(1_ip)

        else if(words(1)=='CODES'.and.exists('BOUND')) then
           !
           ! User-defined codes on boundaries
           !          
           ! kfl_fixbo => kfl_fixbo_qua
           ! bvnat     => bvnat_qua(:,:,1)
           ! call reacod(2_ip)

        else if(words(1)=='CODES') then
           call runend('QUA_REABCS: SPECIFY IF CODES ARE APPLIED TO NODES OR BOUNDARIES')

        else if(words(1)=='INITI') then 
           !
           ! Initial condition
           !
           if(words(2)=='CONST') then
              kfl_inico_qua = 1
           else if(words(2)=='WARMB') then
              kfl_inico_qua = 2
           end if
           ! do ipara=1,10
           !   bvcoe_qua(ipara)=param(ipara+2)
           ! end do

        else if(words(1)=='ONNOD') then
           !
           ! On nodes
           !
           call ecoute('qua_reabcs')
           do while(words(1)/='ENDON')
              ipoin                = int(param(1))
              kfl_fixno_qua(1,ipoin) = int(param(2))
              if(kfl_conbc_qua==0) then
                 kfl_funno_qua(ipoin) = int(param(4))
                 if(kfl_fixno_qua(1,ipoin)==1.and.kfl_funno_qua(ipoin)/=0) then
                    bvess_qua(ipoin,1)=param(3)
                    bvess_qua(ipoin,2)=param(3)
                 else 
                    bvess_qua(ipoin,1)=param(3)
                 end if
              else
                 bvess_qua(ipoin,1)=param(3)
              end if
              call ecoute('qua_reabcs')
           end do

        endif
        call ecoute('qua_reabcs')
     enddo
  end if

end subroutine qua_reabcs
