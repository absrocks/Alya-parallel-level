!------------------------------------------------------------------------
!> @addtogroup Insitu 
!> @ingroup    Postprocess 
!> @{
!> @file    ins_turnon.f90
!> @date    20/04/2017
!> @author  Vishal Mehta
!> @brief   Turnon insitu
!> @details turnon insitu
!> @}
!------------------------------------------------------------------------
subroutine ins_turnon()
  use def_insitu
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE,PAR_SUM,PAR_GATHERV, PAR_ALLGATHER
  use mod_parall,         only : PAR_COMM_CURRENT, PAR_UNIVERSE_SIZE, PAR_MY_UNIVERSE_RANK, PAR_CODE_SIZE
  use def_domain,         only : ndime,meshe,npoin
  use def_kermod,         only : ndivi
  use def_master,         only : INOTMASTER,kfl_paral,npoi1,npoi2,npoi3,IMASTER
  implicit none
  integer(ip)                  :: mynpoin
  integer(ip),pointer          :: dummr(:)
  
  call PAR_COMM_RANK_AND_SIZE(PAR_COMM_CURRENT,PAR_CURRENT_RANK_INS,PAR_CURRENT_SIZE_INS)

  nullify(npoin_all)
  nullify(dummr)
  nullify(advec)
  nullify(elmag)
  nullify(global_num)
  nullify(elmag_buf)
  
  allocate(advec(ndime,npoin))
  allocate(npoin_all(PAR_CURRENT_SIZE_INS))
  npoin_all = 0
  
  NDglobal = npoi1 + ( npoi3 - npoi2 ) + 1
  mynpoin = NDglobal

  call PAR_SUM(NDglobal)
  call PAR_ALLGATHER(mynpoin,npoin_all,1,'IN MY CODE',PAR_COMM_CURRENT)

  if (INOTMASTER) then
     allocate(global_num(mynpoin))
     allocate(elmag(npoin))
     allocate(elmag_buf(mynpoin))
  else
     allocate(global_num(NDglobal))
     allocate(elmag(NDglobal))
  end if
  elmag = 0
  global_num = 0
  advec = 0
  
  
  if ( INOTMASTER ) then
     global_num(1:npoi1) = meshe(ndivi) % lninv_loc(1:npoi1)
     global_num(npoi2:npoi3) = meshe(ndivi) % lninv_loc(npoi2:npoi3)
  end if

  if(INOTMASTER) then
     npoin_all = 0
     call PAR_GATHERV(global_num,dummr,npoin_all,'IN MY CODE')

  else
     call PAR_GATHERV(dummr,global_num,npoin_all,'IN MY CODE')
     
  end if

#ifdef INVIZ
  call launchviz(PAR_MY_UNIVERSE_RANK,PAR_UNIVERSE_SIZE,PAR_CURRENT_RANK_INS,PAR_CURRENT_SIZE_INS)

  if (INOTMASTER) then
     call updatedata(global_num,elmag_buf,mynpoin,ndime)
  else
     call updatedata(global_num,elmag,NDglobal,ndime)
  end if
  
  call updateframe()
#endif
  
end subroutine ins_turnon
