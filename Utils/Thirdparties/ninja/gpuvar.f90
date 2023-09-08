MODULE GPUVARCG
  use def_kintyp,         only : ip,rp
  use def_master,         only : INOTMASTER,IMASTER,INOTSLAVE,ISEQUEN
  use def_master,         only : NPOIN_TYPE
  use def_master,         only : npoi1,npoi2,npoi3
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  IMPLICIT NONE
  integer*8 :: d_val,d_x,d_r,d_p,d_u,d_s,d_off,d_dia
  integer*8 :: d_col,d_row,descr,handlesp,handlebl
  INTEGER :: ti = 1
  real*8 :: tolsta = 1.0e-8_rp
  real*8,allocatable :: dia(:)
END MODULE GPUVARCG

MODULE GPUVARPIPECG
  use def_kintyp,         only : ip,rp
  use def_master,         only : INOTMASTER,IMASTER,INOTSLAVE,ISEQUEN
  use def_master,         only : NPOIN_TYPE
  use def_master,         only : npoi1,npoi2,npoi3
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  IMPLICIT NONE
  integer*8 :: d_val,d_x,d_r,d_w,d_u,d_z,d_q,d_s,d_p,d_dia,d_m,d_n,d_off,str1,str2
  integer*8 :: d_col,d_row,descr,handlesp,handlebl1,handlebl2
  INTEGER :: tiiii = 1
  real*8 :: tolsta = 1.0e-8_rp
  real*8,allocatable :: dia(:)
END MODULE GPUVARPIPECG
 
MODULE GPUVARGMRES
  use def_kintyp,         only : ip,rp
  use def_master,         only : INOTMASTER,IMASTER,INOTSLAVE,ISEQUEN
  use def_master,         only : NPOIN_TYPE
  use def_master,         only : npoi1,npoi2,npoi3
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  IMPLICIT NONE
  integer*8 :: d_val,d_x,d_b,d_Krylov,d_y,d_M,d_r,d_Ax,d_w,d_off
  real*8,allocatable    :: c_Krylov(:,:),c_r(:),c_Ax(:),c_w(:)
  integer*8 :: d_col,d_row,descr,handlesp,handlebl,strcomp,strcopy
  integer*8 :: d_colcsr,d_rowcsr,d_valcsr  
  INTEGER :: tiii = 1
  real*8 :: tolsta = 1.0e-8_rp
  real*8, dimension(:), allocatable :: dia, cs,sn, s, y,forawhile
  real*8, dimension(:,:), allocatable :: H
END MODULE GPUVARGMRES

MODULE GPUVARDEFLCG
  use def_kintyp,         only : ip,rp
  use def_master,         only : INOTMASTER,IMASTER,INOTSLAVE,ISEQUEN
  use def_master,         only : NPOIN_TYPE
  use def_master,         only : npoi1,npoi2,npoi3
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  IMPLICIT NONE
  integer*8 :: d_val,d_x,d_r,d_p,d_u,d_s,d_off,d_dia,d_grpmat,d_lgrp,d_iagro,d_jagro,d_grplhs,d_grprhs
  integer*8 :: d_col,d_row,descr,handlesp,handlebl,strcomp,strcopy
  INTEGER :: tii = 1
  real*8 :: tolsta = 1.0e-8_rp
  real*8,allocatable :: dia(:),grpmat(:),vecallred(:),mu(:),murhs(:)
END MODULE GPUVARDEFLCG

MODULE GPUVARAMGX
  use def_kintyp,         only : ip,rp
  use def_master,         only : INOTMASTER,IMASTER,INOTSLAVE,ISEQUEN
  use def_master,         only : NPOIN_TYPE
  use def_master,         only : npoi1,npoi2,npoi3
  use def_domain,         only : npoin
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  IMPLICIT NONE
  INTEGER :: initflag = 1
  integer*8 :: d_val,d_x,d_r,d_dia
  integer*8 :: d_col,d_row
  integer*8 :: config,resource,mathand,rhshand,unkhand,solhand
  integer*4, allocatable :: ja_global(:)
  
END MODULE GPUVARAMGX

module vargpuprecon
  integer*8             :: d_dia,d_rasmat,d_rasia,d_rasja,d_map,d_permvec
  integer*8             :: handsol,matdescbasezero,matdescbaseone,solinfo,datain,workin,qrbuff,handsp
  integer*8             :: d_csrrasmat,d_csrrasia,d_csrrasja
  integer*4,allocatable :: permvec(:),permbuff(:),map(:),rasia(:),rasja(:)
  INTEGER*8             :: permbuffsz
  integer*4             :: rasprepflag=0,idprecon,nrows,NNZB,nvar,permuteflag=0
  real * 8              :: tolerance = 1.0e-8
end module vargpuprecon

module gpuvarqr
  use def_kintyp,         only : ip,rp
  use def_master,         only : INOTMASTER,IMASTER,INOTSLAVE,ISEQUEN
  use def_master,         only : NPOIN_TYPE
  use def_master,         only : npoi1,npoi2,npoi3
  integer*8             :: d_dia,d_rasmat,d_rasia,d_rasja,d_map,d_permvec
  integer*8             :: handsol,matdesc,solinfo,datain,workin,qrbuff,handsp
  integer*8             :: d_csrrasmat,d_csrrasia,d_csrrasja
  integer*4,allocatable :: permvec(:),permbuff(:),map(:)
  integer*4,allocatable :: rasia(:),rasja(:)
  INTEGER*8             :: permbuffsz
  integer*4             :: rasprepflag=0,idprecon,nrows,NNZB,nvar
  real * 8              :: tolerance = 1.0e-8
end module gpuvarqr

MODULE GPUVARDEFLPIPECG
  use def_kintyp,         only : ip,rp
  use def_master,         only : INOTMASTER,IMASTER,INOTSLAVE,ISEQUEN
  use def_master,         only : NPOIN_TYPE
  use def_master,         only : npoi1,npoi2,npoi3
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  IMPLICIT NONE
  integer*8 :: d_val,d_x,d_r,d_w,d_u,d_z,d_q,d_s,d_p,d_dia,d_m,d_n,d_off,d_grpmat,d_lgrp,d_iagro,d_jagro,d_grplhs,d_grprhs
  integer*8 :: strcopy,strcomp
  integer*8 :: d_col,d_row,descr,handlesp,handlebl1,handlebl2
  INTEGER :: ti = 1
  real*8 :: tolsta = 1.0e-8_rp
  real*8,allocatable :: dia(:),grpmat(:),vecallred(:),mu(:),murhs(:)
END MODULE GPUVARDEFLPIPECG

module nastinvar
  use def_kintyp,         only : ip,rp
  use def_master,         only : INOTMASTER,IMASTER,INOTSLAVE,ISEQUEN
  use def_domain,         only : NDIME,NPOIN,MNODE,NELEM,MGAUS
  use def_nastin,         only : ncomp_nsi
  IMPLICIT NONE
  integer*8     :: dcoord,dltype,dlnods,dgravi_nsi,delmda_gpvol,delmda_gpcar
  integer*8     :: ddt_rho_nsi,dmu_rho_nsi,drhsid,dveloc
end module nastinvar
