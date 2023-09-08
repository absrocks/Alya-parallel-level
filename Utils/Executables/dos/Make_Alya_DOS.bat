@echo off

if "%1" == ""      goto version32
if "%1" == "clean" goto clean
if "%1" == "64"    goto version64
if "%1" == "32"    goto version32

:version32
set LINK_OPT=/c /nopad-source /include:"Sources\kernel\defmod" /include:"Sources\kernel\memory" /include:"Sources\modules\nstinc" /include:"Sources\modules\temper" /include:"Sources\modules\turbul" /include:"Sources\modules\codire"  /include:"Sources\modules\alefor" 
goto compila

:version64
set LINK_OPT=/c /O3 /Qparallel /Qpar_threshold100 /Qpar_report1 /4I8 /4R8 /4L132 /nbs /nopad_source /include:"Sources\kernel\defmod" /include:"Sources\kernel\memory" /include:"Sources\modules\nstinc"  /include:"Sources\modules\temper" /include:"Sources\modules\turbul"  /include:"Sources\modules\codire"   /include:"Sources\modules\alefor" 
goto compila

:clean

echo ----------------------------------------------------------    
echo Deleting objects
echo ----------------------------------------------------------    

del *.obj 
goto end

:compila

echo ----------------------------------------------------------    
echo Compiling modules
echo ----------------------------------------------------------    

ifort  %LINK_OPT% Sources/kernel/defmod/def_domain.f90
ifort  %LINK_OPT% Sources/kernel/defmod/def_inpout.f90
ifort  %LINK_OPT% Sources/kernel/defmod/def_kintyp.f90
ifort  %LINK_OPT% Sources/kernel/defmod/def_master.f90
ifort  %LINK_OPT% Sources/kernel/defmod/def_parame.f90
ifort  %LINK_OPT% Sources/kernel/defmod/def_solver.f90 
ifort  %LINK_OPT% Sources/kernel/memory/Mod_memchk.f90 
ifort  %LINK_OPT% Sources/modules/alefor/def_alefor.f90
ifort  %LINK_OPT% Sources/modules/nstinc/def_nstinc.f90
ifort  %LINK_OPT% Sources/modules/temper/def_temper.f90
ifort  %LINK_OPT% Sources/modules/turbul/def_turbul.f90
ifort  %LINK_OPT% Sources/kernel/master/Mod_iofile.f90
ifort  %LINK_OPT% Sources/kernel/mathru/Mod_gradie.f90
ifort  %LINK_OPT% Sources/kernel/outrut/Mod_output.f90
ifort  %LINK_OPT% Sources/kernel/outrut/Mod_postpr.f90  

echo ----------------------------------------------------------    
echo Compiling kernel
echo ----------------------------------------------------------    

ifort  %LINK_OPT% Sources/kernel/domain/arrind.f90
ifort  %LINK_OPT% Sources/kernel/domain/bouder.f90
ifort  %LINK_OPT% Sources/kernel/domain/bounor.f90
ifort  %LINK_OPT% Sources/kernel/domain/cartbo.f90
ifort  %LINK_OPT% Sources/kernel/domain/cderda.f90
ifort  %LINK_OPT% Sources/kernel/domain/chaord.f90
ifort  %LINK_OPT% Sources/kernel/domain/chenor.f90
ifort  %LINK_OPT% Sources/kernel/domain/conbou.f90
ifort  %LINK_OPT% Sources/kernel/domain/connbo.f90
ifort  %LINK_OPT% Sources/kernel/domain/connpo.f90
ifort  %LINK_OPT% Sources/kernel/domain/countv.f90
ifort  %LINK_OPT% Sources/kernel/domain/crcsol.f90
ifort  %LINK_OPT% Sources/kernel/domain/cresla.f90
ifort  %LINK_OPT% Sources/kernel/domain/cshder.f90
ifort  %LINK_OPT% Sources/kernel/domain/defmsh.f90
ifort  %LINK_OPT% Sources/kernel/domain/Domain.f90
ifort  %LINK_OPT% Sources/kernel/domain/dombou.f90
ifort  %LINK_OPT% Sources/kernel/domain/outdom.f90
ifort  %LINK_OPT% Sources/kernel/domain/elmdel.f90
ifort  %LINK_OPT% Sources/kernel/domain/elmder.f90
ifort  %LINK_OPT% Sources/kernel/domain/elmhes.f90
ifort  %LINK_OPT% Sources/kernel/domain/elmlen.f90
ifort  %LINK_OPT% Sources/kernel/domain/elsest.f  
ifort  %LINK_OPT% Sources/kernel/domain/extcog.f90
ifort  %LINK_OPT% Sources/kernel/domain/extnor.f90
ifort  %LINK_OPT% Sources/kernel/domain/finbou.f90
ifort  %LINK_OPT% Sources/kernel/domain/jacobs.f90
ifort  %LINK_OPT% Sources/kernel/domain/mescek.f90
ifort  %LINK_OPT% Sources/kernel/domain/nortri.f90
ifort  %LINK_OPT% Sources/kernel/domain/nzecof.f90
ifort  %LINK_OPT% Sources/kernel/domain/opfdom.f90
ifort  %LINK_OPT% Sources/kernel/domain/permat.f90
ifort  %LINK_OPT% Sources/kernel/domain/poscog.f90
ifort  %LINK_OPT% Sources/kernel/domain/reagat.f90
ifort  %LINK_OPT% Sources/kernel/domain/reageo.f90
ifort  %LINK_OPT% Sources/kernel/domain/reaset.f90
ifort  %LINK_OPT% Sources/kernel/domain/reaske.f90
ifort  %LINK_OPT% Sources/kernel/domain/reasla.f90
ifort  %LINK_OPT% Sources/kernel/domain/reastr.f90
ifort  %LINK_OPT% Sources/kernel/domain/rulepw.f90
ifort  %LINK_OPT% Sources/kernel/domain/rupclo.f90
ifort  %LINK_OPT% Sources/kernel/domain/rupope.f90
ifort  %LINK_OPT% Sources/kernel/domain/ruqclo.f90
ifort  %LINK_OPT% Sources/kernel/domain/ruqope.f90
ifort  %LINK_OPT% Sources/kernel/domain/rutclo.f90
ifort  %LINK_OPT% Sources/kernel/domain/rutope.f90
ifort  %LINK_OPT% Sources/kernel/domain/setdim.f90
ifort  %LINK_OPT% Sources/kernel/domain/setext.f90
ifort  %LINK_OPT% Sources/kernel/domain/setlbe.f90
ifort  %LINK_OPT% Sources/kernel/domain/shafga.f90
ifort  %LINK_OPT% Sources/kernel/domain/shafun.f90
ifort  %LINK_OPT% Sources/kernel/domain/shaga1.f90
ifort  %LINK_OPT% Sources/kernel/domain/shaga2.f90
ifort  %LINK_OPT% Sources/kernel/domain/shaga3.f90
ifort  %LINK_OPT% Sources/kernel/domain/shape1.f90
ifort  %LINK_OPT% Sources/kernel/domain/shape2.f90
ifort  %LINK_OPT% Sources/kernel/domain/shape3.f90
ifort  %LINK_OPT% Sources/kernel/domain/sslcon.f90
ifort  %LINK_OPT% Sources/kernel/domain/sslcr2.f90
ifort  %LINK_OPT% Sources/kernel/domain/sslcro.f90
ifort  %LINK_OPT% Sources/kernel/domain/sslgct.f90
ifort  %LINK_OPT% Sources/kernel/domain/sslmsh.f90
ifort  %LINK_OPT% Sources/kernel/domain/sslord.f90
ifort  %LINK_OPT% Sources/kernel/domain/sslpry.f90
ifort  %LINK_OPT% Sources/kernel/domain/updmsh.f90

ifort  %LINK_OPT% Sources/kernel/master/blocko.f90
ifort  %LINK_OPT% Sources/kernel/master/Conblk.f90
ifort  %LINK_OPT% Sources/kernel/master/Concou.f90
ifort  %LINK_OPT% Sources/kernel/master/cputab.f90                                                
ifort  %LINK_OPT% Sources/kernel/master/Doiter.f90
ifort  %LINK_OPT% Sources/kernel/master/Endste.f90
ifort  %LINK_OPT% Sources/kernel/master/Begste.f90
ifort  %LINK_OPT% Sources/kernel/master/inirun.f90
ifort  %LINK_OPT% Sources/kernel/master/mediso.f90
ifort  %LINK_OPT% Sources/kernel/master/memunk.f90
ifort  %LINK_OPT% Sources/kernel/master/Alya.f90
ifort  %LINK_OPT% Sources/kernel/master/Newmsh.f90
ifort  %LINK_OPT% Sources/kernel/master/openfi.f90
ifort  %LINK_OPT% Sources/kernel/master/outerr.f90
ifort  %LINK_OPT% Sources/kernel/master/outmem.f90
ifort  %LINK_OPT% Sources/kernel/master/readat.f90
ifort  %LINK_OPT% Sources/kernel/master/Reapro.f90
ifort  %LINK_OPT% Sources/kernel/master/restar.f90
ifort  %LINK_OPT% Sources/kernel/master/rrudat.f90
ifort  %LINK_OPT% Sources/kernel/master/runend.f90
ifort  %LINK_OPT% Sources/kernel/master/setgts.f90
ifort  %LINK_OPT% Sources/kernel/master/Turnof.f90
ifort  %LINK_OPT% Sources/kernel/master/Turnon.f90

ifort  %LINK_OPT% Sources/kernel/mathru/btdbma.f90
ifort  %LINK_OPT% Sources/kernel/mathru/frivel.f90
ifort  %LINK_OPT% Sources/kernel/mathru/funcrd.f90
ifort  %LINK_OPT% Sources/kernel/mathru/funcre.f90
ifort  %LINK_OPT% Sources/kernel/mathru/invert.f90
ifort  %LINK_OPT% Sources/kernel/mathru/invmtx.f90
ifort  %LINK_OPT% Sources/kernel/mathru/matove.f90
ifort  %LINK_OPT% Sources/kernel/mathru/mbmab0.f90
ifort  %LINK_OPT% Sources/kernel/mathru/mbmabt.f90
ifort  %LINK_OPT% Sources/kernel/mathru/mbmatb.f90
ifort  %LINK_OPT% Sources/kernel/mathru/mbvab0.f90
ifort  %LINK_OPT% Sources/kernel/mathru/mbvab1.f90
ifort  %LINK_OPT% Sources/kernel/mathru/mbvab2.f90
ifort  %LINK_OPT% Sources/kernel/mathru/ordena.f90
ifort  %LINK_OPT% Sources/kernel/mathru/rotnod.f90
ifort  %LINK_OPT% Sources/kernel/mathru/sortin.f90
ifort  %LINK_OPT% Sources/kernel/mathru/vecasi.f90
ifort  %LINK_OPT% Sources/kernel/mathru/vecnor.f90
ifort  %LINK_OPT% Sources/kernel/mathru/vecpro.f90
ifort  %LINK_OPT% Sources/kernel/mathru/vecres.f90
ifort  %LINK_OPT% Sources/kernel/mathru/vecuni.f90
ifort  %LINK_OPT% Sources/kernel/mathru/vetoma.f90

ifort  %LINK_OPT% Sources/kernel/outrut/chanum.f90                                                     
ifort  %LINK_OPT% Sources/kernel/outrut/geofem.f90
ifort  %LINK_OPT% Sources/kernel/outrut/geogid.f90
ifort  %LINK_OPT% Sources/kernel/outrut/nulgid.f90
ifort  %LINK_OPT% Sources/kernel/outrut/strcub.f90
ifort  %LINK_OPT% Sources/kernel/outrut/stream.f90
ifort  %LINK_OPT% Sources/kernel/outrut/strloc.f90
ifort  %LINK_OPT% Sources/kernel/outrut/suplot.f90

ifort  %LINK_OPT% Sources/kernel/soldir/fvecdo.f90
ifort  %LINK_OPT% Sources/kernel/soldir/rencon.f90
ifort  %LINK_OPT% Sources/kernel/soldir/renum0.f90
ifort  %LINK_OPT% Sources/kernel/soldir/renum1.f90
ifort  %LINK_OPT% Sources/kernel/soldir/renum2.f90
ifort  %LINK_OPT% Sources/kernel/soldir/renum3.f90
ifort  %LINK_OPT% Sources/kernel/soldir/renum4.f90
ifort  %LINK_OPT% Sources/kernel/soldir/renumb.f90
ifort  %LINK_OPT% Sources/kernel/soldir/renumn.f90
ifort  %LINK_OPT% Sources/kernel/soldir/skyase.f90
ifort  %LINK_OPT% Sources/kernel/soldir/skybak.f90
ifort  %LINK_OPT% Sources/kernel/soldir/skycek.f90
ifort  %LINK_OPT% Sources/kernel/soldir/skydia.f90
ifort  %LINK_OPT% Sources/kernel/soldir/skyini.f90
ifort  %LINK_OPT% Sources/kernel/soldir/skylin.f90
ifort  %LINK_OPT% Sources/kernel/soldir/skylpo.f90
ifort  %LINK_OPT% Sources/kernel/soldir/skyplu.f90
ifort  %LINK_OPT% Sources/kernel/soldir/skyren.f90
ifort  %LINK_OPT% Sources/kernel/soldir/skytri.f90
ifort  %LINK_OPT% Sources/kernel/soldir/Soldir.f90

rem ifort  %LINK_OPT% Sources/kernel/solmum/inimum.f90
rem ifort  %LINK_OPT% Sources/kernel/solmum/matspr.f90
ifort  %LINK_OPT% Sources/kernel/solmum/Solmum.f90

ifort  %LINK_OPT% Sources/kernel/wtools/codfix.f90
ifort  %LINK_OPT% Sources/kernel/wtools/cputim.f90
ifort  %LINK_OPT% Sources/kernel/wtools/listen.f90
ifort  %LINK_OPT% Sources/kernel/wtools/matrea.f90
ifort  %LINK_OPT% Sources/kernel/wtools/vecrea.f90
   
ifort  %LINK_OPT% Sources/kernel/solite/abxblo.f  
ifort  %LINK_OPT% Sources/kernel/solite/amuxso.f  
ifort  %LINK_OPT% Sources/kernel/solite/atmuxs.f  
ifort  %LINK_OPT% Sources/kernel/solite/auxblo.f  
ifort  %LINK_OPT% Sources/kernel/solite/auxdri.f  
ifort  %LINK_OPT% Sources/kernel/solite/bcgsol.f
ifort  %LINK_OPT% Sources/kernel/solite/bcgsta.f
ifort  %LINK_OPT% Sources/kernel/solite/bisblo.f
ifort  %LINK_OPT% Sources/kernel/solite/bisini.f
ifort  %LINK_OPT% Sources/kernel/solite/brkdns.f
ifort  %LINK_OPT% Sources/kernel/solite/cgnrso.f
ifort  %LINK_OPT% Sources/kernel/solite/cgpilu.f
ifort  %LINK_OPT% Sources/kernel/solite/cgsblo.f
ifort  %LINK_OPT% Sources/kernel/solite/cgsolv.f
ifort  %LINK_OPT% Sources/kernel/solite/csrase.f90
ifort  %LINK_OPT% Sources/kernel/solite/dbcgso.f
ifort  %LINK_OPT% Sources/kernel/solite/diablo.f
ifort  %LINK_OPT% Sources/kernel/solite/diapre.f
ifort  %LINK_OPT% Sources/kernel/solite/dqgmre.f
ifort  %LINK_OPT% Sources/kernel/solite/fadmem.f
ifort  %LINK_OPT% Sources/kernel/solite/fcputi.f
ifort  %LINK_OPT% Sources/kernel/solite/fgmres.f
ifort  %LINK_OPT% Sources/kernel/solite/flvecz.f
ifort  %LINK_OPT% Sources/kernel/solite/fomsol.f
ifort  %LINK_OPT% Sources/kernel/solite/frunei.f
ifort  %LINK_OPT% Sources/kernel/solite/frunen.f
ifort  %LINK_OPT% Sources/kernel/solite/frveca.f
ifort  %LINK_OPT% Sources/kernel/solite/fveccl.f
ifort  %LINK_OPT% Sources/kernel/solite/fvecze.f
ifort  %LINK_OPT% Sources/kernel/solite/givens.f
ifort  %LINK_OPT% Sources/kernel/solite/gmrblo.f
ifort  %LINK_OPT% Sources/kernel/solite/gmress.f
ifort  %LINK_OPT% Sources/kernel/solite/ilutso.f
ifort  %LINK_OPT% Sources/kernel/solite/implus.f
ifort  %LINK_OPT% Sources/kernel/solite/itdblo.f
ifort  %LINK_OPT% Sources/kernel/solite/itdriv.f
ifort  %LINK_OPT% Sources/kernel/solite/lusols.f
ifort  %LINK_OPT% Sources/kernel/solite/lutsol.f
ifort  %LINK_OPT% Sources/kernel/solite/maxmem.f
ifort  %LINK_OPT% Sources/kernel/solite/mgsros.f
ifort  %LINK_OPT% Sources/kernel/solite/qsplit.f
ifort  %LINK_OPT% Sources/kernel/solite/reduit.f
ifort  %LINK_OPT% Sources/kernel/solite/Solite.f
ifort  %LINK_OPT% Sources/kernel/solite/stopbi.f
ifort  %LINK_OPT% Sources/kernel/solite/tfqmrs.f
ifort  %LINK_OPT% Sources/kernel/solite/tidycg.f
ifort  %LINK_OPT% Sources/kernel/solite/uppdir.f
ifort  %LINK_OPT% Sources/kernel/solite/veiasi.f
       
echo ----------------------------------------------------------    
echo Compiling Alefor
echo ----------------------------------------------------------    
                                                                                                  
ifort  %LINK_OPT% Sources/modules/alefor/Alefor.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_begite.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_begste.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_concou.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_doiter.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_endste.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_inirot.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_iniunk.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_memall.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_openfi.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_outinf.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_output.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_reabcs.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_reaous.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_reaphy.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_reapro.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_rotate.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_turnof.f90
ifort  %LINK_OPT% Sources/modules/alefor/ale_turnon.f90

echo ----------------------------------------------------------    
echo Compiling Codire
echo ----------------------------------------------------------   

ifort  %LINK_OPT% Sources/modules/codire/Codire.f90   

echo ----------------------------------------------------------    
echo Compiling Nstinc
echo ----------------------------------------------------------    
                                                                                                  
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_assmat.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_assrhs.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_bounod.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_elmcon.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_elmdip.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_elmdir.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_elmesv.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_elmgat.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_elmlap.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_elmorp.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_elmpre.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_elmsto.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_eosvel.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_matfve.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_matpre.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_ortpro.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_prebcs.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsf_solite.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_autbcs.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_bcntoe.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_begite.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_begste.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_bounod.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_bouset.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_codfix.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_conblk.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_concou.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_corner.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_cvgunk.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_doiter.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_elmtss.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_endite.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_endste.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_inisol.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_iniunk.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_intbcs.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_intrst.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_memall.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_memsol.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_newmsh.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_openda.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_openfi.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_outbcs.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_outerr.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_outinf.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_output.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_outset.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_outtan.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_promsh.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_reabcs.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_reanut.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_reaous.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_reaphy.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_reapro.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_restar.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_tistep.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_turnof.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_turnon.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_updbcs.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_updfor.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_updmsh.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_updtss.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsi_updunk.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_assmat.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_assrhs.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_bouave.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_bougat.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_bouopb.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_bouope.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_bouwal.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmchl.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmdir.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmexf.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmgat.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmmul.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmope.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmpre.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmpro.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmrc1.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmrc2.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmrco.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmrcp.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmrdi.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmres.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmrtr.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmrvp.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmsch.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmsto.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmtem.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmtes.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_elmvsg.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_ifconf.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_matrix.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_rotmat.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_rotunk.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsm_solite.f90
ifort  %LINK_OPT% Sources/modules/nstinc/nsp_solite.f90
ifort  %LINK_OPT% Sources/modules/nstinc/Nstinc.f90

echo ----------------------------------------------------------    
echo Compiling Temper
echo ----------------------------------------------------------    
                                                                                                  
ifort  %LINK_OPT% Sources/modules/temper/tem_assmat.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_assrhs.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_bcntoe.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_begite.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_begste.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_bouope.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_bouset.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_bouwal.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_conblk.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_concou.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_cvgunk.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_doiter.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmchl.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmdif.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmdir.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmgat.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmope.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmpro.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmres.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmset.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmshc.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmtes.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmtss.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_elmvel.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_endite.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_endste.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_exaerr.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_exasol.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_inisol.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_iniunk.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_matrix.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_memall.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_newmsh.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_openfi.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_outbcs.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_outerr.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_outhfl.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_outinf.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_output.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_outset.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_radpos.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_radvuf.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_reabcs.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_reanut.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_reaous.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_reaphy.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_reapro.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_restar.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_solite.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_tistep.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_turnof.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_turnon.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_updbcs.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_updtss.f90
ifort  %LINK_OPT% Sources/modules/temper/tem_updunk.f90
ifort  %LINK_OPT% Sources/modules/temper/Temper.f90

echo ----------------------------------------------------------    
echo Compiling Turbul
echo ----------------------------------------------------------    
                                                                                                  
ifort  %LINK_OPT% Sources/modules/turbul/tur_addarr.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_assmat.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_assrhs.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_begite.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_begste.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_conblk.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_concou.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_cvgunk.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_doiter.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_elmchl.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_elmdir.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_elmfu1.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_elmga1.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_elmgl1.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_elmop1.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_elmpr1.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_elmsg1.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_elmts1.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_endite.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_endste.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_inisol.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_iniunk.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_matrix.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_memall.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_newmsh.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_openfi.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_outerr.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_outinf.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_output.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_reabcs.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_reanut.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_reaous.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_reaphy.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_reapro.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_restar.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_solite.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_tistep.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_turnof.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_turnon.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_updbcs.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_updedd.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_updtss.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_updunk.f90
ifort  %LINK_OPT% Sources/modules/turbul/tur_updwbc.f90
ifort  %LINK_OPT% Sources/modules/turbul/Turbul.f90

echo ----------------------------------------------------------    
echo Linking objects and creating Alya.exe
echo ----------------------------------------------------------    
                          
ifort /exe:Alya.exe *.obj

:end
