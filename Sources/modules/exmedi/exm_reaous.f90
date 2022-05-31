subroutine exm_reaous
   !-----------------------------------------------------------------------
   !****f* Exmedi/exm_reaous
   ! NAME 
   !    exm_reaous
   ! DESCRIPTION
   !    This routine reads the output strategy 
   ! USES
   !    listen
   ! USED BY
   !    exm_turnon
   !***
   !-----------------------------------------------------------------------
   use      def_parame
   use      def_inpout
   use      def_master
   use      def_domain
   use      mod_memchk

   use      def_exmedi
   use mod_ecoute, only :  ecoute

   implicit none
   integer(ip) idumy
   
   if( INOTSLAVE ) then
      
      !
      ! Reach the section
      !
      rewind(lisda)
      call ecoute('exm_reaous')
      do while(words(1)/='OUTPU')
         call ecoute('exm_reaous')
      end do
      !
      ! Begin to read data
      !
      do while(words(1)/='ENDOU')
         call ecoute('exm_reaous')

         call posdef(2_ip,idumy)
         
         if(words(1)=='PARAM') then
         ! Postprocess general parameters
            call ecoute('exm_reaous')
            do while(words(1)/='ENDPA')
               if(words(1)=='ISOCH') then
                  call runend('EXM_REAOUS: ISOCHRONES THRESHOLD MUST BE WRITTEN IN PROPERTIES SECTION')
               end if
               call ecoute('exm_reaous')
            end do
         
         else if (words(1) == "SAVEC") then
            if (exists("CELLM")) then
               ! SAVE_CONVERGENCE CELL_MODEL -- saves ohara.mXXcXX files with curves for different 
               !concentrations during the ODE iterations for each material and cell type
               kfl_save_convergence_cellmodel = 1_ip 
            end if

         else if (words(1) == "DUMPI") then
            ! DUMP_INITIAL_CONDITIONS -- saves initial conditions to be hardcoded for a cell model 
            !concentrations during the ODE iterations for each material and cell type
            kfl_save_init_cellmodel = 1_ip 

         end if
      end do
   end if

end subroutine exm_reaous
       
   