-------------------------------------------------------------
-----------------------ALYA STANDARD COMPILATION------------
-------------------------------------------------------------
1.Go into the Executables/unix and select the config file corresponding to your machine from configure.in folder:
 - config_ifort.in (recomended for marenostrum 3)
 - config_gfortran.in
 - config_IBM_xlf.in

2.Copy the selected file with new name "config.in" into the unix folder:
"cp configure.in/config_ifort.in config.in"

3.Edit the file config.in, select the options you want, not all the options are compatible between them.
¡¡IMPORTANT!! Load the mandatory modules in the system for the selected options.

4.Run the configure command:
./configure [-g | -x] [all | module1  module2 ... | service1 service2 ...]
Example: ./configure -x nastin parall
¡¡IMPORTANT!! To run Alya in Parallel add the "parall" service in the configure comand.

5.If you want to use METIS run the command:
make metis4

6.Compile Alya
make -j num_processors

-------------------------------------------------------------
-----------------------ALYA MAKE OPTIONS---------------------
-------------------------------------------------------------

 - make help       print the help
 - make alya2pos   Compile alya2pos tool, to transform alya output in standard formats
 - make clean      Remove all binaries and objects directories
 - make cleandocs  Remove all alya documentation files
 - make docs       Creates alya documentation. Open in a browser /Documents/doxygen/html/index.html
 - make tutorial   Creates alya tutorial documentation.
 - make testSuite  Runs the testsuite locally with the configure in TestSuite/testSuite.xml
 - make libple     Compiles the libple library in Thirparties.
 - make libplepp   Compiles the libplepp library in Thirdparties.
 - make metis4     Compiles metis4 library in Thirparties.
 - make ninja      Compiles ninja library in Thirparties, it's a GPU solver.

-------------------------------------------------------------
-----------------------ALYA ONLINE DOCUMENTATION-------------
-------------------------------------------------------------
Overview:
http://bsccase02.bsc.es/alya/overview/

Documentation:
http://bsccase02.bsc.es/alya/

Tutorials:
http://bsccase02.bsc.es/alya/tutorial/index.html

Bibliography:
http://bsccase02.bsc.es/alya/overview/bibliography/bibliography.html

Gallery:
https://www.youtube.com/playlist?list=PLbABsNMD2jhy4CPd5iH_AkZRhmdMslfXD

