rem 
rem  Generate a pdf file form the Latex info file of Alya
rem
set name=%1
"C:\Archivos de programa\gnuplot\bin\wgnuplot.exe" %name%-latex.plt
rem pdflatex  %name%-latex.tex
latex %name%-latex.tex
latex %name%-latex.tex
dvips -o %name%-latex.ps %name%-latex
ps2pdf %name%-latex.ps 
