# This makefile is designed to compile and link Fortran code 
# for the 2D simulation in directory Flow2d/Code
# 
# Notes 1) This version allows for the use of dbx to control 
#          program during execution for debugging. 
#       2) The executable code is called `main'.
#       3) The mainline should be in the file main.f

# Set Fortran compiler flags for STOKES
FFLAGS = -ansi -c -C

# List of object files to be linked:

LIST  = main.o\
acosh.o\
adcont.o\
addisp.o\
bccalc.o\
bcelim.o\
bndc.o\
cargrd.o\
coeff.o\
cfppty.o\
conv.o\
crst.o\
dhat.o\
difphi.o\
factr.o\
gradp.o\
grdgeo.o\
hoconv.o\
inital.o\
input.o\
makgrd.o\
masflx.o\
meana.o\
null.o\
pptys.o\
resid.o\
sgen.o\
src.o\
srcuvp.o\
tres.o\
type.o\
uhat.o\
uvptrs.o\
vhat.o\
weight.o\
what.o\
assign.o\
solt.o\
soluvp.o\
WATsitB/globm.o\
WATsitB/watsoln.o\
WATsitB/d1mach.o\
WATsitB/dblas1.o\
WATsitB/watsitB.o


# Rule for making an executable called `main'.
main : $(LIST)  
	f77 $(LIST) -o main

# Rule for creating an object file from its dependent fortran file.
%.o : %.f  
	f77 $(FFLAGS) $*.f -o $*.o
#
#********************************************************************
# The Fortran check task

# This stage illustrates how to use ftnchek to check the sample code
# that was compiled above. This version sends output for each file to
# its own .lis file and sends global messages to check.lis. Note that
# routines inside libraries are NOT checked.

# The Fortran checking is invoked by issuing the command:
# make check.lis
# It is assumed that the makefile is stored in the file `makefile'

#******************************************************************
# Template for check.lis

# Add project files to the make suffixes
.SUFFIXES: .prj

# Create the list of project filenames from the object file list
LISTP= $(LIST:.o=.prj)

# Rule for checking complete program for global errors
check.lis: $(LISTP)
	/opt/Ftnchek/ftnchek-2.10.1/ftnchek -output=check.lis $(LISTP) 

# Rule for creating a .prj file from its dependent .f file
.f.prj:
	/opt/Ftnchek/ftnchek-2.10.1/ftnchek -declare -noextern -library -project \
-output=$*.lis $<


