OBJ=  f2allm.o f2pyth.o fdilut.o gmc_random.o mkasym.o mkf2.o mkr.o pyth_xsec.o pythia6220.o r1990.o r1998.o radgen.o radgen_event.o radgen_init.o main.o mycernlib.o  

# -lstdc++ for Gagik's stuff
#
FOR   =  -lstdc++ -lnsl -lcrypt -ldl -lgfortran
CERNLIBS =  -L/apps/cernlib/x86_64_rhel6_4.7.2/2005/lib -lmathlib  -lpacklib
claspyth : $(OBJ)
	 gfortran  -o -fno-automatic -ffixed-line-length-none -fno-second-underscore -o claspyth $(OBJ) $(FOR)  
$(OBJ) : %.o: %.F
	gfortran   -DLinux -fno-automatic -ffixed-line-length-none -fno-second-underscore  -I. -I./include -c $< -o $@  

clean:
	rm -f claspyth   $(OBJ)



##g77  -O2 -fno-automatic -finit-local-zero -ffixed-line-length-none -fno-second-underscore \
##        -DLinux \
##        -I. -I./ -I/group/clas/builds/release-4-14/packages/include -I/group/clas/builds/release-4-14/packages/inc_derived -I/apps/tcl/include -I/usr/X11R6/include -c \
##        aac.F  -o /home/avakian/w6/tmp/obj/LinuxRHEL3/generator/aac.o


















