#FILENAME: makefile

FC = gfortran

GDB = -g

OBJS = num_mod.o ang_mod.o dvr_mod.o build_dq_nw.o build_grids.o \
	dvreval.o gq.o quads.o setupfls.o build_KEmats.o get_grid_param.o \
	plot_dvr.o cdiag.o main.o

main: $(OBJS)
	$(FC)  -o  ../bin/xoneE $(OBJS) -llapack

.f90.o: 
	$(FC) -c $< 

.f.o:
	$(FC) -c $< 

.SUFFIXES: .f90

clean: 
	rm -f *~
	rm -f *.o
	rm -f *.mod
