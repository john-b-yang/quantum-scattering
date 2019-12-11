CXX=clang++
CFLAGS := -O3
CFLAGS += -std=c++11

odir = ./o.d
ddir = ./d.d
fortObjDir = o.fortObj

FC = gfortran

SRCFILES:= $(notdir $(wildcard *.cpp))
FORTfSRC := $(notdir $(wildcard *.f))
FORT90SRC := $(notdir $(wildcard *.f90))
OBJS:=$(patsubst %.cpp,$(odir)/%.o,$(SRCFILES))
FORTOBJSF=$(patsubst %.f,$(fortObjDir)/%.o,$(FORTfSRC))
FORTOBJSF90=$(patsubst %.f90,$(fortObjDir)/%.o,$(FORT90SRC))
FORTOBJS = $(FORTOBJSF) $(FORTOBJSF90)

DEPS:=$(patsubst $(odir)/%.o,$(ddir)/%.d, $(OBJS))

$(odir)/%.o:%.cpp GNUmakefile
	mkdir -p $(odir);$(CXX) -c $(CFLAGS) $< -o $@
	mkdir -p $(ddir);$(CXX) -MM $(CFLAGS) $< | sed '1s/^/o.d\//' > $*.d;mv $*.d $(ddir)

$(fortObjDir)/%.o:%.f90 GNUmakefile
	mkdir -p $(fortObjDir);$(FC) -c -O3 $< -o $@

$(fortObjDir)/%.o:%.f GNUmakefile
	mkdir -p $(fortObjDir);$(FC) -c -O3 $< -o $@

testFort:GNUmakefile $(FORTOBJS) $(OBJS) 
	$(CXX) $(CFLAGS) $(FORTOBJS) $(OBJS) -lgfortran -o testFort.exe

-include $(DEPS)

clean:
	rm -r *.exe $(odir) $(ddir) $(fortObjDir)

listsrc:
	@echo $(SRCFILES) $(FORTSRC)
listobj:
	@echo $(OBJS)
listdep:
	@echo $(DEPS)
listfort:
	@echo $(FORTSRC)
	@echo $(FORTOBJS)
