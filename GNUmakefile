#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

SRC = $(HOME)/src
FORTDIR = $(HOME)/src/fortran

CXX = clang++
odir = $(SRC)/objects
ddir = $(SRC)/dependicies

FC = gfortran
fortObjDir = $(FORTDIR)/objects

VPATH = $(HOME) $(SRC) $(FORTDIR)
CXXFLAGS = -Wall -I$(SRC) -std=c++11

CPPSRC := $(notdir $(wildcard $(SRC)/*.cpp))
CPPOBJS := $(patsubst %.cpp, $(odir)/%.o, $(CPPSRC))

FORTFSRC := $(notdir $(wildcard $(FORTDIR)/*.f))
FORT90SRC := $(notdir $(wildcard $(FORTDIR)/*.f90))
FORTFOBJS := $(patsubst %.f, $(fortObjDir)/%.o, $(FORTFSRC))
FORT90OBJS := $(patsubst %.f90, $(fortObjDir)/%.o, $(FORT90SRC))
FORTOBJS = $(FORTFOBJS) $(FORT90OBJS)

<<<<<<< HEAD
# DEPS := $(patsubst $(odir)/%.o, $(ddir)/%.d, $(CPPOBJS))
=======
DEPS := $(patsubst $(odir)/%.o, $(ddir)/%.d, $(CPPOBJS))
>>>>>>> zachs_branch

$(odir)/%.o:%.cpp GNUmakefile
	mkdir -p $(odir);$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
	mkdir -p $(ddir);$(CXX) -MM $(CFLAGS) $< | sed '1s/^/o.d\//' > $*.d;mv $*.d $(ddir)

$(fortObjDir)/%.o:%.f GNUmakefile
	mkdir -p $(fortObjDir); $(FC) -c -O3 $< -o $@

$(fortObjDir)/%.o:%.f90 GNUmakefile
	mkdir -p $(fortObjDir); $(FC) -c -O3 $< -o $@

listsrc:
	@echo $(CPPSRC) $(FORTFSRC) $(FORT90SRC)

listcpp:
	@echo $(CPPSRC)
	@echo $(CPPOBJS)

listfort:
	@echo $(FORTFSRC) $(FORT90SRC)
	@echo $(FORTOBJS)

listdep:
	@echo $(DEPS)

<<<<<<< HEAD
clean:
	rm -r $(CPPOBJS) $(CPPOBJS:.o=.d) $(FORTOBJS) *.exe 
=======
cleantree:
	rm -r $(CPPOBJS) $(CPPOBJS:.o=.d) $(FORTOBJS) $(DEPS) 
>>>>>>> zachs_branch
