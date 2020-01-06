HOME = ./..

include $(HOME)/GNUmakefile
# CFLAGS += -O3 
FFLAGS = -lgfortran
CFLAGS += -g 

testFEMDVR: $(CPPOBJS) $(FORTOBJS) testFEMDVR.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FFLAGS) -o testFEMDVR.exe $^

-include $(OBJS:.o=.d)
