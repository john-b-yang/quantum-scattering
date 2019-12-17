HOME = ./..

include $(HOME)/GNUmakefile
CFLAGS += -O3 
FFLAGS = -lgfortran
# CFLAGS += -g 

testQuad: $(CPPOBJS) $(FORTOBJS) testQuad.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FFLAGS) -o testQuad.exe $^

clean:
	rm *.exe

-include $(OBJS:.o=.d)
