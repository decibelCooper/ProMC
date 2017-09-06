# S.Chekanov, D. Blyth


ifndef PROMC
$(error PROMC env variable is not set. Run setup.sh first)
endif

ifndef ROOTSYS 
$(error ROOTSYS env variable is not set. Install ROOT first)
endif


include ${PROMC}/etc/config.mk
#include ${ROOTSYS}/etc/Makefile.arch

# Root variables
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGTTLIBS   = $(shell root-config --glibs)
# Assign or add variables
CXXFLAGS     += $(ROOTCFLAGS)
LIBS         += $(ROOTLIBS) -lEG
LIBS         += -L$(PROMC)/lib -lpromc -lprotoc -lprotobuf -lprotobuf-lite -lcbook -lz


INCLUDE1= -I./inc -I./
INCLUDE2= -I./src
INCLUDE3= -I$(PROMC)/include -I$(PROMC)/src

Tasks:     clean pcsim2promc

SOURCE_FILES := pcsim2promc.cpp
objects       = $(patsubst %.cpp,%.o,$(SOURCE_FILES))

%.o: %.cpp
	$(CXX) $(OPT) $(CXXFLAGS) $(INCLUDE1) $(INCLUDE2) $(INCLUDE3) -o $@ -c $<

pcsim2promc: $(objects)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@
clean:
	        @rm -f *.o *~ pcsim2promc; echo "Clear.." 
