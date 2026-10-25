# AMPT v1.26t9b/v2.26t9b + open-charm framework (charm_mod.f, coal_new.f).  Same flags as the stock Makefile's
# gfortran build; LHAPDF 6 (CTEQ6.1M) and the official EPS09s code are linked in.  scripts/build.sh sets the paths.
FC = gfortran
FFLAGS = -O2 -std=legacy -ffixed-line-length-132
CXX = g++
LHAPDF_PREFIX ?= ../../deps/lhapdf
EPS09S_DIR ?= ../../deps/eps09s
LIBS = -L$(LHAPDF_PREFIX)/lib -lLHAPDF -lstdc++ -Wl,-rpath,$(LHAPDF_PREFIX)/lib

OBJS = amptsub.o art1f.o hijing1.383_ampt.o hipyset1.35.o zpc.o main.o linana.o charm_mod.o coal_new.o eps09s_f.o

all: ampt

%.o: %.f hqctl.inc
	$(FC) $(FFLAGS) -c $< -o $@

eps09s_f.o: eps09s_f.cpp $(EPS09S_DIR)/eps09s.cpp
	$(CXX) -O2 -I$(EPS09S_DIR) -c $< -o $@

ampt: $(OBJS)
	$(FC) $(FFLAGS) -o ampt $(OBJS) $(LIBS)

clean:
	rm -f *.o ampt

.PHONY: all clean
