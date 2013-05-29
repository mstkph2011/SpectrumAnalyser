#vpath %.h ../
#vpath %.cxx ../
#vpath %.h $(ROOTSYS)/include

INCLFLAGS  = -I$(ROOTSYS)/include -I../
LDLIBS     = $(shell root-config --glibs) \
                -L/$(ROOTSYS)/lib -lMathCore -lSpectrum 

#LDLIBS     = $(shell root-config --glibs)
# -lMinuit


SOURCES = $(wildcard *.cxx)
HEADERS = $(wildcard *.h)
OBJECTS = $(patsubst %.cxx,%.o,$(wildcard *.cxx))

CC = g++
CFLAGS = -O2 -DNDEBUG -Wall


all:                    SpecTest

SpecTest: SpecTest.cxx THypGeSpectrumAnalyser.o THypGePeakFitFunction.o
	$(CC) $(CFLAGS) $(INCLFLAGS) $(LDLIBS) SpecTest.cxx THypGeSpectrumAnalyser.o THypGePeakFitFunction.o $(LDLIBS) -o SpecTest

THypGeSpectrumAnalyser.o: THypGeSpectrumAnalyser.cxx THypGeSpectrumAnalyser.h
	$(CC) $(CFLAGS) $(INCLFLAGS) -c THypGeSpectrumAnalyser.cxx

THypGePeakFitFunction.o: THypGePeakFitFunction.cxx THypGePeakFitFunction.h
	$(CC) $(CFLAGS) $(INCLFLAGS) -c THypGePeakFitFunction.cxx
.PHONY :                clean

clean :
			rm $(OBJECTS) *.o *.d.* SpecTest

