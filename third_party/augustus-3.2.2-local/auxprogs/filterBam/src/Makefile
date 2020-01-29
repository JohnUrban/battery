#	Makefile for filterBam
#
# 	Created: 11-November-2011
#	Last modified: 4-March-2012


SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleAlignment.cc \
			printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
# BAMTOOLS = 
INCLUDES = -I$(BAMTOOLS)/include/ -Iheaders -I./bamtools
LIBS = $(BAMTOOLS)/lib/libbamtools.a -lz
CFLAGS = -std=c++0x
VPATH = functions

all : $(PROGRAM) $(OBJECTS) CHECKBAM BIN

BIN : $(PROGRAM) CHECKBAM
	cp filterBam ../../../bin/filterBam

CHECKBAM:
	@if [ -z "${BAMTOOLS}" ]; then \
		echo '[Makefile]: $${BAMTOOLS} has not been set.'; \
		exit 101; \
	elif [ "${BAMTOOLS}" ]; then \
		echo "filterBam compiled with BAMTOOLS=${BAMTOOLS}"; \
	fi

$(PROGRAM) : $(OBJECTS)
	$(LINK.cc) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS) 

$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CFLAGS) $(CPPFLAGS) -c $^ -o $@ $(INCLUDES)


clean:
	rm -f *~ $(OBJECTS)
