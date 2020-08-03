CXX=g++
CXXFLAGS=-fopenmp -fexceptions
LDFLAGS=-fopenmp -lboost_program_options

BINDIR=bin
build_folder := $(shell mkdir -p $(BINDIR))

RELEASE_OBJS=HMMIO.o commandLineParser.o dataLoad.o getPrH.o lnsnum.o main.o master.o masterUtils.o options.o processChr.o reestimator.o stateLD.o stateTemplate.o test.o viterbi.o
DEBUG_OBJS=HMMIO.o commandLineParser.o dataLoad.o getPrH.o lnsnum.o main.o master.o masterUtils.o options.o processChr.o reestimator.o stateLD.o stateTemplate.o test.o viterbi.o
OUTPUT_RESULTS_OBJS=HMMIO.o dataLoad.o lnsnum.o masterUtils.o options.o outputResults.o stateLD.o stateTemplate.o viterbi.o writeResults.o

CFLAGS_DEBUG=-g
CFLAGS_RELEASE=-fomit-frame-pointer -O3 -ffast-math -floop-interchange -floop-strip-mine -floop-block -funroll-loops
CFLAGS_PRECOMPILE=-O2

all: release output_results

release: CFLAGS_EXTRA=$(CFLAGS_RELEASE)
release: $(RELEASE_OBJS)
	$(CXX) $(LDFLAGS) $(RELEASE_OBJS) -o $(BINDIR)/NucHMM-learn
	strip $(BINDIR)/NucHMM-learn

debug: CFLAGS_EXTRA=$(CFLAGS_DEBUG)
debug: $(DEBUG_OBJS)
	$(CXX) $(LDFLAGS) $(DEBUG_OBJS) -o $(BINDIR)/NucHMM-learn-debug

output_results: $(OUTPUT_RESULTS_OBJS)
	$(CXX) $(LDFLAGS) $(OUTPUT_RESULTS_OBJS) -o $(BINDIR)/NucHMM-output_results
	strip $(BINDIR)/NucHMM-output_results

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CFLAGS_EXTRA) -o $@ -c $<

clean: clean-objs
	rm -f $(BINDIR)/*

clean-objs:
	rm -f *.o

