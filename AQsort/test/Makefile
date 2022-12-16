CXX = g++
CXXFLAGS = -O0 -g 
#CXXFLAGS = -O3 -DNDEBUG
#CXXFLAGS +=-fopenmp
CXXFLAGS += -I../include
CXXFLAGS += -std=c++11

PROGRAMS = spm_fo spm_lambda

.PHONY: all
all: $(PROGRAMS)

$(PROGRAMS): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -f $(PROGRAMS)
