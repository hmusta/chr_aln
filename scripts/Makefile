CXX ?= g++
CXXFLAGS ?= -std=c++20 -O3 -Wall -Wextra

BINARIES := compute_edit_distance syri_to_cigar

all: $(BINARIES)

%: %.cpp helpers.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(BINARIES)

.PHONY: all clean
