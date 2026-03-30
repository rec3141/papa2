# Makefile for building standalone CPU-only libdada2.so (Python-compatible)

CXX = g++
CXXFLAGS = -O3 -fPIC -std=c++11 -DNO_RCPP -DNDEBUG -march=native -fopenmp -Wno-format -flto

# Source files for standalone build (excludes Rmain.cpp, RcppExports.cpp,
# taxonomy.cpp, chimera.cpp, evaluate.cpp, filter.cpp)
CSRCS = src/derep.c src/loess.c

SRCS = src/dada2_capi.cpp \
       src/paired_capi.cpp \
       src/taxonomy_capi.cpp \
       src/chimera_capi.cpp \
       src/cluster.cpp \
       src/containers.cpp \
       src/pval.cpp \
       src/error.cpp \
       src/kmers.cpp \
       src/misc.cpp \
       src/nwalign_endsfree.cpp \
       src/nwalign_vectorized.cpp

OBJS = $(SRCS:.cpp=.o)
LDFLAGS = -lm

.PHONY: all clean

all: libdada2.so

COBJS = $(CSRCS:.c=.o)

libdada2.so: $(OBJS) $(COBJS)
	$(CXX) -shared -o $@ $^ $(LDFLAGS) -fopenmp -lz

src/%.o: src/%.c
	gcc -O3 -fPIC -DNDEBUG -march=native -c $< -o $@

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f src/*.o libdada2.so

# Print config
info:
	@echo "CXX: $(CXX)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
