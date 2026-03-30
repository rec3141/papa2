# Makefile for building standalone CPU-only libpapa2.so (Python-compatible)

CXX ?= g++
CC ?= gcc
CXXFLAGS += -fPIC -std=c++11 -DNO_RCPP -DNDEBUG -fopenmp -Wno-format
CFLAGS += -fPIC -DNDEBUG

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

.PHONY: all clean docs

all: libpapa2.so

COBJS = $(CSRCS:.c=.o)

libpapa2.so: $(OBJS) $(COBJS)
	$(CXX) -shared -o $@ $^ $(LDFLAGS) -fopenmp -lz

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f src/*.o libpapa2.so

# Build documentation with mkdocs-material
docs:
	mkdocs build

# Print config
info:
	@echo "CXX: $(CXX)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
