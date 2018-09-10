#

# Sub directories containing source code, except for the main programs
SUBDIRS := include src include/readpaf include/rapidjson/include include/zstr/src

#
# Set libraries, paths, flags and options
#

#Basic flags every build needs
LIBS = -lz -llzma -lbz2  -lpthread
CXXFLAGS ?= -g -O3 -Wall -Wextra
CXXFLAGS += -std=c++11 
CFLAGS ?= -O3 -std=c99
CXX ?= g++
CC ?= gcc

# Change the value of HDF5, EIGEN, or HTS below to any value to disable compilation of bundled code
HTS?=install

# Default to build and link the libhts submodule
ifeq ($(HTS), install)
    HTS_LIB=./include/htslib/libhts.a
    HTS_INCLUDE=-I./include/htslib -I./include/htslib/htslib
else
    # Use system-wide htslib
    HTS_LIB=
    HTS_INCLUDE=
    LIBS += -lhts
endif

# Include the src subdirectories
NP_INCLUDE=$(addprefix -I./, $(SUBDIRS))

# Add include flags
CPPFLAGS += $(NP_INCLUDE) $(HTS_INCLUDE)

# Main programs to build
PROGRAM=preqclr

all: $(PROGRAM)

#
# Build libhts
#
htslib/libhts.a:
	cd include/htslib && make || exit 255


#
# Source files
#

# Find the source files by searching subdirectories
CPP_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.cpp))
C_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.c))
EXE_SRC=./src/main/preqclr.cpp

# Automatically generated object names
CPP_OBJ=$(CPP_SRC:.cpp=.o)
C_OBJ=$(C_SRC:.c=.o)

# Compile objects
.cpp.o:
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) -fPIC $<

.c.o:
	$(CC) -o $@ -c $(CFLAGS) $(CPPFLAGS) -fPIC $<

# Link main executable
$(PROGRAM): ./src/main/preqclr.o $(CPP_OBJ) $(C_OBJ) 
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) -fPIC $< $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(LIBS) $(LDFLAGS)


clean:
	rm -f $(PROGRAM) $(CPP_OBJ) $(C_OBJ) src/main/preqclr.o src/sequence.o
