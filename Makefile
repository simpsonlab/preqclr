#

# Sub directories containing source code, except for the main programs
SUBDIRS := include src include/readpaf include/rapidjson/include include/zstr/src include/spoa 

#
# Set libraries, paths, flags and options
#

#Basic flags every build needs
LIBS = -lz -Llib/ -lspoa
CXXFLAGS ?= -g -O3
CXXFLAGS += -std=c++11 
CFLAGS ?= -O3 -std=c99
CXX ?= g++
CC ?= gcc

# Include the src subdirectories
NP_INCLUDE=$(addprefix -I./, $(SUBDIRS))

# Add include flags
CPPFLAGS += $(NP_INCLUDE)

# Main programs to build
PROGRAM=preqclr

all: $(PROGRAM)

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
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) -fPIC $< $(CPP_OBJ) $(C_OBJ) $(LIBS) $(LDFLAGS)


clean:
	rm -f $(PROGRAM) $(CPP_OBJ) $(C_OBJ) src/main/preqclr.o src/sequence.o
