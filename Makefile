default: native

# Identify all directory locations
EMP_DIR = ../Empirical

TARGET = Mutation-Study

CXX = c++

# Specify sets of compilation flags to use
FLAGS_version = -std=c++23
FLAGS_warn    = -Wall -Wextra -Wno-unused-function -Woverloaded-virtual -pedantic
FLAGS_include = -I$(EMP_DIR)/include/
FLAGS_main    = $(FLAGS_version) $(FLAGS_warn) $(FLAGS_include) -fopenmp

FLAGS_OPT     = $(FLAGS_main) -O3 -DNDEBUG 
FLAGS_DEBUG  = $(FLAGS_main) -g -DEMP_TRACK_MEM # EMP_TRACK_MEM

native: FLAGS = $(FLAGS_OPT)
native: $(TARGET)

debug: FLAGS = $(FLAGS_DEBUG)
debug: $(TARGET)

$(TARGET): main.cpp
	$(CXX) $(FLAGS) main.cpp -o $(TARGET)

# 'make new' runs 'clean' then 'native'
new: clean native

clean:
	rm -rf $(TARGET)