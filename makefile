CFLAGS = -m64 -fopenmp -O3 -std=c++11 -w

## Comment the line below if you want to use your own CUBA version 
#WITH_CUBA=true

CUBA_HEADER_DIR=/usr/local/include
CUBA_LIBRARY_DIR=/usr/local/include
CUBA_LIBFLAG=-lcuba
CUBA_EXEC=

ifdef WITH_CUBA
CUBA_HEADER_DIR=$(PWD)/Cuba
CUBA_LIBRARY_DIR=$(PWD)/Cuba
CUBA_LIBFLAG=$(PWD)/Cuba/libcuba.a
CUBA_EXEC=cd Cuba && ./configure --prefix=$(PWD)/ && make 
endif

LHAPDF_HEADER_DIR=/usr/local/include
LHAPDF_LIBRARY_DIR=/usr/local/lib

XS_FILES =$(patsubst ./XS/%.cpp, ./XS/%.o, $(wildcard ./XS/*.cpp))
SOURCE =$(patsubst ./src/%.cpp, ./src/%.o, $(wildcard ./src/*.cpp))


LFLAGS = -lm -lLHAPDF $(CUBA_LIBFLAG)
LIB_FLAGS = -L$(LHAPDF_LIBRARY_DIR) -L$(CUBA_LIBRARY_DIR)
INC_FLAGS = -I$(LHAPDF_HEADER_DIR) -I$(CUBA_HEADER_DIR)


all: main

main : % : %.o $(SOURCE) $(XS_FILES)
	$(CUBA_EXEC)
	$(CXX) $(CFLAGS) $(SOURCE) $(XS_FILES) $<  -o $@ $(INC_FLAGS) $(LIB_FLAGS) $(LFLAGS)
	rm $@.o


%.o : %.cpp
	$(CXX) $(CFLAGS) $(INC_FLAGS) -c $< -o $@


clean:
	rm	-rf *.out *.o *.dSYM ./src/*.o ./XS/*.o main Cuba/libcuba.a
	cd Cuba && make clean

cleansrc:
	rm 	-rf *.out *.o ./src/*.o  *.dSYM
