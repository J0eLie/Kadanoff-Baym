CC = icpc
CFLAGS = -diag-disable=10441 -O3 -std=c++11
#-D_DEBUG
#OMP = -qopenmp
OMP =

INC = -I$(MKLROOT)/include -I$(MKLROOT)/include/fftw/include

LIB = -L$(MKLROOT)/lib -qmkl -L$(MKLROOT)/interfaces/fftw3xc -lfftw3xc_intel 

SRC = grid.cpp \
			gf.cpp \
			coulomb.cpp \
			sigma.cpp \
			collision.cpp \
			propagate.cpp \
			output.cpp \
			fft.cpp

OBJ = $(SRC:%.cpp=%.o)

all : KadanoffBaym

KadanoffBaym : main.o $(OBJ)
	$(CC) $(CFLAGS) $(OMP) $(LIB) $(^F) -o $@

%.o : %.cpp
	$(CC) $(CFLAGS) $(OMP) $(INC) -c $(<F) -o $(@F)

clean :
	rm -f KadanoffBaym *.o
