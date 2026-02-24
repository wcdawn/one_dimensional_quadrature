# Makefile

CXX = g++
CFLAGS = -Wall -Wextra -std=c++23 -pedantic -O3

OBJ = main.o \
	  quadrature1d.o \
	  quadrature_gauss_legendre.o

EXEC = quad.x

.PHONY : default all clean

default : all

all : $(EXEC)

$(EXEC) : $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^

%.o : %.cpp
	$(CXX) $(CFLAGS) -c $<

../gauss_legendre/gauss_legendre.txt : ../gauss_legendre/gauss_legendre_quadrature.py
	$(MAKE) -C ../gauss_legendre

main.o : quadrature1d.hpp quadrature_gauss_legendre.hpp
quadrature1d.o : quadrature1d.hpp
quadrature_gauss_legendre.o : quadrature1d.hpp quadrature_gauss_legendre.hpp ../gauss_legendre/gauss_legendre.txt

clean:
	$(MAKE) -C ../gauss_legendre clean
	$(RM) $(OBJ)
	$(RM) $(EXEC)
