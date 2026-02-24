# Makefile

PYTHON = python

.PHONY : default all format clean

default : all

all : gauss_legendre.txt

gauss_legendre.txt: gauss_legendre_quadrature.py
	$(PYTHON) gauss_legendre_quadrature.py

format :
	black gauss_legendre_quadrature.py

clean:
	$(RM) gauss_legendre.txt legendre.png
