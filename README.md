# One-Dimensional Quadrature

Basic investigations of one dimensional quadratures.

# Computing Quadarture Weights and Abscissae with Arbitrary Precision

I have often found it useful in my work to perform one-dimensional numerical integration with Gauss-Legendre quadrature series. Typically, it is necessary to pre-compute the quadrature weights and abscissae as these are the roots of high-order polynomials and cannot be expressed simply for quadrature orders N > 5 (approximately).

It is reasonably easy to find numerical values for these Gauss-Legendre weights and abscissae online (my favorite [here](https://pomax.github.io/bezierinfo/legendre-gauss.html)). However, I was curious how much work it would be to compute the points with arbitrary precision. I was able to find a method that works, but it is really quite slow.

# Method

Of course, the absciassae for the Gauss-Legendre quadrature for order N are the zeroes of the Nth Legendre polynomial. The weights then have a simple experssion as a function of the derivatives of the Legendre polynomials.

By using the [`mpmath`](https://mpmath.org/) library in Python, it is possible to perform numerical computations with arbitrary decimal precision (at a cost, of course). Using recursion relationships to define the Legendre polynomials, I was able to numerically evaluate them in arbitrary precision without loss of accuracy. However, the use of arbitrary precision in recursive Python functions makes for a very slow evaluation of these quantities.

One of the slowest (worst) parts of my solution algorithm for the zeroes of the Legendre polynomials is using a bisection search to find them. I begin by laying out a uniform grid and doing bisection searches within each interval until I find the correct number of zeros (known a priori).

# Results

I have tested this up to N = 16 Gauss-Legendre quadrature. It takes a few seconds to build the quadrature points. I then wrote a C++ function to test them out numerically. It is pretty easy to see the dobule-precision floating point numbers saturate quickly.

To try for yourself:
1. Clone the repository.
2. `cd src`
3. `make`
4. `./quad.x`

```
begin

ARGV:
./quad.x

n=1 int=2.0000000000000000e+00 err=-3.1705803038420699e-01
n=2 int=1.6758236553899863e+00 err=7.1183142258066656e-03
n=3 int=1.6830035477269170e+00 err=-6.1578111123949952e-05
n=4 int=1.6829416886959734e+00 err=2.8091981962852230e-07
n=5 int=1.6829419704071922e+00 err=-7.9139916842052571e-10
n=6 int=1.6829419696142791e+00 err=1.5139001163788635e-12
n=7 int=1.6829419696157950e+00 err=-1.9984014443252818e-15
n=8 int=1.6829419696157930e+00 err=0.0000000000000000e+00
n=9 int=1.6829419696157930e+00 err=0.0000000000000000e+00
n=10 int=1.6829419696157930e+00 err=0.0000000000000000e+00
n=11 int=1.6829419696157930e+00 err=0.0000000000000000e+00
n=12 int=1.6829419696157928e+00 err=2.2204460492503131e-16
n=13 int=1.6829419696157930e+00 err=0.0000000000000000e+00
n=14 int=1.6829419696157932e+00 err=-2.2204460492503131e-16
n=15 int=1.6829419696157935e+00 err=-4.4408920985006262e-16
n=16 int=1.6829419696157930e+00 err=0.0000000000000000e+00

end
```

# Future Work

It seems to me that the bottleneck for this method is the use of a bisection search on a uniform grid of intervals to find the zeroes. I strongly suspect that with a bit more knowledge about the behaviour of the Legendre polynomials, it should be possible to use a far superior algorithm to find the zeros of the Legendre polynomials. The method would be much faster subsequently.
