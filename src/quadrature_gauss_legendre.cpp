#include "quadrature_gauss_legendre.hpp"

void Quadrature_gauss_legendre::populate(int order)
{
  try
  {
    points = quad.at(order-1);
  }
  catch (...)
  {
    std::cerr << "Unacceptable Gauss-Legendre order: " << std::to_string(order) << std::endl;
    std::abort();
  }
}

const std::vector<std::vector<Quadrature_point>> Quadrature_gauss_legendre::quad = std::vector<std::vector<Quadrature_point>>{
#include "../gauss_legendre/gauss_legendre.txt"
};

double Quadrature_gauss_legendre::integrate(const std::function<double(double)> & f, const double xlo, const double xhi) const
{
  double xsum{0.0};
  for (const auto & qp : points)
  {
    const double xi{0.5*(xhi-xlo)*qp.x + 0.5*(xhi+xlo)}; // coordinate transform
    xsum += qp.w * f(xi);
  }
  return (xhi-xlo)*0.5*xsum;
}
