#include <boost/math/special_functions/erf.hpp>

extern "C" {
  double cpp_erf(double p) {
    return boost::math::erf(p);
  }

  double cpp_erf_inv(double p) {
    return boost::math::erf_inv(p);
  }
}

