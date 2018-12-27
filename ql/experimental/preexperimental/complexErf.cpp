#include "complexErf.hpp"

#include <ql/math/comparison.hpp>
#include <ql/math/distributions/normaldistribution.hpp>

#include <boost/math/special_functions/erf.hpp>

namespace QuantLib {
std::complex<double> erf(const std::complex<double>& z, const std::size_t order) {
    const double x = z.real();
    const double y = z.imag();
    double erfr = boost::math::erf(x);
    double erfi = 0.0;
    const double emxs = std::exp(-x * x);
    if (!close_enough(x, 0.0)) {
        erfr += emxs / (2.0 * M_PI * x) * (1.0 - cos(2.0 * x * y));
        erfi += emxs / (2.0 * M_PI * x) * sin(2.0 * x * y);
    } else {
        erfi += y / (M_PI);
    }
    double rr = 0.0, ri = 0.0;
    for (int n = 1; n <= order; ++n) {
        double nd = static_cast<double>(n);
        rr += exp(-0.25 * nd * nd) / (nd * nd + 4.0 * x * x) *
              (2.0 * x - 2.0 * x * cosh(nd * y) * cos(2.0 * x * y) + nd * sinh(n * y) * sin(2.0 * x * y));
        ri += exp(-0.25 * nd * nd) / (nd * nd + 4.0 * x * x) *
              (2.0 * x * cosh(nd * y) * sin(2.0 * x * y) + nd * sinh(nd * y) * cos(2.0 * x * y));
    }
    rr *= 2.0 / M_PI * emxs;
    ri *= 2.0 / M_PI * emxs;
    return std::complex<double>(erfr + rr, erfi + ri);
}
} // namespace QuantLib
