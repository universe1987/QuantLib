/*! \file complexErf.hpp
    \brief complex error function
*/

#ifndef quantlib_complexErf_hpp
#define quantlib_complexErf_hpp

#include <complex>

namespace QuantLib {
/* Reference: Benhamou, Eric and Croissant, Olivier, Local Time for the SABR Model:
   Connection with the 'Complex' Black Scholes and Application to CMS and Spread Options (October 2007).
   Available at SSRN: https://ssrn.com/abstract=1064461 */
std::complex<double> erf(const std::complex<double>& z, const std::size_t order = 10);
} // namespace QuantLib

#endif
