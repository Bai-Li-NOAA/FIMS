/**
 * @file fims_math.hpp
 * @brief TODO: provide a brief description.
 * @copyright This file is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project. See LICENSE in the source
 * folder for reuse information.
 */
#ifndef FIMS_MATH_HPP
#define FIMS_MATH_HPP

// note: this is modeling platform specific, must be controlled by
// preprocessing macros
#include <cmath>
#include <random>
#include <sstream>

#include "../interface/interface.hpp"
#include "fims_vector.hpp"

namespace fims_math {
#ifdef STD_LIB

/**
 * @brief The exponential function.
 *
 * @param x value to exponentiate. Please use fims_math::exp<double>(x) if x is
 * an integer.
 * @return the exponentiated value
 */
template <class Type>
inline const Type exp(const Type &x) {
  return std::exp(x);
}

/**
 * @brief The natural log function (base e)
 * @param x the value to take the log of. Please use fims_math::log<double>(x)
 * if x is an integer.
 * @return
 */
template <class Type>
inline const Type log(const Type &x) {
  return std::log(x);
}

template <class Type>
inline const Type cos(const Type &x) {
  return std::cos(x);
}

template <class Type>
inline const Type sqrt(const Type &x) {
  return std::sqrt(x);
}

template <class Type>
inline const Type pow(const Type &x, const Type &y) {
  return std::pow(x, y);
}

template <class Type>
inline const Type lgamma(const Type &x) {
  return std::lgamma(x);
}
#endif

#ifdef TMB_MODEL

/**
 * @brief The exponential function.
 * The code cannot be tested using the compilation flag
 * -DTMB_MODEL through CMake and Google Test
 * @param x value to exponentiate. Please use fims_math::exp<double>(x) if x is
 * an integer.
 * @return the exponentiated value
 */
template <class Type>
inline const Type exp(const Type &x) {
  using ::exp;
  return exp(x);
}

template <>
inline const double exp(const double &x) {
  return std::exp(x);
}

/**
 * @brief The natural log function (base e)
 * The code cannot be tested using the compilation flag
 * -DTMB_MODEL through CMake and Google Test.
 * @param x the value to take the log of. Please use fims_math::log<double>(x)
 * if x is an integer.
 * @return The natural log of the value.
 */
template <class Type>
inline const Type log(const Type &x) {
  return log(x);
}

template <>
inline const double log(const double &x) {
  return std::log(x);
}

template <class Type>
inline const Type cos(const Type &x) {
  return cos(x);
}

template <>
inline const double cos(const double &x) {
  return std::cos(x);
}

template <class Type>
inline const Type sqrt(const Type &x) {
  return sqrt(x);
}

template <>
inline const double sqrt(const double &x) {
  return std::sqrt(x);
}

template <class Type>
inline const Type pow(const Type &x, const Type &y) {
  return pow(x, y);
}

template <>
inline const double pow(const double &x, const double &y) {
  return std::pow(x, y);
}

template <class Type>
inline const Type lgamma(const Type &x) {
  using ::lgamma;
  return lgamma(x);
}

template <>
inline const double lgamma(const double &x) {
  return std::lgamma(x);
}

#endif

/**
 * @brief The general logistic function
 *
 * \f$ \frac{1.0}{ 1.0 + exp(-1.0 * slope (x - inflection_point))} \f$
 *
 * @param inflection_point the inflection point of the logistic function
 * @param slope the slope of the logistic function
 * @param x the index the logistic function should be evaluated at
 * @return
 */
template <class Type>
inline const Type logistic(const Type &inflection_point, const Type &slope,
                           const Type &x) {
  return static_cast<Type>(1.0) /
         (static_cast<Type>(1.0) +
          exp(Type(-1.0) * slope * (x - inflection_point)));
}

/**
 * @brief A logit function for bounding of parameters
 *
 * \f$ -\mathrm{log}(b-x) + \mathrm{log}(x-a) \f$
 * @param a lower bound
 * @param b upper bound
 * @param x the parameter in bounded space
 * @return the parameter in real space
 *
 */
template <class Type>
inline const Type logit(const Type &a, const Type &b, const Type &x) {
  return -fims_math::log(b - x) + fims_math::log(x - a);
}

/**
 * @brief An inverse logit function for bounding of parameters
 *
 * \f$ a+\frac{b-a}{1+\mathrm{exp}(-\mathrm{logit}(x))}\f$
 * @param a lower bound
 * @param b upper bound
 * @param logit_x the parameter in real space
 * @return the parameter in bounded space
 *
 */
template <class Type>
inline const Type inv_logit(const Type &a, const Type &b, const Type &logit_x) {
  return a + (b - a) / (static_cast<Type>(1.0) + fims_math::exp(-logit_x));
}

/**
 * @brief The general double logistic function
 *
 * \f$ \frac{1.0}{ 1.0 + exp(-1.0 * slope_{asc} (x - inflection_point_{asc}))}
 * \left(1-\frac{1.0}{ 1.0 + exp(-1.0 * slope_{desc} (x -
 * inflection_point_{desc}))} \right)\f$
 *
 * @param inflection_point_asc the inflection point of the ascending limb of the
 * double logistic function
 * @param slope_asc the slope of the ascending limb of the double logistic
 * function
 * @param inflection_point_desc the inflection point of the descending limb of
 * the double logistic function, where inflection_point_desc >
 * inflection_point_asc
 * @param slope_desc the slope of the descending limb of the double logistic
 * function
 * @param x the index the logistic function should be evaluated at
 * @return
 */

template <class Type>
inline const Type double_logistic(const Type &inflection_point_asc,
                                  const Type &slope_asc,
                                  const Type &inflection_point_desc,
                                  const Type &slope_desc, const Type &x) {
  return (static_cast<Type>(1.0)) /
         (static_cast<Type>(1.0) +
          exp(Type(-1.0) * slope_asc * (x - inflection_point_asc))) *
         (static_cast<Type>(1.0) -
          (static_cast<Type>(1.0)) /
              (static_cast<Type>(1.0) +
               exp(Type(-1.0) * slope_desc * (x - inflection_point_desc))));
}

/**
 *
 * Used when x could evaluate to zero, which will result in a NaN for
 * derivative values.
 *
 * Evaluates:
 *
 * \f$ (x^2+C)^.5 \f$
 *
 * @param x value to keep positive
 * @param C default = 1e-5
 * @return
 */
template <class Type>
const Type ad_fabs(const Type &x, Type C = 1e-5) {
  return sqrt((x * x) + C);
}

/**
 *
 * Returns the minimum between a and b in a continuous manner using:
 *
 * (a + b - fims_math::ad_fabs(a - b))*.5;
 * Reference: \ref fims_math::ad_fabs()
 *
 * This is an approximation with minimal error.
 *
 * @param a
 * @param b
 * @param C default = 1e-5
 * @return
 */

template <typename Type>
inline const Type ad_min(const Type &a, const Type &b, Type C = 1e-5) {
  return (a + b - fims_math::ad_fabs(a - b, C)) * static_cast<Type>(0.5);
}

/**
 * Returns the maximum between a and b in a continuous manner using:
 *
 * (a + b + fims_math::ad_fabs(a - b)) *.5;
 * Reference: \ref fims_math::ad_fabs()
 * This is an approximation with minimal error.
 *
 * @param a
 * @param b
 * @param C default = 1e-5
 * @return
 */
template <typename Type>
inline const Type ad_max(const Type &a, const Type &b, Type C = 1e-5) {
  return (a + b + fims_math::ad_fabs(a - b, C)) * static_cast<Type>(.5);
}

/**
 * Sum elements of a vector
 *
 * @brief
 *
 * @param v A vector of constants.
 * @return A single numeric value.
 */
template <class T>
T sum(const std::vector<T> &v) {
  T ret = 0.0;
  for (int i = 0; i < v.size(); i++) {
    ret += v[i];
  }
  return ret;
}

/**
 * Sum elements of a vector
 *
 * @brief
 *
 * @param v A vector of constants.
 * @return A single numeric value.
 */
template <class T>
T sum(const fims::Vector<T> &v) {
  T ret = 0.0;
  for (int i = 0; i < v.size(); i++) {
    ret += v[i];
  }
  return ret;
}

/**
 * @brief The double normal function for selectivity
 *
 *
 * @param age_peak_sel_start age at which selectivity=1 starts, or p1
 * @param width_peak_sel width of "top" in which selectivity=1, or p2; 
 * determines the age at which selectivity=1 ends
 * @param slope_asc slope of the ascending section, or p3
 * @param slope_desc slope of the descending seciton, or p4
 * @param sel_age_zero_logit selectivity at age-0 (parameterized in logit 
 * space), or p5
 * @param sel_age_A_logit selectivity at age A (parameterized in logit space), 
 * or p6
 * @param x the index the logistic function should be evaluated at
 * @return
 */

template <class Type>
inline const Type double_normal(//const Type nages, //Option B
                              const Type &age_peak_sel_start,
                              const Type &width_peak_sel,
                              const Type &slope_asc,
                              const Type &slope_desc,
                              const Type &sel_age_zero_logit,
                              const Type &sel_age_A_logit,
                              const Type &x) {
  // Creating a bunch of placeholder variables for convenience;
    // Plan to remove and improve code efficiency later
  // ?s:
    // Am I using static_cast correctly? Trying to specify fixed value of 1, etc
      // Do I need static_cast for Type<2.0>?
    // Should use fims_math::inv_logit here instead, w/ a=0 and b=1
    // Am I using fims_math::pow() correctly?
  const Type max_age = Type(12.0); // Option A
  //const Type max_age = nages // - static_cast<Type>(1.0); // Option B 
  const Type sel_age_zero = static_cast<Type>(1.0) / 
    (static_cast<Type>(1.0) + exp(Type(-1.0) * sel_age_zero_logit));
  const Type sel_age_A = static_cast<Type>(1.0) / 
    (static_cast<Type>(1.0) + exp(Type(-1.0) * sel_age_A_logit));
  const Type gamma = age_peak_sel_start + static_cast<Type>(1.0) + 
    (Type(0.99) * max_age - age_peak_sel_start - static_cast<Type>(1.0)) / 
    (static_cast<Type>(1.0) + exp(Type(-1.0) * width_peak_sel));
  const Type alpha_a = sel_age_zero +
    (static_cast<Type>(1.0) - sel_age_zero) *
    (exp(Type(-1.0) * pow<Type>((x - age_peak_sel_start), Type(2.0)) / 
      exp(slope_asc)) - 
      exp(Type(-1.0) * fims_math::pow<Type>(age_peak_sel_start, Type(2.0)) / 
      exp(slope_asc))) / 
    (static_cast<Type>(1.0) - exp(Type(-1.0) * 
      pow<Type>(age_peak_sel_start, Type(2.0)) / 
      exp(slope_asc)));
  const Type beta_a = static_cast<Type>(1.0) + 
    (sel_age_A - static_cast<Type>(1.0)) * 
    (exp(Type(-1.0) * (pow((x - gamma), Type(2.0))) / 
      exp(slope_desc)) - static_cast<Type>(1.0)) / 
    (exp(Type(-1.0) * (pow((max_age - gamma), Type(2.0))) / 
      exp(slope_desc)) - static_cast<Type>(1.0));
  const Type j_one_a = pow((static_cast<Type>(1.0) + 
    exp(Type(-20.0) * (x - age_peak_sel_start) / 
      (static_cast<Type>(1.0) + fims_math::ad_fabs(x - age_peak_sel_start)))), 
    Type(-1.0));
  const Type j_two_a = pow((static_cast<Type>(1.0) + 
    exp(Type(-20.0) * (x - gamma) / 
      (static_cast<Type>(1.0) + fims_math::ad_fabs(x - gamma)))), 
    Type(-1.0));
  return alpha_a * (static_cast<Type>(1.0) - j_one_a) +
    j_one_a * ((static_cast<Type>(1.0) - j_two_a) +
    j_two_a * beta_a);
}

}  // namespace fims_math

#endif /* FIMS_MATH_HPP */
