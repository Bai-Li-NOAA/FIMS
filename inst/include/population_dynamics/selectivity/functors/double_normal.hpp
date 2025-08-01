/**
 * @file double_normal.hpp
 * @brief Declares the DoubleNormalSelectivity class which implements the
 * double normal selectivity function in Methot and Wetzel (2013).
 * @copyright This file is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project. See LICENSE in the source
 * folder for reuse information.
 */
#ifndef POPULATION_DYNAMICS_SELECTIVITY_DOUBLE_NORMAL_HPP
#define POPULATION_DYNAMICS_SELECTIVITY_DOUBLE_NORMAL_HPP

// #include "../../../interface/interface.hpp"
#include "../../../common/fims_math.hpp"
#include "../../../common/fims_vector.hpp"
#include "../../../common/information.hpp"//Added in case I can run GetNages() here
#include "selectivity_base.hpp"

namespace fims_popdy {

/**
 * @brief DoubleNormalSelectivity class that constructs and runs the double normal
 * function, defined below.
 */
template <typename Type>
struct DoubleNormalSelectivity : public SelectivityBase<Type> {
  size_t nages;         /*< the number of ages in the model*/
  fims::Vector<Type> age_peak_sel_start; /**< age at which selectivity=1
            starts, or p1 */
  fims::Vector<Type> width_peak_sel; /**< width of "top" in which selectivity=1, 
            or p2; determines the age at which selectivity=1 ends */
  fims::Vector<Type> slope_asc; /**< slope of the ascending section, or p3;
            ?: Is it problematic to use the same param name as dbl logistic? */
  fims::Vector<Type> slope_desc; /**< slope of the descending seciton, or p4;
            ?: Is it problematic to use the same param name as dbl logistic? */
  fims::Vector<Type> sel_age_zero_logit; /** selectivity at age0 (parameterized
            in logit space), or p5 */
  fims::Vector<Type> sel_age_A_logit; /** selectivity at age A 
            (parameterized in logit space), or p6 */

  DoubleNormalSelectivity() : SelectivityBase<Type>() {}

  virtual ~DoubleNormalSelectivity() {}

  /**
   * @brief Method of the double normal class that implements the
   * double normal function, provided below.
   *
   * UPDATE
   *
   * @param x  The independent variable in the double normal function (e.g.,
   * age or size in selectivity).
   */
  virtual const Type evaluate(const int nages,
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
      // Should use fims_math::inv_logit here instead, w/ a=0 and b=1
      // Am I using fims_math::pow() correctly?
    this->nages = nages; // assign value to nages (borrow from fleet.hpp syntax)
    const Type sel_age_zero = static_cast<Type>(1.0) / 
      (static_cast<Type>(1.0) + exp(Type(-1.0) * sel_age_zero_logit[0]));
    const Type sel_age_A = static_cast<Type>(1.0) / 
      (static_cast<Type>(1.0) + exp(Type(-1.0) * sel_age_A_logit[0]));
    const Type gamma = age_peak_sel_start[0] + static_cast<Type>(1.0) + 
      (Type(0.99) * nages - age_peak_sel_start[0] - static_cast<Type>(1.0)) / 
      (static_cast<Type>(1.0) + exp(Type(-1.0) * width_peak_sel[0]));
    const Type alpha_a = sel_age_zero +
      (static_cast<Type>(1.0) - sel_age_zero) *
      (exp(Type(-1.0) * pow<Type>((x - age_peak_sel_start[0]), static_cast<Type>(2.0)) / 
        exp(slope_asc[0])) - 
        exp(fims_math::pow<Type>(age_peak_sel_start[0], static_cast<Type>(2.0)) / 
        exp(slope_asc[0]))) / 
      (static_cast<Type>(1.0) - exp(pow<Type>(age_peak_sel_start[0], static_cast<Type>(2.0)) / 
        exp(slope_asc[0])));
    const Type beta_a = static_cast<Type>(1.0) + 
      (sel_age_A - static_cast<Type>(1.0)) * 
      (exp(Type(-1.0) * (pow((x - gamma), static_cast<Type>(2.0))) / 
        exp(slope_desc[0])) - static_cast<Type>(1.0)) / 
      (exp(Type(-1.0) * (pow((x - gamma), static_cast<Type>(2.0))) / 
        exp(slope_desc[0])) - static_cast<Type>(1.0));
    const Type j_one_a = pow((static_cast<Type>(1.0) + 
      exp(Type(-20.0) * (x - age_peak_sel_start[0]) / 
        (static_cast<Type>(1.0) + abs(x - age_peak_sel_start[0])))), 
      Type(-1.0));
    const Type j_two_a = pow((static_cast<Type>(1.0) + 
      exp(Type(-20.0) * (x - gamma) / 
        (static_cast<Type>(1.0) + abs(x - gamma)))), 
      Type(-1.0));
    return alpha_a * (static_cast<Type>(1.0) - j_one_a) +
      j_one_a * ((static_cast<Type>(1.0) - j_two_a) +
      j_two_a * beta_a);
  }

  /**
   * @brief Method of the double normal selectivity class that implements the
   * double normal function, provided below.
   *
   * UPDATE
   *
   * @param x  The independent variable in the double normal function (e.g.,
   * age or size in selectivity).
   * @param pos Position index, e.g., which year.
   */
  virtual const Type evaluate(const int nages,
                              const Type &age_peak_sel_start,
                              const Type &width_peak_sel,
                              const Type &slope_asc,
                              const Type &slope_desc,
                              const Type &sel_age_zero_logit,
                              const Type &sel_age_A_logit,
                              const Type &x,
                              size_t pos) {
    // Creating a bunch of placeholder variables for convenience;
    // Plan to remove and improve code efficiency later
    // ?s:
      // Am I using static_cast correctly? Trying to specify fixed value of 1, etc
      // Should use fims_math::inv_logit here instead, w/ a=0 and b=1
      // Am I using fims_math::pow() correctly?

    this->nages = nages; // assign value to nages (borrow from fleet.hpp syntax)
    const Type sel_age_zero = static_cast<Type>(1.0) / 
      (static_cast<Type>(1.0) + exp(Type(-1.0) * sel_age_zero_logit.get_force_scalar(pos)));
    const Type sel_age_A = static_cast<Type>(1.0) / 
      (static_cast<Type>(1.0) + exp(Type(-1.0) * sel_age_A_logit.get_force_scalar(pos)));
    const Type gamma = age_peak_sel_start.get_force_scalar(pos) + static_cast<Type>(1.0) + 
      (Type(0.99) * nages - age_peak_sel_start.get_force_scalar(pos) - static_cast<Type>(1.0)) / 
      (static_cast<Type>(1.0) + exp(Type(-1.0) * width_peak_sel.get_force_scalar(pos)));
    const Type alpha_a = sel_age_zero +
      (static_cast<Type>(1.0) - sel_age_zero) *
      (exp(Type(-1.0) * pow<Type>((x - age_peak_sel_start.get_force_scalar(pos)), static_cast<Type>(2.0)) / 
        exp(slope_asc.get_force_scalar(pos))) - 
        exp(fims_math::pow<Type>(age_peak_sel_start.get_force_scalar(pos), static_cast<Type>(2.0)) / 
        exp(slope_asc.get_force_scalar(pos)))) / 
      (static_cast<Type>(1.0) - exp(pow<Type>(age_peak_sel_start.get_force_scalar(pos), static_cast<Type>(2.0)) / 
        exp(slope_asc.get_force_scalar(pos))));
    const Type beta_a = static_cast<Type>(1.0) + 
      (sel_age_A - static_cast<Type>(1.0)) * 
      (exp(Type(-1.0) * (pow((x - gamma), static_cast<Type>(2.0))) / 
        exp(slope_desc.get_force_scalar(pos))) - static_cast<Type>(1.0)) / 
      (exp(Type(-1.0) * (pow((x - gamma), static_cast<Type>(2.0))) / 
        exp(slope_desc.get_force_scalar(pos))) - static_cast<Type>(1.0));
    const Type j_one_a = pow((static_cast<Type>(1.0) + 
      exp(Type(-20.0) * (x - age_peak_sel_start.get_force_scalar(pos)) / 
        (static_cast<Type>(1.0) + abs(x - age_peak_sel_start.get_force_scalar(pos))))), 
      Type(-1.0));
    const Type j_two_a = pow((static_cast<Type>(1.0) + 
      exp(Type(-20.0) * (x - gamma) / 
        (static_cast<Type>(1.0) + abs(x - gamma)))), 
      Type(-1.0));
    return alpha_a * (static_cast<Type>(1.0) - j_one_a) +
      j_one_a * ((static_cast<Type>(1.0) - j_two_a) +
      j_two_a * beta_a);
  }
};

}  // namespace fims_popdy

#endif /* POPULATION_DYNAMICS_SELECTIVITY_DOUBLE_NORMAL_HPP */
