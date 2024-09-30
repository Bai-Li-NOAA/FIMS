/*
 *
 * This File is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project. See LICENSE in the
 * source folder for reuse information.
 *
 * PopulationBase  file
 * The purpose of this file is to declare the PopulationBase class
 * which is the base class for all population functors.
 *
 * DEFINE guards for population module outline to define the
 * population hpp file if not already defined.
 */
#ifndef POPULATION_DYNAMICS_POPULATION_BASE_HPP
#define POPULATION_DYNAMICS_POPULATION_BASE_HPP

#include "../../../common/model_object.hpp"

namespace fims_popdy {

/** @brief Base class for all population functors.
 *
 * @tparam Type The type of the population functor.
 */

template <typename Type>
struct PopulationBase : public fims_model_object::FIMSObject<Type> {
  // id_g is the ID of the instance of the PopulationBase class.
  // this is like a memory tracker.
  // Assigning each one its own ID is a way to keep track of
  // all the instances of the PopulationBase class.
  static uint32_t
      id_g; /**< The ID of the instance of the PopulationBase class */

  /** @brief Constructor.
   */
  PopulationBase() {
    // increment id of the singleton population class
    this->id = PopulationBase::id_g++;
  }

  virtual ~PopulationBase() {}

  /**
   * @brief Calculates the population.
   */
  virtual const Type Evaluate() = 0;
};

// default id of the singleton population class
template <typename Type>
uint32_t PopulationBase<Type>::id_g = 0;

}  // namespace fims_popdy

#endif /* POPULATION_DYNAMICS_POPULATION_BASE_HPP */
