/*
 * File:   population.hpp
 *
 * This File is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project. See LICENSE in the
 * source folder for reuse information.
 *
 * Population module file
 * The purpose of this file is to define the Statistical Catch at Age Population class and its fields
 * and methods.
 *
 *
 */
#ifndef FIMS_POPULATION_DYNAMICS_SPM_HPP
#define FIMS_POPULATION_DYNAMICS_SPM_HPP

#include "../../common/model_object.hpp"
#include "../fleet/fleet.hpp"
#include "../../interface/interface.hpp"
#include "../population/functors/population_base.hpp"

namespace fims_popdy {

    template <typename Type>
struct SPM : public fims_model_object::FIMSObject<Type> {
  static uint32_t id_g; /*!< reference id for population object*/
  size_t nyears;        /*!< total number of years in the fishery*/
  size_t nfleets;       /*!< total number of fleets in the fishery*/


  fims::Vector<double> years;     /*!< vector of years for referencing*/
  fims::Vector<Type> r;           /*!< estimated parameter: intrinsic rate of population growth*/
  fims::Vector<Type> K;           /*!< estimated parameter: population carrying capacity*/
  fims::Vector<Type> m;           /*!< estimated parameter: production shape*/
  fims::Vector<Type> psi;         /*!< estimated parameter: initial depletion*/
  fims::Vector<Type> depletion_P;   /*!< Derived quantity: annual depletion*/
  fims::Vector<Type> expected_index; /*!< Expected values: index of abundance*/

  #ifdef TMB_MODEL
  ::objective_function<Type>
      *of;  // :: references global namespace, defined in src/FIMS.cpp,
            // available anywhere in the R package
    #endif


  SPM() { this->id = SPM::id_g++; }
  
  void Initialize(nyears){
    this->nyears = nyears;
    nfleets = fleets.size();
    expected_index.resize(nyears * nfleets);
    depletion_P.resize(nyears);
    expected_depletion_P.resize(nyears);
    years.resize(nyears);
    r.resize(1);
    K.resize(1);
    m.resize(1);
    psi.resize(1);
  }

  void Prepare() {
     for (size_t fleet = 0; fleet < this->fleets.size(); fleet++) {
      this->fleets[fleet]->Prepare();
    }
    std::fill(expected_index.begin(), expected_index.end(), 0.0);
    std::fill(depletion_P.begin(), depletion_P.end(), 0.0);
    std::fill(r.begin(), r.end(), 0.0);
    std::fill(K.begin(), K.end(), 0.0);
    std::fill(m.begin(), m.end(), 0.0);
    std::fill(psi.begin(), psi.end(), 0.0);
  }

  /**
 * @brief Calculates the yearly depletion
 */
void calculateExpectedDepletionP(year){
     // Get the previous year's depletion value
      expected_depletion_previous_year = this->expected_depletion_P[year-1];

      // This is the Pella-Tomlinson production model expression.
      this->expected_depletion_P = expected_depletion_previous_year +
            (r/(m-1.0))*expected_depletion_previous_year*(1.0-std::pow(expected_depletion_previous_year,this->m-1.0)) -
            obs_catch[year-1]/K;
}

};
}