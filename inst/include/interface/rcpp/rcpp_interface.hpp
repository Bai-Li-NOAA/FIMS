/*
 * File:   rcpp_interface.hpp
 *
 *
 * This File is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project. See LICENSE file for reuse
 * information.
 *
 *
 */
#ifndef FIMS_INTERFACE_RCPP_INTERFACE_HPP
#define FIMS_INTERFACE_RCPP_INTERFACE_HPP

#include "rcpp_objects/rcpp_fishing_mortality.hpp"
#include "rcpp_objects/rcpp_fleet.hpp"
#include "rcpp_objects/rcpp_growth.hpp"
#include "rcpp_objects/rcpp_maturity.hpp"
#include "rcpp_objects/rcpp_natural_mortality.hpp"
#include "rcpp_objects/rcpp_population.hpp"
#include "rcpp_objects/rcpp_recruitment.hpp"
#include "rcpp_objects/rcpp_selectivity.hpp"
#include "rcpp_objects/rcpp_tmb_dnorm_distribution.hpp"

/**
 * @brief Create the TMB model object and add interface objects to it.
 */
template<typename Type>
bool CreateTMBModel() {
  for (size_t i = 0; i < FIMSRcppInterfaceBase::fims_interface_objects.size();
       i++) {
    FIMSRcppInterfaceBase::fims_interface_objects[i]->add_to_fims_tmb<Type>();
  }

  // base model
  std::shared_ptr<fims::Information<Type> > d0 =
      fims::Information<Type>::GetInstance();
  d0->CreateModel();

  return true;
}

RCPP_EXPOSED_CLASS(Parameter)
RCPP_MODULE(fims) {
  Rcpp::function("CreateTMBModel", &CreateTMBModel);

  Rcpp::class_<Parameter>("Parameter")
      .constructor()
      .constructor<double>()
      .constructor<Parameter>()
      .field("value", &Parameter::value)
      .field("min", &Parameter::min)
      .field("max", &Parameter::max)
      .field("is_random_effect", &Parameter::is_random_effect)
      .field("estimated", &Parameter::estimated);

  Rcpp::class_<BevertonHoltRecruitmentInterface>("BevertonHoltRecruitment")
      .constructor()
      .field("steep", &BevertonHoltRecruitmentInterface::steep)
      .field("rzero", &BevertonHoltRecruitmentInterface::rzero)
      .field("phizero", &BevertonHoltRecruitmentInterface::phizero)
      .method("get_id", &BevertonHoltRecruitmentInterface::get_id);

  Rcpp::class_<LogisticSelectivityInterface>("LogisticSelectivity")
      .constructor()
      .field("median", &LogisticSelectivityInterface::median)
      .field("slope", &LogisticSelectivityInterface::slope)
      .method("get_id", &LogisticSelectivityInterface::get_id);

  Rcpp::class_<DnormDistributionsInterface>("TMBDnormDistribution")
      .constructor()
      .method("get_id", &DnormDistributionsInterface::get_id)
      .method("evaluate", &DnormDistributionsInterface::evaluate)
      .field("x", &DnormDistributionsInterface::x)
      .field("mean", &DnormDistributionsInterface::mean)
      .field("sd", &DnormDistributionsInterface::sd);

  Rcpp::class_<EWAAGrowthInterface>("EWAAgrowth")
      .constructor()
      .field("ages", &EWAAGrowthInterface::ages)
      .field("weights", &EWAAGrowthInterface::weights)
      .method("evaluate", &EWAAGrowthInterface::evaluate);
}

#endif /* RCPP_INTERFACE_HPP */
