/*
 * File:   rcpp_selectivity.hpp
 *
 * Author: Matthew Supernaw
 * National Oceanic and Atmospheric Administration
 * National Marine Fisheries Service
 * Email: matthew.supernaw@noaa.gov
 *
 * Created on May 31, 2022 at 12:04 PM
 *
 * This File is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project.
 *
 * This software is a "United States Government Work" under the terms of the
 * United States Copyright Act.  It was written as part of the author's official
 * duties as a United States Government employee and thus cannot be copyrighted.
 * This software is freely available to the public for use. The National Oceanic
 * And Atmospheric Administration and the U.S. Government have not placed any
 * restriction on its use or reproduction.  Although all reasonable efforts have
 * been taken to ensure the accuracy and reliability of the software and data,
 * the National Oceanic And Atmospheric Administration and the U.S. Government
 * do not and cannot warrant the performance or results that may be obtained by
 * using this  software or data. The National Oceanic And Atmospheric
 * Administration and the U.S. Government disclaim all warranties, express or
 * implied, including warranties of performance, merchantability or fitness
 * for any particular purpose.
 *
 * Please cite the author(s) in any work or product based on this material.
 *
 */
#ifndef FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_SELECTIVITY_HPP
#define FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_SELECTIVITY_HPP

#include "../../../population_dynamics/selectivity/selectivity.hpp"
#include "rcpp_interface_base.hpp"

/****************************************************************
 * Selectivity Rcpp interface                                   *
 ***************************************************************/
/**
 * @brief SelectivityInterfaceBase class should be inherited to
 * define different Rcpp interfaces for each possible Selectivity function
 * */
class SelectivityInterfaceBase : public FIMSRcppInterfaceBase {
 public:
  static uint32_t id_g;
  uint32_t id;
  static std::map<uint32_t, SelectivityInterfaceBase*> selectivity_objects;

  SelectivityInterfaceBase() {
    this->id = SelectivityInterfaceBase::id_g++;
    SelectivityInterfaceBase::selectivity_objects[this->id] = this;
    FIMSRcppInterfaceBase::fims_interface_objects.push_back(this);
  }

  virtual uint32_t get_id() = 0;
};

uint32_t SelectivityInterfaceBase::id_g = 1;
std::map<uint32_t, SelectivityInterfaceBase*>
    SelectivityInterfaceBase::selectivity_objects;

/**
 * @brief Rcpp interface for logistic selectivity as an S4 object. To
 * instantiate from R: logistic_selectivity <- new(fims$logistic_selectivity)
 */
class LogisticSelectivityInterface : public SelectivityInterfaceBase {
 public:
  Parameter median; /**< the index value at which the response reaches .5 */
  Parameter slope;  /**< the width of the curve at the median */

  LogisticSelectivityInterface() : SelectivityInterfaceBase() {}

  virtual ~LogisticSelectivityInterface() {}

  /** @brief returns the id for the logistic selectivity interface */
  virtual uint32_t get_id() { return this->id; }

  /** @brief this adds the parameter values and derivatives to the TMB model
   * object */
   template<typename Type> 
  virtual bool add_to_fims_tmb() {
    std::shared_ptr<fims::Information<Type> > d0 =
        fims::Information<Type>::GetInstance();

    std::shared_ptr<fims::LogisticSelectivity<Type> > ls0 =
        std::make_shared<fims::LogisticSelectivity<Type> >();

    // set relative info
    ls0->id = this->id;
    ls0->median = this->median.value;
    if (this->median.estimated) {
      if (this->median.is_random_effect) {
        d0->RegisterRandomEffect(ls0->median);
      } else {
        d0->RegisterParameter(ls0->median);
      }
    }
    ls0->slope = this->slope.value;
    if (this->slope.estimated) {
      if (this->slope.is_random_effect) {
        d0->RegisterRandomEffect(ls0->slope);
      } else {
        d0->RegisterParameter(ls0->slope);
      }
    }

    // add to Information
    d0->selectivity_models[ls0->id] = ls0;

   
    return true;
  }
};

#endif