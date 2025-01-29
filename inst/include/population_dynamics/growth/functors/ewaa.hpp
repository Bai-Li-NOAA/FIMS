/**
 * @file ewaa.hpp
 * @brief Declares the growth functor class which is the base class for all
 * growth functors.
 * @copyright This file is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project. See LICENSE in the source
 * folder for reuse information.
 */
#ifndef POPULATION_DYNAMICS_GROWTH_EWAA_HPP
#define POPULATION_DYNAMICS_GROWTH_EWAA_HPP

//#include "../../../interface/interface.hpp"
#include <map>

#include "growth_base.hpp"

namespace fims_popdy {

    /**
     *  @brief EWAAgrowth class that returns the EWAA function value.
     */
    template <typename Type>
    struct EWAAgrowth : public GrowthBase<Type> {
        // add submodule class members here
        // these include parameters of the submodule
        // a map looks up values based on a reference key
        // in this case, our key is age (first double), and
        //  the value is the weight at that age (second double)
        std::map<int, std::map<double, double> > ewaa; /**<map of doubles for EWAA values by age,
          where age starts at zero > */
        typedef typename std::map<int, std::map<double, double> >::iterator year_iterator; /**< Iterator for ewaa at year map object > */

        typedef typename std::map<double, double>::iterator weight_iterator; /**< Iterator for ewaa map object > */

        EWAAgrowth() : GrowthBase<Type>() {
        }

        virtual ~EWAAgrowth() {
        }

        /**
         * @brief Returns the weight at age a (in kg) from the input vector.
         *
         * @param a  age of the fish, the age vector must start at zero
         */
        virtual const Type evaluate(int year, const double& a) {
            year_iterator it = this->ewaa.find(year);
            if (it != this->ewaa.end()) {
                weight_iterator jt = (*it).second.find(a);
                if (jt == (*it).second.end()) {
                    return 0.0;
                }else{
                    return (*jt).second;
                }
            }
            return 0.0;
        }
    };
} // namespace fims_popdy
#endif /* POPULATION_DYNAMICS_GROWTH_EWAA_HPP */