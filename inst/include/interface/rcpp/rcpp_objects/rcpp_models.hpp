#ifndef FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_MODELS_HPP
#define FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_MODELS_HPP

#include <set>
#include "../../../common/def.hpp"
#include "../../../models/fisheries_models.hpp"
#include "../../../utilities/fims_json.hpp"
#include "rcpp_population.hpp"

#include "rcpp_interface_base.hpp"

class FisheryModelInterfaceBase : public FIMSRcppInterfaceBase
{
public:
    /**
     * @brief The static id of the FleetInterfaceBase object.
     */
    static uint32_t id_g;
    /**
     * @brief The local id of the FleetInterfaceBase object.
     */
    uint32_t id;
    /**
     * @brief The map associating the IDs of FleetInterfaceBase to the objects.
     * This is a live object, which is an object that has been created and lives
     * in memory.
     */
    static std::map<uint32_t, FisheryModelInterfaceBase *> live_objects;

    /**
     * @brief The constructor.
     */
    FisheryModelInterfaceBase()
    {
        this->id = FisheryModelInterfaceBase::id_g++;
        /* Create instance of map: key is id and value is pointer to
        FleetInterfaceBase */
        FisheryModelInterfaceBase::live_objects[this->id] = this;
    }

    /**
     * @brief Construct a new Data Interface Base object
     *
     * @param other
     */
    FisheryModelInterfaceBase(const FisheryModelInterfaceBase &other) : id(other.id)
    {
    }

    /**
     * @brief The destructor.
     */
    virtual ~FisheryModelInterfaceBase()
    {
    }

    /**
     * @brief Get the ID for the child fleet interface objects to inherit.
     */
    virtual uint32_t get_id() = 0;
};
// static id of the FleetInterfaceBase object
uint32_t FisheryModelInterfaceBase::id_g = 1;
// local id of the FleetInterfaceBase object map relating the ID of the
// FleetInterfaceBase to the FleetInterfaceBase objects
std::map<uint32_t, FisheryModelInterfaceBase *> FisheryModelInterfaceBase::live_objects;

class CatchAtAgeInterface : public FisheryModelInterfaceBase
{
    std::set<uint32_t> population_ids;
    typedef typename std::set<uint32_t>::iterator population_id_iterator;

public:
    CatchAtAgeInterface() : FisheryModelInterfaceBase()
    {
    }

    void AddPopulation(uint32_t id)
    {
        this->population_ids.insert(id);
    }

    virtual uint32_t get_id()
    {
        return this->id;
    }

    void Show()
    {
        std::shared_ptr<fims_info::Information<double>> info =
            fims_info::Information<double>::GetInstance();

        fims_popdy::CatchAtAge<double> *model = (fims_popdy::CatchAtAge<double> *)info->models_map[this->get_id()].get();
        //        model->Show();
        std::cout << this->to_json(); // fims::JsonParser::PrettyFormatJSON(model->ToJSON());

        std::ofstream o("test.json");
        o << this->to_json();
        o.close();
    }

    virtual void finalize()
    {
    }

    std::string to_json()
    {

        typename std::map<uint32_t, std::shared_ptr<PopulationInterfaceBase>>::iterator pit;
        std::vector<uint32_t> pop_ids(this->population_ids.begin(), this->population_ids.end());

        std::stringstream ss;
        std::set<uint32_t> fleet_ids; // all fleets in the model
        ss << "{\n";
        ss << "\"model\" : \"catch_at_age\",\n";
        ss << "\"id\" : " << this->get_id() << ",\n";
        ss << "\"populations\" : [\n";
        // loop through populations for this model
        std::vector<std::string> pop_strings;
        for (size_t p = 0; p < pop_ids.size(); p++)
        {

            pit = PopulationInterfaceBase::live_objects.find(pop_ids[p]);
            if (pit != PopulationInterfaceBase::live_objects.end())
            {
                PopulationInterface *pop = (PopulationInterface *)(*pit).second.get();
                fleet_ids.insert(pop->fleet_ids->begin(), pop->fleet_ids->end());
                pop_strings.push_back((*pit).second->to_json());
                // if (p == pop_ids.size() - 1)
                // {
                //     ss << (*pit).second->to_json();
                // }
                // else
                // {
                //     ss << (*pit).second->to_json() << ",\n";
                // }
            }
        }
        if (pop_strings.size() > 0)
        {
            for (size_t i = 0; i < pop_strings.size() - 1; i++)
            {
                ss << pop_strings[i] << ",\n";
            }
            ss << pop_strings[pop_strings.size() - 1] << "\n";
        }
        ss << "],\n";
        ss << "\"fleets\" : [\n";
        typename std::map<uint32_t, std::shared_ptr<FleetInterfaceBase>>::iterator fit;
        // all fleets encapuslated in this model run
        std::vector<uint32_t> fids(fleet_ids.begin(), fleet_ids.end());
        // // loop through fleets for this model
        for (size_t f = 0; f < fids.size(); f++)
        {
            fit = FleetInterfaceBase::live_objects.find(fids[f]);
            if (fit != FleetInterfaceBase::live_objects.end())
            {
                if (f == fids.size() - 1)
                {
                    ss << (*fit).second->to_json();
                }
                else
                {
                    ss << (*fit).second->to_json() << ",\n";
                }
            }
        }

        ss << "]\n}";

        return fims::JsonParser::PrettyFormatJSON(ss.str());
    }

#ifdef TMB_MODEL

    template <typename Type>
    bool add_to_fims_tmb_internal()
    {
        std::shared_ptr<fims_info::Information<Type>> info =
            fims_info::Information<Type>::GetInstance();

        std::shared_ptr<fims_popdy::CatchAtAge<Type>> model = std::make_shared<fims_popdy::CatchAtAge<Type>>();

        population_id_iterator it;

        for (it = this->population_ids.begin(); it != this->population_ids.end(); ++it)
        {
            model->AddPopulation((*it));
        }

        // add to Information
        info->models_map[this->get_id()] = model;

        return true;
    }

    virtual bool add_to_fims_tmb()
    {
        FIMS_INFO_LOG("adding CAA model object to TMB");
        this->add_to_fims_tmb_internal<TMB_FIMS_REAL_TYPE>();
        this->add_to_fims_tmb_internal<TMB_FIMS_FIRST_ORDER>();
        this->add_to_fims_tmb_internal<TMB_FIMS_SECOND_ORDER>();
        this->add_to_fims_tmb_internal<TMB_FIMS_THIRD_ORDER>();

        return true;
    }

#endif
};

#endif