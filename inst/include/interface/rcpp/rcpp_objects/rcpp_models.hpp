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
    static std::map<uint32_t, std::shared_ptr<FisheryModelInterfaceBase>> live_objects;

    /**
     * @brief The constructor.
     */
    FisheryModelInterfaceBase()
    {
        this->id = FisheryModelInterfaceBase::id_g++;
        /* Create instance of map: key is id and value is pointer to
        FleetInterfaceBase */
        // FisheryModelInterfaceBase::live_objects[this->id] = this;
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

    virtual std::string to_json()
    {
        return "std::string to_json() not yet implemented.";
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
std::map<uint32_t, std::shared_ptr<FisheryModelInterfaceBase>> FisheryModelInterfaceBase::live_objects;

class CatchAtAgeInterface : public FisheryModelInterfaceBase
{
    std::shared_ptr<std::set<uint32_t>> population_ids;
    typedef typename std::set<uint32_t>::iterator population_id_iterator;

public:
    CatchAtAgeInterface() : FisheryModelInterfaceBase()
    {
        this->population_ids = std::make_shared<std::set<uint32_t>>();
        std::shared_ptr<CatchAtAgeInterface> caa = std::make_shared<CatchAtAgeInterface>(*this);
        FIMSRcppInterfaceBase::fims_interface_objects.push_back(caa);
        FisheryModelInterfaceBase::live_objects[this->id] = caa;
    }

    /**
     * @brief Construct a new Catch At Age Interface object
     *
     * @param other
     */
    CatchAtAgeInterface(const CatchAtAgeInterface &other) : FisheryModelInterfaceBase(other),
                                                            population_ids(other.population_ids)
    {
    }

    void AddPopulation(uint32_t id)
    {
        this->population_ids->insert(id);
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
        model->Show();
        // std::cout << this->to_json(); // fims::JsonParser::PrettyFormatJSON(model->ToJSON());

        // std::ofstream o("test.json");
        // o << this->to_json();
        // o.close();
    }

    virtual void finalize()
    {
    }

    std::string population_to_json(PopulationInterface *population_interface)
    {

        std::stringstream ss;

        std::shared_ptr<fims_info::Information<double>> info =
            fims_info::Information<double>::GetInstance();

        typename fims_info::Information<double>::population_iterator pit;

        pit = info->populations.find(population_interface->get_id());

        if (pit != info->populations.end())
        {
            std::shared_ptr<fims_popdy::Population<double>> &pop = (*pit).second;

            fims::Vector<double> &derived_ssb = pop->derived_quantities["spawning_biomass"];
            fims::Vector<double> &derived_naa = pop->derived_quantities["numbers_at_age"];
            fims::Vector<double> &derived_biomass = pop->derived_quantities["biomass"];
            fims::Vector<double> &derived_recruitment = pop->derived_quantities["recruitment"];

            // ToDo: add list of fleet ids operating on this population
            ss << "{\n";
            ss << " \"name\" : \"Population\",\n";

            ss << " \"type\" : \"population\",\n";
            ss << " \"tag\" : \"" << population_interface->name << "\",\n";
            ss << " \"id\": " << population_interface->id << ",\n";
            ss << " \"recruitment_id\": " << population_interface->recruitment_id << ",\n";
            ss << " \"growth_id\": " << population_interface->growth_id << ",\n";
            ss << " \"maturity_id\": " << population_interface->maturity_id << ",\n";

            ss << " \"parameters\": [\n{\n";
            ss << " \"name\": \"log_M\",\n";
            ss << " \"id\":" << population_interface->log_M.id_m << ",\n";
            ss << " \"type\": \"vector\",\n";
            ss << " \"values\": " << population_interface->log_M << "\n},\n";

            ss << "{\n";
            ss << "  \"name\": \"log_init_naa\",\n";
            ss << "  \"id\":" << population_interface->log_init_naa.id_m << ",\n";
            ss << "  \"type\": \"vector\",\n";
            ss << "  \"values\":" << population_interface->log_init_naa << " \n}],\n";

            ss << " \"derived_quantities\": [{\n";
            ss << "  \"name\": \"ssb\",\n";
            ss << " \"dimensions\" : [" << this->make_dimensions(1, population_interface->nyears + 1) << "],";
            ss << "  \"values\":[";
            if (derived_ssb.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_ssb.size() - 1; i++)
                {
                    ss << derived_ssb[i] << ", ";
                }
                ss << derived_ssb[derived_ssb.size() - 1] << "]\n";
            }
            ss << " },\n";

            ss << "{\n";
            ss << "   \"name\": \"naa\",\n";
            ss << " \"dimensions\" : [[" << this->make_dimensions(1, population_interface->nyears) << "],[" << this->make_dimensions(population_interface->ages[0], population_interface->ages[population_interface->ages.size() - 1], population_interface->nyears + 1) << "]],";
            ss << "   \"values\":[";
            if (derived_naa.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_naa.size() - 1; i++)
                {
                    ss << derived_naa[i] << ", ";
                }
                ss << derived_naa[derived_naa.size() - 1] << "]\n";
            }
            ss << " },\n";

            ss << "{\n";
            ss << "   \"name\": \"biomass\",\n";
            ss << " \"dimensions\" : [" << population_interface->make_dimensions(1, population_interface->nyears + 1) << "],";
            ss << "   \"values\":[";
            if (derived_biomass.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_biomass.size() - 1; i++)
                {
                    ss << derived_biomass[i] << ", ";
                }
                ss << derived_biomass[derived_biomass.size() - 1] << "]\n";
            }
            ss << " },\n";

            ss << "{\n";
            ss << "   \"name\": \"recruitment\",\n";
            ss << " \"dimensions\" : [" << this->make_dimensions(1, population_interface->nyears + 1) << "],";
            ss << "   \"values\":[";
            if (derived_recruitment.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_recruitment.size() - 1; i++)
                {
                    ss << derived_recruitment[i] << ", ";
                }
                ss << derived_recruitment[derived_recruitment.size() - 1] << "]\n";
            }
            ss << " }\n]\n";

            ss << "}";
        }
        else
        {
            ss << "{\n";
            ss << " \"name\" : \"Population\",\n";

            ss << " \"type\" : \"population\",\n";
            ss << " \"tag\" : \"" << population_interface->get_id() << " not found in Information.\",\n}";

#warning Add error log here
        }
        return ss.str();
    }

    std::string fleets_to_json(FleetInterface *fleet_interface)
    {
        std::stringstream ss;

        std::shared_ptr<fims_info::Information<double>> info =
            fims_info::Information<double>::GetInstance();

        typename fims_info::Information<double>::fleet_iterator fit;

        fit = info->fleets.find(fleet_interface->get_id());

        if (fit != info->fleets.end())
        {

            std::shared_ptr<fims_popdy::Fleet<double>> &fleet = (*fit).second;

            fims::Vector<double> &derived_caa = fleet->derived_quantities["catch_numbers_at_age"];
            fims::Vector<double> &derived_cal = fleet->derived_quantities["catch_numbers_at_length"];
            fims::Vector<double> &derived_proportion_cnaa = fleet->derived_quantities["proportion_catch_numbers_at_age"];
            fims::Vector<double> &derived_proportion_cnal = fleet->derived_quantities["proportion_catch_numbers_at_length"];
            fims::Vector<double> &derived_index = fleet->derived_quantities["expected_index"];
            fims::Vector<double> &derived_catch = fleet->derived_quantities["expected_catch"];
            fims::Vector<double> &derived_cwaa = fleet->derived_quantities["catch_weight_at_age"];
            fims::Vector<double> &derived_age_comp = fleet->derived_quantities["age_composition"];
            fims::Vector<double> &derived_length_comp = fleet->derived_quantities["length_composition"];

            ss << "{\n";
            ss << " \"name\" : \"Fleet\",\n";

            ss << " \"type\" : \"fleet\",\n";
            ss << " \"tag\" : \"" << fleet_interface->name << "\",\n";
            ss << " \"id\": " << fleet_interface->id << ",\n";
            ss << " \"is_survey\": " << fleet_interface->is_survey << ",\n";
            ss << " \"nlengths\": " << fleet_interface->nlengths.get() << ",\n";
            ss << "\"parameters\": [\n";
            ss << "{\n";
            ss << " \"name\": \"log_Fmort\",\n";
            ss << " \"id\":" << fleet_interface->log_Fmort.id_m << ",\n";
            ss << " \"type\": \"vector\",\n";
            ss << " \"values\": " << fleet_interface->log_Fmort << "\n},\n";
        
            ss << " {\n";
            ss << " \"name\": \"log_q\",\n";
            ss << " \"id\":" << fleet_interface->log_q.id_m << ",\n";
            ss << " \"type\": \"vector\",\n";
            ss << " \"values\": " << fleet_interface->log_q << "\n},\n";
            if (fleet_interface->nlengths > 0) {
              ss << " {\n";
              ss << " \"name\": \"age_length_conversion_matrix\",\n";
              ss << " \"id\":" << fleet_interface->age_length_conversion_matrix.id_m << ",\n";
              ss << " \"type\": \"vector\",\n";
              ss << " \"values\": " << fleet_interface->age_length_conversion_matrix << "\n}\n";
            }
            ss << " \"derived_quantities\": [{\n";
            ss << "  \"name\": \"cnaa\",\n";
            ss << " \"dimensions\" : [" << this->make_dimensions(1, fleet->nyears + 1) << "],";
            ss << "  \"values\":[";
            if(derived_caa.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_caa.size() - 1; i++)
                {
                    ss << derived_caa[i] << ", ";
                }
                ss << derived_caa[derived_caa.size() - 1] << "]\n";
            }
            ss << " },\n";
            ss << " {\n";   
            ss << "  \"name\": \"cnal\",\n";
            ss << "  \"values\":[";
            if (derived_cal.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_cal.size() - 1; i++)
                {
                    ss << derived_cal[i] << ", ";
                }
                ss << derived_cal[derived_cal.size() - 1] << "]\n";
            }
            ss << " },\n";
            ss << " {\n";
            ss << "  \"name\": \"cwaa\",\n";
            ss << "  \"values\":[";
            if (derived_cwaa.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_cwaa.size() - 1; i++)
                {
                    ss << derived_cwaa[i] << ", ";
                }
                ss << derived_cwaa[derived_cwaa.size() - 1] << "]\n";
            }
            ss << " },\n";
            ss << " {\n";
            ss << "  \"name\": \"proportion_catch_numbers_at_age\",\n";
            ss << "  \"values\":[";
            if (derived_proportion_cnaa.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_proportion_cnaa.size() - 1; i++)
                {
                    ss << derived_proportion_cnaa[i] << ", ";
                }
                ss << derived_proportion_cnaa[derived_proportion_cnaa.size() - 1] << "]\n";
            }
            ss << " },\n";
            ss << " {\n";
            ss << "  \"name\": \"proportion_catch_numbers_at_length\",\n";
            ss << "  \"values\":[";
            if (derived_proportion_cnal.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_proportion_cnal.size() - 1; i++)
                {
                    ss << derived_proportion_cnal[i] << ", ";
                }
                ss << derived_proportion_cnal[derived_proportion_cnal.size() - 1] << "]\n";
            }
            ss << " },\n";
            ss << " {\n";
            ss << "  \"name\": \"expected_index\",\n";
            ss << "  \"values\":[";
            if (derived_index.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_index.size() - 1; i++)
                {
                    ss << derived_index[i] << ", ";
                }
                ss << derived_index[derived_index.size() - 1] << "]\n";
            }
            ss << " },\n";
            ss << " {\n";
            ss << "  \"name\": \"expected_catch\",\n";
            ss << "  \"values\":[";
            if (derived_catch.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_catch.size() - 1; i++)
                {
                    ss << derived_catch[i] << ", ";
                }
                ss << derived_catch[derived_catch.size() - 1] << "]\n";
            }
            ss << " },\n";
            ss << " {\n";
            ss << "  \"name\": \"age_composition \",\n";
            ss << "  \"values\":[";
            if (derived_age_comp.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_age_comp.size() - 1; i++)
                {
                    ss << derived_age_comp[i] << ", ";
                }
                ss << derived_age_comp[derived_age_comp.size() - 1] << "]\n";
            }
            ss << " },\n";
            ss << " {\n";
            ss << "  \"name\": \"length_composition \",\n";
            ss << "  \"values\":[";
            if (derived_length_comp.size() == 0)
            {
                ss << "]\n";
            }
            else
            {
                for (size_t i = 0; i < derived_length_comp.size() - 1; i++)
                {
                    ss << derived_length_comp[i] << ", ";
                }
                ss << derived_length_comp[derived_length_comp.size() - 1] << "]\n";
            }
            ss << " }\n]\n";
            ss << "}";
            
        }
        else
        {
            ss << "{\n";
            ss << " \"name\" : \"Fleet\",\n";

            ss << " \"type\" : \"fleet\",\n";
            ss << " \"tag\" : \"" << fleet_interface->get_id() << " not found in Information.\",\n}";
        }
        return ss.str();
    }

    virtual std::string
    to_json()
    {

        this->Show();
        typename std::map<uint32_t, std::shared_ptr<PopulationInterfaceBase>>::iterator pit;
        std::vector<uint32_t> pop_ids(this->population_ids->begin(), this->population_ids->end());

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
                pop_strings.push_back(this->population_to_json(pop));
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
                    ss << this->fleets_to_json((FleetInterface *)(*fit).second.get()) << "\n";
                }
                else
                {
                    ss << this->fleets_to_json((FleetInterface *)(*fit).second.get()) << ",\n";
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

        for (it = this->population_ids->begin(); it != this->population_ids->end(); ++it)
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