#ifndef SCHEDULER_ILP_H
#define SCHEDULER_ILP_H

#include <boost/date_time.hpp>
#include <boost/optional.hpp>
#include <tuple>
#include <utility>
#include <vector>

#include "Algorithm/FocusCorners.h"
#include "Initializer.h"
#include "Misc/Constants.h"
#include "Misc/StationEndposition.h"
#include "Misc/Subnetting.h"
#include "Scan/Subcon.h"
#include "Scheduler.h"
#include "Station/Network.h"


namespace VieVS {
/**
 * @class Scheduler
 * @brief this is the VLBI scheduling class which is responsible for the scan selection and the creation of the
 * schedule
 *
 * @author Matthias Schartner
 * @date 28.06.2017
 */
class SchedulerILP : public Scheduler {
    friend class Output;
public:
    /**
    * @brief constructor
    * @author Matthias Schartner
    *
    * @param init initializer
    * @param path path to VieSchedpp.xml file
    * @param fname file name
    */
    SchedulerILP( Initializer &init, std::string path, std::string fname ) : 
        Scheduler(init, path, fname) { /* STUB */ }


    /**
    * @brief constructor
    * @author Matthias Schartner
    *
    * @param name session name
    * @param path session path
    * @param network_ station network
    * @param sourceList source list
    * @param scans list of scans
    * @param xml VieSchedpp.xml file
    */
    SchedulerILP(std::string name, std::string path, Network network, SourceList sourceList,
            std::vector<Scan> scans,
            boost::property_tree::ptree xml, std::shared_ptr<ObservingMode> obsModes = nullptr ) : 
        Scheduler(name, path, network, sourceList, scans, xml, obsModes) { /* STUB */ }
};
}
#endif // SCHEDULER_ILP_H