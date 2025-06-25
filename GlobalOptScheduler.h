/*
 *  VieSched++ Very Long Baseline Interferometry (VLBI) Scheduling Software
 *  Copyright (C) 2018  Matthias Schartner
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file GlobalOptScheduler.h
 * @brief class GlobalOptScheduler
 *
 * @author Hank Lewis
 * @date 24.05.2025
 */

#ifndef GLOBAL_OPT_SCHEDULER_H
#define GLOBAL_OPT_SCHEDULER_H

#include <boost/date_time.hpp>
#include <boost/optional.hpp>
#include <tuple>
#include <utility>
#include <vector>


#include "Initializer.h"
// #include "Misc/Constants.h"
// #include "Misc/StationEndposition.h"
// #include "Misc/Subnetting.h"
// #include "Scan/Subcon.h"
#include "Station/Network.h"
#include "Scheduler.h"


namespace VieVS {
/**
 * @class Scheduler
 * @brief this is an overload of the original Scheduler class that implements
 * globally optimal scheduling using integer linear programming (ILP)
 *
 * @author Hank Lewis
 * @date 24.05.2025
 */
class GlobalOptScheduler : public Scheduler {
    friend class Output;
public:
    /**
    * @brief constructor
    * @author Hank Lewis
    *
    * @param init initializer
    * @param path path to VieSchedpp.xml file
    * @param fname file name
    */
    GlobalOptScheduler( Initializer &init, std::string path, std::string fname ) : 
        Scheduler(init, path, fname) { /* STUB */ }


    /**
    * @brief constructor
    * @author Hank Lewis
    *
    * @param name session name
    * @param path session path
    * @param network_ station network
    * @param sourceList source list
    * @param scans list of scans
    * @param xml VieSchedpp.xml file
    */
    GlobalOptScheduler( std::string name, std::string path, Network network, SourceList sourceList,
            std::vector<Scan> scans,
            boost::property_tree::ptree xml, std::shared_ptr<ObservingMode> obsModes = nullptr ) : 
        Scheduler(name, path, network, sourceList, scans, xml, obsModes) { /* STUB */ }

    /**
     * @brief main function that starts the scheduling
     * @author Hank Lewis
     */
    void start() noexcept;
};
}
#endif // GLOBAL_OPT_SCHEDULER_H