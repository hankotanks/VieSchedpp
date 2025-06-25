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
#ifdef WITH_GUROBI
#include <gurobi_c++.h>
#endif
#include <memory>
#include <tuple>
#include <utility>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <numeric>
#include <algorithm>

#include "Initializer.h"
#include "Misc/Constants.h"
#include "Misc/StationEndposition.h"
#include "Misc/Subnetting.h"
#include "Misc/TimeSystem.h"
#include "Misc/VieVS_NamedObject.h"
#include "Scan/Subcon.h"
#include "Source/AbstractSource.h"
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
private:
    void initialize() noexcept;
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
        Scheduler(init, path, fname) { GlobalOptScheduler::initialize(); }

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
        Scheduler(name, path, network, sourceList, scans, xml, obsModes) { 
            GlobalOptScheduler::initialize(); 
        }

    /**
     * @brief destructor
     * @author Hank Lewis
    */
    ~GlobalOptScheduler() {
#ifdef WITH_GUROBI
        delete env_;
        delete model_;
#endif
    }

    /**
     * @brief main function that starts the scheduling
     * @author Hank Lewis
     */
    void start() noexcept override;

    /**
     * @brief this function creates a subcon with all scans, times and scores
     * @author Hank Lewis
     *
     * @param subnetting true if subnetting is allowed, false otherwise
     * @param type scan type
     * @param endposition required endposition
     * @return subcon with all information
     */
    Subcon createSubcon( const std::shared_ptr<Subnetting> &subnetting, Scan::ScanType type,
                        const boost::optional<StationEndposition> &endposition = boost::none ) noexcept override;


    /**
     * @brief constructs all visible scans
     * @author Hank Lewis
     *
     * @param type scan type
     * @param endposition required endposition
     * @param doNotObserveSourcesWithinMinRepeat consider scans (with reduced weight) if they are within min repeat time
     * @return subcon with all visible single source scans
     */
    Subcon allVisibleScans( Scan::ScanType type, const boost::optional<StationEndposition> &endposition = boost::none,
                        bool doNotObserveSourcesWithinMinRepeat = true ) noexcept override;


    /**
     * @brief updates the selected next scans to the schedule
     * @author Hank Lewis
     *
     * @param scan best possible next scans
     * @param of outstream file object
     */
    void update( Scan &scan, std::ofstream &of ) noexcept override;


    /**
     * @brief updates and prints the number of all considered scans
     * @author Hank Lewis
     *
     * @param n1scans number of single source scans
     * @param n2scans number of subnetting scans
     * @param depth recursion depth
     * @param of outstream file object
     */
    void consideredUpdate( unsigned long n1scans, unsigned long n2scans, int depth, std::ofstream &of ) noexcept override;


    /**
     * @brief statistics output
     * @author Hank Lewis
     *
     * @param of output stream
     */
     void statistics( std::ofstream &of ) override;


    /**
     * @brief schedule high impact scans
     * @author Hank Lewis
     *
     * @param himp high impact scan descriptor
     * @param of outstream object
     */
    void highImpactScans( HighImpactScanDescriptor &himp, std::ofstream &of ) override;

    /**
     * @brief schedule fringeFinder blocks
     * @author Hank Lewis
     *
     * @param of outstream object
     */
    void calibratorBlocks( std::ofstream &of ) override;

    /**
     * @brief schedule fringeFinder blocks
     * @author Hank Lewis
     *
     * @param of outstream object
     */
    void parallacticAngleBlocks( std::ofstream &of ) override;

    /**
     * @brief schedule fringeFinder blocks
     * @author Hank Lewis
     *
     * @param of outstream object
     */
    void differentialParallacticAngleBlocks( std::ofstream &of ) override;


    /**
     * @brief checks the schedule with an independend methode
     * @author Hank Lewis
     *
     * @param of outstream file object
     */
    bool checkAndStatistics( std::ofstream &of ) noexcept override;

    /**
     * @brief check if there is a satellite too close to a scan
     * @author Hank Lewis
     */
    void checkSatelliteAvoidance() override;
    
#ifdef WITH_GUROBI
private:
    const GRBVar& getX(const unsigned int t, 
        std::shared_ptr<AbstractSource> src, const Station& sta);
    const GRBVar& getY(const unsigned int t, std::shared_ptr<AbstractSource> src);
private:
    std::unordered_map<unsigned long, std::size_t> sta2idx_;
    std::unordered_map<unsigned long, std::size_t> src2idx_;
    GRBEnv* env_ = nullptr;
    GRBModel* model_ = nullptr;
    std::vector<GRBVar> x_;
    std::vector<GRBVar> y_;
#endif
};
}
#endif // GLOBAL_OPT_SCHEDULER_H