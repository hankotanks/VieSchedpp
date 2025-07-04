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
 * @file VieSchedpp.h
 * @brief class VieSchedpp
 *
 * @author Matthias Schartner
 * @date 02.07.2018
 */

#ifndef VIEVS_SCHEDULER_H
#define VIEVS_SCHEDULER_H

// clang-format off
#include "Simulator/Simulator.h"
#include "Simulator/Solver.h"
// clang-format on

#include <boost/format.hpp>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

#include "Initializer.h"
#include "Misc/CalibratorBlock.h"
#include "Misc/HighImpactScanDescriptor.h"
#include "ObservingMode/Mode.h"
#include "Output/Output.h"
#include "Scheduler.h"
#include "XML/ParameterSettings.h"
#include "Scheduler.h"
#include "GlobalOptScheduler.h"
#ifdef VIESCHEDPP_LOG
#include <boost/log/attributes/named_scope.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace VieVS {

/**
 * @class VieSchedpp
 * @brief this is the VieSched++ core class
 *
 * this class starts the scheduling process.
 * It monitores the creation of all schedules and all threads.
 *
 * @author Matthias Schartner
 * @date 02.07.2018
 */

class VieSchedpp {
   public:
    /**
     * @brief constructor
     * @author Matthias Schartner
     */
    VieSchedpp() = default;


    /**
     * @brief constructor
     * @author Matthias Schartner
     *
     * @param inputFile path and name of input file (usually called VieSchedpp.xml)
     */
    explicit VieSchedpp( const std::string &inputFile );


    /**
     * @brief start VieSched++
     * @author Matthias Schartner
     */
    void run(const bool ILP = false);


   private:
    std::string inputFile_;            ///< VieSchedpp.xml file
    std::string path_;                 ///< path to VieSchedpp.xml file
    std::string sessionName_;          ///< session name
    boost::property_tree::ptree xml_;  ///< content of VieSchedpp.xml file
    unsigned long nsta_;

    SkdCatalogReader skdCatalogs_;                                          ///< sked catalogs
    std::vector<VieVS::MultiScheduling::Parameters> multiSchedParameters_;  ///< list of all multi scheduling parameters

    /**
     * @brief read sked catalogs
     * @author Matthias Schartner
     */
    void readSkdCatalogs();


    /**
     * @brief initialize parallel processing
     * @author Matthias Schartner
     */
    void multiCoreSetup();


    /**
     * @brief initialize log files
     * @author Matthias Schartner
     */
    void init_log();

    /**
     * @brief summarize simulation result
     * @author Matthias Schartner
     *
     * @param output flag if output should be printed to console
     * @return score per version
     */
    std::map<int, double> summarizeSimulationResult( const Network &network, const SourceList &srclist,
                                                     bool output = true );

    /**
     * @brief get priority values from xml file
     * @author Matthias Schartner
     *
     * @param type output type (mfe or rep)
     * @param header list of header file entries
     * @return list of tuples (first = name, second = index, third = priority value) of priority elements
     */
    std::vector<std::tuple<std::string, int, double>> getPriorityCoefficients( const std::string &type,
                                                                               const Network &network,
                                                                               const SourceList &srclist,
                                                                               const std::vector<std::string> &header );

    /**
     * @brief list best schedules based on simulations
     * @author Matthias Schartner
     *
     * @param of out file stream
     * @param type output type (mfe or rep)
     * @param storage values per schedule (first = version, second = list of values)
     * @param priorityLookup list of tuples (first = name, second = index, third = priority value) of priority elements
     * @return list of scores per session (key = version number, value = score)
     */
    std::map<int, double> listBest( std::ofstream &of, const std::string &type,
                                    const std::map<int, std::vector<double>> &storage,
                                    const std::vector<std::tuple<std::string, int, double>> &priorityLookup );

    /**
     * @brief output best schedules based on simulations
     * @author Matthias Schartner
     *
     * @param mfe_costs mean formal error costs
     * @param rep_costs repeatability costs
     * @param storage values per schedule (first = version, second = list of values)
     * @param output flag if output should be printed to console
     * @return score per version
     */
    std::map<int, double> printRecommendation( const std::map<int, double> &mfe_costs,
                                               const std::map<int, double> &rep_costs,
                                               const std::map<int, std::vector<double>> &storage, bool output = true );
};
}  // namespace VieVS

#endif  // VIEVS_SCHEDULER_H
