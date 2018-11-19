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
 * @file Scheduler.h
 * @brief class Scheduler
 *
 * @author Matthias Schartner
 * @date 29.06.2017
 */

#ifndef SCHEDULER_H
#define SCHEDULER_H
#include <vector>
#include <boost/date_time.hpp>
#include <utility>
#include <tuple>
#include <boost/optional.hpp>

#include "Station/Network.h"
#include "Initializer.h"
#include "Scan/Subcon.h"
#include "Misc/Constants.h"
#include "Misc/StationEndposition.h"
#include "Misc/Subnetting.h"

namespace VieVS{
    /**
     * @class Scheduler
     * @brief this is the VLBI scheduling class which is responsible for the scan selection and the creation of the
     * schedule
     *
     * @author Matthias Schartner
     * @date 28.06.2017
     */
    class Scheduler: public VieVS_NamedObject {
        friend class Output;

    public:

        /**
         * @brief general parameters used for scheduling
         * @author Matthias Schartner
         */
        struct Parameters {
            boost::optional<Subnetting> subnetting = boost::none; ///< subnetting parameters
            double subnettingMinNSta = 0.60; /// minimum number of subnetting station percent (deprecated)
            bool fillinmodeDuringScanSelection = true; ///< flag if fillin modes are allowed
            bool fillinmodeInfluenceOnSchedule = true; ///< fillin modes scans influence schedule if set to true
            bool fillinmodeAPosteriori = false; ///< fillin mode a posteriori

            bool idleToObservingTime = true; ///< idle to observing time

            bool andAsConditionCombination = true; ///< condition combination model
            unsigned int currentIteration = 0; ///< current iteration number
            unsigned int maxNumberOfIterations = 999; ///< max number of iterations
            unsigned int numberOfGentleSourceReductions = 0; ///< number of gentle source reductions
            unsigned int minNumberOfSourcesToReduce = 0; ///< min number of sources to reduce

            bool writeSkyCoverageData = false; ///< flag if sky coverage data should be printed to file
        };

        /**
         * @brief pre calculated values (deprecated)
         * @author Matthias Schartner
         */
        struct PreCalculated{
            std::vector<std::vector<int>> subnettingSrcIds; ///< list of all available second sources in subnetting
        };

        /**
         * @brief constructor
         * @author Matthias Schartner
         *
         * @param init initializer
         * @param path path to VieSchedpp.xml file
         * @param session name
         */
        Scheduler(Initializer &init, std::string path, std::string fname);

        /**
         * @brief constructor
         * @author Matthias Schartner
         *
         * @param name session name
         * @param network_ station network
         * @param sources source list
         * @param scans list of scans
         * @param xml VieSchedpp.xml file
         */
        Scheduler(std::string name, Network network_, std::vector<Source> sources,
                  std::vector<Scan> scans, boost::property_tree::ptree xml);

        /**
         * @brief main function that starts the scheduling
         * @author Matthias Schartner
         */
        void start() noexcept;

        /**
         * @brief this function creates a subcon with all scans, times and scores
         * @author Matthias Schartner
         *
         * @param subnetting true if subnetting is allowed, false otherwise
         * @param type scan type
         * @param endposition required endposition
         * @return subcon with all information
         */
        Subcon createSubcon(const boost::optional<Subnetting> &subnetting, Scan::ScanType type,
                            const boost::optional<StationEndposition> &endposition = boost::none) noexcept;

        /**
         * @brief constructs all visible scans
         * @author Matthias Schartner
         *
         * @param endposition required endposition
         * @return subcon with all visible single source scans
         */
        Subcon allVisibleScans(Scan::ScanType type, const boost::optional<StationEndposition> &endposition= boost::none) noexcept;


        /**
         * @brief updates the selected next scans to the schedule
         * @author Matthias Schartner
         *
         * @param scan best possible next scans
         * @param of outstream file object
         */
        void update(Scan &scan, std::ofstream &of) noexcept;

        /**
         * @brief updates and prints the number of all considered scans
         * @author Matthias Schartner
         *
         * @param n1scans number of single source scans
         * @param n2scans number of subnetting scans
         * @param depth recursion depth
         * @param of outstream file object
         */
        void consideredUpdate(unsigned long n1scans, unsigned long n2scans, int depth, std::ofstream &of) noexcept;

        /**
         * @brief statistics output
         * @author Matthias Schartner
         *
         * @param of output stream
         */
        void statistics(std::ofstream &of);

        /**
         * @brief schedule high impact scans
         * @author Matthias Schartner
         *
         * @param himp high impact scan descriptor
         * @param of outstream object
         */
        void highImpactScans(HighImpactScanDescriptor &himp, std::ofstream &of);

        /**
         * @brief checks the schedule with an independend methode
         * @author Matthias Schartner
         *
         * @param of outstream file object
         */
        bool checkAndStatistics(std::ofstream &of) noexcept;

        /**
         * @brief get all sources
         * @author Matthias Schartner
         *
         * @return all sources
         */
        const std::vector<Source> &getSources() const noexcept{
            return sources_;
        }

        /**
         * @brief get station network
         * @author Matthias Schartner
         *
         * @return station network
         */
        const Network &getNetwork() const noexcept{
            return network_;
        }

        /**
         * @brief get all scans
         * @author Matthias Schartner
         *
         * @return all scans
         */
        const std::vector<Scan> &getScans() const noexcept{
            return scans_;
        }


    private:
        static unsigned long nextId; ///< next id for this object type
        std::string path_; ///< path to VieSchedpp.xml directory

        boost::property_tree::ptree xml_; ///< content of VieSchedpp.xml file

        std::vector<Source> sources_; ///< all sources
        Network network_; ///< station network
        std::shared_ptr<const ObservingMode> obsModes_ = nullptr; ///< observing modes
        std::shared_ptr<const Mode> currentObservingMode_ = nullptr; ///< current observing mode
        std::vector<Scan> scans_; ///< all scans in schedule

        Parameters parameters_; ///< general scheduling parameters
        PreCalculated preCalculated_; ///< pre calculated values

        unsigned long nSingleScansConsidered = 0; ///< considered single source scans
        unsigned long nSubnettingScansConsidered = 0; ///< considered subnetting scans
        unsigned long nObservationsConsidered = 0; ///< considered baselines

        boost::optional<HighImpactScanDescriptor> himp_;
        boost::optional<MultiScheduling::Parameters> multiSchedulingParameters_;

        /**
         * @brief start recursive scan selection
         * @author Matthias Schartner
         *
         * @param endTime end time of this recursion
         * @param of outstream object
         * @param type scan type
         * @param opt_endposition required endposition
         * @param subcon precalculated subcon
         * @param depth recursion depth
         */
        void startScanSelection(unsigned int endTime, std::ofstream &of, Scan::ScanType type,
                                boost::optional<StationEndposition> &opt_endposition, boost::optional<Subcon> &subcon,
                                int depth);


        /**
         * @brief checks if some parameters need to be changed
         * @author Matthias Schartner
         *
         * @param time current time in seconds since start
         * @param output flag if output to log file is required
         * @param of outstream file object
         * @param tagalong schedule tagalong scans
         * @return true if a hard break was found
         */
        bool checkForNewEvents(unsigned int time, bool output, std::ofstream &of, bool tagalong) noexcept;

        /**
         * @brief output of source overview
         * @author Matthias Schartner
         *
         * @param outstream object
         */
        void listSourceOverview(std::ofstream &of) noexcept;

        /**
         * @brief start tagalong mode
         * @author Matthias Schartner
         *
         * @param station tagalong station
         * @param of outstream object
         */
        void startTagelongMode(Station &station, std::ofstream &of);

        /**
         * @brief check optimization conditions
         * @author Matthias Schartner
         *
         * @param of outstream object
         */
        bool checkOptimizationConditions(std::ofstream &of);

        /**
         * @brief change station availability
         * @author Matthias Schartner
         *
         * @param endposition required endpositions
         * @param change change type
         */
        void changeStationAvailability(const boost::optional<StationEndposition> &endposition,
                                   StationEndposition::change change);

        /**
         * @brief start scan selection between scans
         * @author Matthias Schartner
         *
         * for example fillin mode a posteriori
         *
         * @param duration duration (deprecated?)
         * @param of outstream object
         * @param type scan type
         * @param output output flag
         * @param ignoreTagalong ignore tagalong modes flag
         */
        void startScanSelectionBetweenScans(unsigned int duration, std::ofstream &of, Scan::ScanType type, bool output=false, bool ignoreTagalong=false);

        /**
         * @brief reset all events to time zero
         * @author Matthias Schartner
         *
         * @param of outstream object
         */
        void resetAllEvents(std::ofstream &of);

        /**
         * @brief ignore tagalong parameter
         * @author Matthias Schartner
         */
        void ignoreTagalongParameter();

        /**
         * @brief transform idle time to observing time
         * @author Matthias Schartner
         *
         * @param ts time stamp
         * @param of outstrem object
         */
        void idleToScanTime(Timestamp ts, std::ofstream &of);

        /**
         * @brief sort schedule
         * @author Matthias Schartner
         *
         * @param ts time stamp
         */
        void sortSchedule(Timestamp ts = Timestamp::start);

        /**
         * @brief sort schedule based on station
         * @author Matthias Schartner
         *
         * @param staid station id
         * @param ts time stamp
         */
        void sortSchedule(unsigned long staid, Timestamp ts = Timestamp::start);

        /**
         * @brief calibrator update
         * @author Matthias Schartner
         *
         * @param bestScans scheduled scans
         * @param prevHighElevationScores previouse high elevation scores
         * @param prevLowElevationScores previouse low elevation scores
         * @param highestElevations highest elevation scores
         * @param lowestElevations lowest elevation scores
         * @return true if no more calibrator scans are needed, otherwise false
         */
        bool calibratorUpdate(const std::vector<Scan> &bestScans,
                              std::vector<double> &prevHighElevationScores, std::vector<double> &prevLowElevationScores,
                              std::vector<double> &highestElevations, std::vector<double> &lowestElevations);

        /**
         * @brief write calibrator statistics
         * @author Matthias Schartner
         *
         * @param of output stream object
         * @param highestElevations highest elevations scheduled so far in calibrator block
         * @param lowestElevations lowest elevations scheduled so far in calibrator block
         */
        void writeCalibratorStatistics(std::ofstream &of, std::vector<double> &highestElevations,
                                       std::vector<double> &lowestElevations);

        /**
         * @brief write calibrator block header
         * @author Matthias Schartner
         *
         * @param of outstream object
         */
        void writeCalibratorHeader(std::ofstream &of);

    };
}
#endif /* SCHEDULER_H */

