#include "GlobalOptScheduler.h"
#include "Misc/TimeSystem.h"
#include <iostream>

namespace VieVS {
#ifdef WITH_GUROBI
    void GlobalOptScheduler::start() noexcept {
        std::cout << "[sourceList_]" << std::endl;
        for(const auto& src : sourceList_.refSources()) std::cout << src->getName() << std::endl;

        std::cout << "[network_]" << std::endl; 
        for(const auto& sta : network_.refStations()) std::cout << sta.getName() << std::endl;

        std::cout << "[TimeSystem]" << std::endl;
        std::cout << "mjdStart: " << TimeSystem::mjdStart << std::endl;
        std::cout << "startTime: " << TimeSystem::startTime << std::endl;
        std::cout << "endTime: " << TimeSystem::endTime << std::endl;
        std::cout << "duration: " << TimeSystem::duration << std::endl;
    }
#endif
} // namespace VieVS

/* 
 * int version_;                 ///< version
 * std::string path_;            ///< path to VieSchedpp.xml directory
 * 
 * boost::property_tree::ptree xml_;  ///< content of VieSchedpp.xml file
 * 
 * SourceList sourceList_;                                       ///< session source list
 * Network network_;                                             ///< station network
 * std::shared_ptr<const ObservingMode> obsModes_ = nullptr;     ///< observing modes
 * std::shared_ptr<const Mode> currentObservingMode_ = nullptr;  ///< current observing mode
 * std::vector<Scan> scans_;                                     ///< all scans in schedule
 * 
 * Parameters parameters_;        ///< general scheduling parameters
 * PreCalculated preCalculated_;  ///< pre calculated values
 * 
 * unsigned long nSingleScansConsidered = 0;      ///< considered single source scans
 * unsigned long nSubnettingScansConsidered = 0;  ///< considered subnetting scans
 * unsigned long nObservationsConsidered = 0;     ///< considered baselines
 * 
 * boost::optional<HighImpactScanDescriptor> himp_;                          ///< high impact scan descriptor
 * std::vector<CalibratorBlock> calib_;                                      ///< fringeFinder impact scan descriptor
 * boost::optional<MultiScheduling::Parameters> multiSchedulingParameters_;  ///< multi scheduling paramters
 */