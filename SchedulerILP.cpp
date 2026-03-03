#include "SchedulerILP.h"
#include <exception>

namespace {
    const unsigned int DEFAULT_BLOCK_LENGTH = 30;

    unsigned int getBlockLength(const VieVS::Network& network) {
        auto it = std::max_element(
            network.getStations().begin(), 
            network.getStations().end(), 
            [](const VieVS::Station& a, const VieVS::Station& b) {
                return a.getPARA().minScan < b.getPARA().minScan;
            });

        return (it == network.getStations().end()) ? DEFAULT_BLOCK_LENGTH : it->getPARA().minScan;
    }
} // private

namespace VieVS {
    void SchedulerILP::initialize() noexcept {
        try {
            model_ = new Model(network_, sourceList_, getBlockLength(network_));
        } catch(...) {
            model_ = nullptr;
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Proceeding with default scheduler";
#else
            std::cout << "[info] Proceeding with default scheduler";
#endif  
        }
    }
    
    void SchedulerILP::start() noexcept {
        Scheduler::start();
    }

    void SchedulerILP::update( Scan &scan, std::ofstream &of ) noexcept {
        Scheduler::update(scan, of);
    }

    void SchedulerILP::statistics( std::ofstream &of ) {
        Scheduler::statistics(of);
    }
} // namespace VieVS
