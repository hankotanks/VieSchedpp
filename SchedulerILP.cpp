#include "SchedulerILP.h"

#include <exception>
#include <memory>

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

#define DEFAULT_WINDOW_BLOCK_COUNT 10

namespace VieVS {
    void SchedulerILP::initialize() noexcept {
        unsigned int blockLength = getBlockLength(network_);
        // TODO: This should be tunable
        unsigned int windowLength = blockLength * DEFAULT_WINDOW_BLOCK_COUNT;
        try {
            model_ = new Model(network_, sourceList_, blockLength, windowLength);
        } catch(GRBException& e) {
            (void) e;
            model_ = nullptr;
        } catch(const std::exception& e) {
            model_ = nullptr;
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( error ) << e.what();
#else
            std::cout << "[error] " << e.what();
#endif  
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Proceeding with default scheduler";
#else
            std::cout << "[info] Proceeding with default scheduler";
#endif  
        }
    }
    
    void SchedulerILP::start() noexcept {
        Scheduler::start();
        model_->optimize(scans_);
    }

    void SchedulerILP::update( Scan &scan, std::ofstream &of ) noexcept {
        Scheduler::update(scan, of);
    }

    void SchedulerILP::statistics( std::ofstream &of ) {
        Scheduler::statistics(of);
    }
} // namespace VieVS
