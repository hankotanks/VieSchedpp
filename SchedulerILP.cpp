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
        std::vector<Scan> scansOptimal = model_->optimize(scans_);
        if(scansOptimal.empty()) {
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Reverting to previous solution";
#else
            std::cout << "[info] Reverting to previous solution";
#endif
            return;
        }

        scans_.clear();
        scans_.insert(scans_.end(), scansOptimal.begin(), scansOptimal.end());

        sortSchedule( Timestamp::start );

        // TODO: remove this debug block
#if 0
        for(const Station& s : network_.getStations()) {
            size_t count = 0;
            std::cout << s.getName() << ": ";
            for(const Scan& scan : scans_) {
                if(scan.findIdxOfStationId(s.getId())) {
                    std::cout << scan.getSourceId() << "[" << scan.getTimes().getScanTime() / 60 << "], ";
                    ++count;
                }
            }
            std::cout << std::endl;
            std::cout << s.getName() << ": " << count << std::endl;
        }
        std::cout << "length: " << scans_.size() << std::endl;
#endif 

        updateObservingTimes();

        // check if there was an error during the session
        std::ofstream of;
        of.open( path_ + "ilp_statistics_check.log" );
        if ( !checkAndStatistics( of ) ) {
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( error ) << "error(s) found while checking the schedule";
#else
            cout << "[error] error(s) found while checking the schedule";
#endif
        }
        of.close();

        sortSchedule( Timestamp::start );
    }

    void SchedulerILP::update( Scan &scan, std::ofstream &of ) noexcept {
        Scheduler::update(scan, of);
    }

    void SchedulerILP::statistics( std::ofstream &of ) {
        Scheduler::statistics(of);
    }
} // namespace VieVS
