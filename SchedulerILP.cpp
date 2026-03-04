#include "SchedulerILP.h"

#include <exception>
#include <memory>

#include "Source/AbstractSource.h"

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

        // compute optimal scans from the ILP model
        std::vector<Scan> scansOptimal = model_->optimize(scans_);
        if(scansOptimal.empty()) {
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Reverting to previous solution";
#else
            std::cout << "[info] Reverting to previous solution";
#endif
            return;
        }

        // replace old scans with them
        scans_.clear();
        scans_.insert(scans_.end(), scansOptimal.begin(), scansOptimal.end());

        // clear all the statistics from the greedy run
        for(Station& s : network_.refStations()) {
            s.setStatistics(Station::Statistics());
        }
        for(std::shared_ptr<VieVS::AbstractSource> q : sourceList_.refSources()) {
            q->setStatistics(AbstractSource::Statistics());
        }

        sortSchedule( Timestamp::start );

        std::string fileName = getName() + "_iteration_ilp.txt";
        std::ofstream of;
        of.open( path_ + fileName );

        // TODO: This currently causes rare slew time violations
#if 0
        if ( parameters_.idleToObservingTime ) {
            switch ( ScanTimes::getAlignmentAnchor() ) {
                case ScanTimes::AlignmentAnchor::start: {
                    idleToScanTime( Timestamp::end, of );
                    break;
                }
                case ScanTimes::AlignmentAnchor::end: {
                    idleToScanTime( Timestamp::start, of );
                    break;
                }
                case ScanTimes::AlignmentAnchor::individual: {
                    idleToScanTime( Timestamp::end, of );
                    idleToScanTime( Timestamp::start, of );
                    break;
                }
                default:
                    break;
            }
        }
#endif

        updateObservingTimes();

        // check if there was an error during the session (and append to the output log)
        if ( !checkAndStatistics( of ) ) {
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( error ) << "error(s) found while checking the schedule";
#else
            cout << "[error] error(s) found while checking the schedule";
#endif
        }
        of.close();
    }

    void SchedulerILP::update( Scan &scan, std::ofstream &of ) noexcept {
        Scheduler::update(scan, of);
    }

    void SchedulerILP::statistics( std::ofstream &of ) {
        Scheduler::statistics(of);
    }
} // namespace VieVS
