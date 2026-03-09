#include "SchedulerILP.h"

#include <exception>
#include <memory>
#include <stdexcept>

#include "Output/Output.h"
#include "Scheduler.h"
#include "Source/AbstractSource.h"

namespace {
    unsigned int getBlockLength(const VieVS::Network& network) {
        auto it = std::max_element(
            network.getStations().begin(), 
            network.getStations().end(), 
            [](const VieVS::Station& a, const VieVS::Station& b) {
                return a.getPARA().minScan < b.getPARA().minScan;
            });

        return it->getPARA().minScan;
    }
} // private

namespace VieVS {
    void SchedulerILP::initialize() noexcept {
        unsigned int blockLength = getBlockLength(network_);
        
        try {
            unsigned int windowLength;
            windowLength = xml_.get<unsigned int>( "VieSchedpp.general.ilpOptimizationWindow" );
            if(windowLength < 3 * blockLength) {
                throw std::runtime_error("Length of optimization window must be >= 3 times the minimum scan length");
            }
            // initialize the model
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

        for(SkyCoverage& sky : network_.refSkyCoverages()) {
            sky.calculateSkyCoverageScores();
        }

        // save initial solution for later output
        VieVS::Scheduler initial(this);

        // compute optimal scans from the ILP model
        std::vector<Scan> scansOptimal = model_->optimize(scans_);
        if(scansOptimal.empty()) {
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Reverting to previous solution";
#else
            std::cout << "[info] Reverting to previous solution";
#endif
            return;
        } else {
            VieVS::Output output(initial);
            output.writeStatistics(statisticsOf_);
            ++version_;
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Output initial solution's statistics";
#else
            std::cout << "[info] Output initial solution's statistics";
#endif
        }

        // replace old scans with them
        scans_.clear();
        scans_.insert(scans_.end(), scansOptimal.begin(), scansOptimal.end());
        sortSchedule( Timestamp::start );

#if 0
        for(const Scan& scan : scans_) {
            std::cout << "[" << scan.getId() << "] source: " << scan.getSourceId();
            for(unsigned long sId : scan.getStationIds()) {
                if(auto idx = scan.findIdxOfStationId(sId)) {
                    std::cout << ", " << network_.getStation(sId).getName() << 
                        " (" << scan.getTimes().getObservingTime(*idx) << 
                        " for " << scan.getTimes().getObservingDuration(*idx) << ")";
                }
            }
            std::cout << std::endl;
        }
#endif

        // clear all the statistics from the greedy run
        for(Station& s : network_.refStations()) {
            s.setStatistics(Station::Statistics());
        }
        for(std::shared_ptr<VieVS::AbstractSource> q : sourceList_.refSources()) {
            q->setStatistics(AbstractSource::Statistics());
        }

        const std::map<unsigned long, unsigned long>& staids2skyCoverageId = network_.getStaid2skyCoverageId();
        for(const Station& s : network_.getStations()) {
            SkyCoverage& sky = network_.refSkyCoverage(staids2skyCoverageId.at(s.getId()));
            sky.clearObservations();

            for(const Scan& scan : scans_) {
                if(auto i = scan.findIdxOfStationId(s.getId())) {
                    const PointingVector& pv0 = scan.getPointingVector(*i);
                    sky.update(pv0);
                }   
            }

            sky.calculateSkyCoverageScores();
        }

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
