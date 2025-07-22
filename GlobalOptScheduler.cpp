#include "GlobalOptScheduler.h"

namespace {
    const unsigned int MIN_SCAN_DEFAULT = 30;

    unsigned int getMinScan(const VieVS::Network& network) {
        auto it = std::max_element(
            network.getStations().begin(), 
            network.getStations().end(), 
            [](const VieVS::Station& a, const VieVS::Station& b) {
                return a.getPARA().minScan < b.getPARA().minScan;
            });

        return (it == network.getStations().end()) ? \
            MIN_SCAN_DEFAULT : it->getPARA().minScan;
    }
}

namespace VieVS {
#ifdef WITH_GUROBI
    void GlobalOptScheduler::initialize() noexcept {
        minScan_ = getMinScan(network_);
        blockCount_ = TimeSystem::duration / minScan_;

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "src. count: " << sourceList_.getNSrc();
        BOOST_LOG_TRIVIAL( info ) << "sta. count: " << network_.getNSta();
        BOOST_LOG_TRIVIAL( info ) << "scan block duration [s]: " << minScan_;
        BOOST_LOG_TRIVIAL( info ) << "time block count: " << blockCount_;
#else
        std::cout << "[info] src. count: " << sourceList_.getNSrc();
        std::cout << "[info] sta. count: " << network_.getNSta();
        std::cout << "[info] scan block duration [s]: " << minScan_;
        std::cout << "[info] time block count: " << blockCount_;
#endif

        for(const auto src : sourceList_.getSources())
            src2idx_.insert(std::make_pair(src->getId(), src2idx_.size()));
        for(const auto& sta : network_.getStations()) {
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << sta.getName() << " has observed for " << sta.getTotalObservingTime() << " seconds";
#else
            std::cout << "[info] " << sta.getName() << " has observed for " << sta.getTotalObservingTime() << " seconds";
#endif
            sta2idx_.insert(std::make_pair(sta.getId(), sta2idx_.size()));
        }
            
        try {
            env_ = new GRBEnv(true);
            env_->start();
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Started GRB environment";
#else
            std::cout << "[info] Started GRB environment";
#endif
        } catch (GRBException& e) {
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( error ) << "Gurobi Exception (" << e.getErrorCode() << "): " << e.getMessage();
#else
            std::cout << "[error] Gurobi Exception (" << e.getErrorCode() << "): " << e.getMessage();
#endif
        }

        model_ = new GRBModel(*env_);
        
        // add variables
        for(unsigned int i = 0; i < blockCount_; ++i) {
            for(const auto src : sourceList_.getSources()) {
                y_.emplace_back(model_->addVar(0.0, 1.0, 0.0, GRB_BINARY));

                for(const Station& sta : network_.getStations()) {
                    x_.emplace_back(model_->addVar(0.0, 1.0, 0.0, GRB_BINARY));
                }
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Initialized variables: " << x_.size() << "(x), " << y_.size() << "(y)";
#else
        std::cout << "[info] Initialized variables: " << x_.size() << "(x), " << y_.size() << "(y)";
#endif

        // add constraints

        // each station can only observe one source at a time
        for(unsigned int i = 0; i < blockCount_; ++i) {
            for(const Station& sta : network_.getStations()) {
                GRBLinExpr expr;
                for(const auto src : sourceList_.getSources()) 
                    expr += GlobalOptScheduler::getX(i, src, sta, true);

                model_->addConstr(expr <= 1);
            }
        }

        // must be sufficient time to slew between two targets
        for(Station& sta : network_.refStations())
        for(const auto currSrc : sourceList_.getSources()) for(const auto nextSrc : sourceList_.getSources())
        for(unsigned int currT = 0, nextT; currT < blockCount_ - 1; ++currT) for(nextT = currT + 1; nextT < blockCount_; ++nextT) {
            unsigned int slewT = GlobalOptScheduler::calculateSlewTime(sta, currSrc, nextSrc, currT * minScan_, nextT * minScan_);
            if(nextT - currT >= 1 + slewT / minScan_) continue;

            GRBLinExpr expr = GlobalOptScheduler::getX(currT, currSrc, sta, true) + \
                GlobalOptScheduler::getX(nextT, nextSrc, sta, true);

            model_->addConstr(expr <= 1);
        }

        // a valid scan requires >= 2 participating stations
        for(unsigned int i = 0; i < blockCount_; ++i) {
            for(const auto src : sourceList_.getSources()) {
                GRBLinExpr expr;
                for(const Station& sta : network_.getStations())
                    expr += GlobalOptScheduler::getX(i, src, sta, true);

                model_->addConstr(expr <= (GlobalOptScheduler::getY(i, src, true) * network_.getNSta()));
                model_->addConstr(expr >= GlobalOptScheduler::getY(i, src, true) * 2);
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added constraints";
#else
        std::cout << "[info] Added constraints";
#endif

        min_ = model_->addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER);
        for(const Station& sta : network_.getStations()) {
            GRBLinExpr expr;
            for(unsigned int i = 0; i < blockCount_; ++i) {
                for(const auto src : sourceList_.getSources()) {
                    expr += GlobalOptScheduler::getX(i, src, sta, true);
                }
            }

            model_->addConstr(min_ <= expr);
        }

        // objective function
        GRBLinExpr obj { min_ };
        model_->setObjective(obj, GRB_MAXIMIZE);

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Configured objective function";
#else
        std::cout << "[info] Configured objective function";
#endif
    }

    const GRBVar& GlobalOptScheduler::getX(unsigned int t, 
        const std::shared_ptr<const AbstractSource> src, const Station& sta, const bool tIdx) {
        
        if(!tIdx) t /= TimeSystem::duration * getMinScan(network_);
        auto idx = t * sourceList_.getNSrc() * network_.getNSta() + \
            src2idx_[src->getId()] * network_.getNSta() + sta2idx_[sta.getId()];

        return x_[idx];
    }

    const GRBVar& GlobalOptScheduler::getY(unsigned int t, 
        const std::shared_ptr<const AbstractSource> src, const bool tIdx) {
        
        if(!tIdx) t /= TimeSystem::duration * getMinScan(network_);
        auto idx = t * sourceList_.getNSrc() + src2idx_[src->getId()];

        return y_[idx];
    }

    unsigned int GlobalOptScheduler::calculateSlewTime(
        Station& sta, 
        const std::shared_ptr<const AbstractSource> currSrc, 
        const std::shared_ptr<const AbstractSource> nextSrc,
        const double currT,
        const double nextT) {
        if(currSrc->getId() == nextSrc->getId()) return 0;
        
        PointingVector currVec(sta.getId(), currSrc->getId());
        PointingVector nextVec(sta.getId(), nextSrc->getId());

        currVec.setTime(currT);
        nextVec.setTime(nextT);
        sta.calcAzEl_rigorous(currSrc, currVec);
        sta.calcAzEl_rigorous(nextSrc, nextVec);

        PointingVector tempVec(nextVec);
        if(!sta.isVisible(tempVec)) return std::numeric_limits<unsigned int>::max();

        return sta.getAntenna().slewTime(currVec, nextVec);
    }
    
    void GlobalOptScheduler::start() noexcept {
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Starting optimization";
#else
        std::cout << "[info] Starting optimization";
#endif
        model_->optimize();

        std::vector<unsigned int> eols(network_.getNSta(), 0);

        for(unsigned int i = 0; i < blockCount_; ++i) {
            std::vector<bool> eolsMask(network_.getNSta(), false);
            std::vector<std::vector<PointingVector>> pvs;
            for(auto src : sourceList_.refSources()) {
                std::vector<PointingVector> pvCurr;
                for(Station& sta : network_.refStations()) {
                    // skip if no observation is made
                    if(GlobalOptScheduler::getX(i, src, sta, true).get(GRB_DoubleAttr_X) < 0.5) continue;
                    
                    // station participates in scan at this timestep
                    PointingVector pv(sta.getId(), src->getId());
                    pv.setTime(i * minScan_);
                    sta.calcAzEl_rigorous(src, pv);
                    pvCurr.push_back(std::move(pv));

                    eolsMask[sta2idx_[sta.getId()]] = true;
                }

                if(pvCurr.size()) {
                    pvs.push_back(std::move(pvCurr));
                }
            }

            for(std::vector<PointingVector>& pvSub : pvs) {
                if(pvSub.empty()) continue;
                std::vector<unsigned int> eolsCurr = std::vector<unsigned int>(pvSub.size(), i * minScan_);
                // for(const PointingVector& pv : pvSub) eolsCurr.emplace_back(eols[sta2idx_[pv.getStaid()]]);

                std::vector<PointingVector> pvSubEnd(pvSub);
                for(PointingVector& pvCurr : pvSubEnd) {
                    pvCurr.setTime(pvCurr.getTime() + minScan_);
                    network_.refStation(pvCurr.getStaid()).calcAzEl_rigorous(sourceList_.getSource(pvCurr.getSrcid()), pvCurr);
                }

                Scan scanCurr(pvSub, eolsCurr, Scan::ScanType::standard);
                
                scanCurr.referenceTime().setObservingStarts(i * minScan_);
                scanCurr.setFixedScanDuration(minScan_);
                scanCurr.setPointingVectorsEndtime(pvSubEnd);

                sourceList_.refSource(scanCurr.getSourceId())->update(scanCurr.getNSta(), scanCurr.getNObs(), minScan_, true);
                scans_.push_back(scanCurr);
            }

            for(size_t j = 0; j < eolsMask.size(); ++j) if(eolsMask[j]) eols[j] = i * minScan_;
        }
    }
#else
    void GlobalOptScheduler::initialize() noexcept { /* STUB */ }
    void GlobalOptScheduler::start() noexcept {
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Running standard Scheduler instead of ILP program";
        BOOST_LOG_TRIVIAL( info ) << "This might be because Gurobi wasn't found during configuration";
#else
        std::cout << "[info] Running standard Scheduler instead of ILP program";
        std::cout << "[info] Gurobi was probably not found during CMake configuration";
#endif
        Scheduler::start();
    }
#endif

    Subcon GlobalOptScheduler::createSubcon( const std::shared_ptr<Subnetting> &subnetting, Scan::ScanType type,
                         const boost::optional<StationEndposition> &endposition ) noexcept {
        return Scheduler::createSubcon(subnetting, type, endposition);
    }

    Subcon GlobalOptScheduler::allVisibleScans( Scan::ScanType type, const boost::optional<StationEndposition> &endposition,
                            bool doNotObserveSourcesWithinMinRepeat ) noexcept {
        return Scheduler::allVisibleScans(type, endposition, doNotObserveSourcesWithinMinRepeat);
    }

    void GlobalOptScheduler::update( Scan &scan, std::ofstream &of ) noexcept {
        return Scheduler::update(scan, of);
    }

    void GlobalOptScheduler::consideredUpdate( unsigned long n1scans, unsigned long n2scans, int depth, std::ofstream &of ) noexcept {
        return Scheduler::consideredUpdate(n1scans, n2scans, depth, of);
    }

    void GlobalOptScheduler::statistics( std::ofstream &of ) {
        return Scheduler::statistics(of);
    }

    void GlobalOptScheduler::highImpactScans( HighImpactScanDescriptor &himp, std::ofstream &of ) {
        return Scheduler::highImpactScans(himp, of);
    }

    void GlobalOptScheduler::calibratorBlocks( std::ofstream &of ) {
        return Scheduler::calibratorBlocks(of);
    }

    void GlobalOptScheduler::parallacticAngleBlocks( std::ofstream &of ) {
        return Scheduler::parallacticAngleBlocks(of);
    }

    void GlobalOptScheduler::differentialParallacticAngleBlocks( std::ofstream &of ) {
        return Scheduler::differentialParallacticAngleBlocks(of);
    }

    bool GlobalOptScheduler::checkAndStatistics( std::ofstream &of ) noexcept {
        return Scheduler::checkAndStatistics(of);
    }

    void GlobalOptScheduler::checkSatelliteAvoidance() {
        return Scheduler::checkSatelliteAvoidance();
    }
} // namespace VieVS
