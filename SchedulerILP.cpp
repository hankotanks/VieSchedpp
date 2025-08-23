#include "SchedulerILP.h"
#include "Misc/Subnetting.h"

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

#define SKY_COVERAGE_CELL_COUNT 13

namespace VieVS {
#ifdef WITH_GUROBI
    void SchedulerILP::initialize() noexcept {
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
            sta2idx_.insert(std::make_pair(sta.getId(), sta2idx_.size()));
            sta2pv0_.insert(std::make_pair(sta.getId(), sta.getCurrentPointingVector()));
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

        for(const Station& sta : network_.getStations()) {
            for(std::size_t i = 0; i < SKY_COVERAGE_CELL_COUNT; ++i)
                z_.emplace_back(model_->addVar(0.0, 1.0, 0.0, GRB_BINARY));
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
                    expr += SchedulerILP::getX(i, src, sta);

                model_->addConstr(expr <= 1);
            }
        }

        // must be sufficient time to slew between two targets
        for(Station& sta : network_.refStations())
        for(const auto currSrc : sourceList_.getSources()) for(const auto nextSrc : sourceList_.getSources())
        for(unsigned int currT = 0, nextT; currT < blockCount_ - 1; ++currT) for(nextT = currT + 1; nextT < blockCount_; ++nextT) {
            unsigned int slewT = SchedulerILP::calculateSlewTime(sta, currSrc, nextSrc, currT * minScan_, nextT * minScan_);
            if(nextT - currT >= 1 + slewT / minScan_) continue;

            GRBLinExpr expr = SchedulerILP::getX(currT, currSrc, sta) + \
                SchedulerILP::getX(nextT, nextSrc, sta);

            model_->addConstr(expr <= 1);
        }

        // a valid scan requires >= 2 participating stations
        for(unsigned int i = 0; i < blockCount_; ++i) {
            for(const auto src : sourceList_.getSources()) {
                GRBLinExpr expr;
                for(const Station& sta : network_.getStations())
                    expr += SchedulerILP::getX(i, src, sta);

                model_->addConstr(expr <= (SchedulerILP::getY(i, src) * network_.getNSta()));
                model_->addConstr(expr >= SchedulerILP::getY(i, src) * 2);
            }
        }

        // keep sky coverage cells up to date
        for(Station& sta : network_.refStations()) {
            for(std::size_t i = 0; i < SKY_COVERAGE_CELL_COUNT; ++i) {
                GRBLinExpr expr;
                for(unsigned int j = 0; j < blockCount_; ++j) {
                    for(const auto src : sourceList_.getSources()) {
                        if(SchedulerILP::calculateCell(j, src, sta) == i) 
                            expr += SchedulerILP::getX(j, src, sta);
                    }
                }

                model_->addConstr(SchedulerILP::getZ(sta, i) <= expr);
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added constraints";
#else
        std::cout << "[info] Added constraints";
#endif

        min_ = model_->addVar(0.0, SKY_COVERAGE_CELL_COUNT, 0.0, GRB_CONTINUOUS);
        for(const Station& sta : network_.getStations()) {
            GRBLinExpr expr;
            for (std::size_t i = 0; i < SKY_COVERAGE_CELL_COUNT; ++i) expr += SchedulerILP::getZ(sta, i);
            model_->addConstr(min_ <= expr);
        }

        GRBLinExpr obj { min_ };
        model_->setObjective(obj, GRB_MAXIMIZE);

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Configured objective function";
#else
        std::cout << "[info] Configured objective function";
#endif
    }

    const GRBVar& SchedulerILP::getX(
        unsigned int t, 
        const std::shared_ptr<const AbstractSource> src, const Station& sta) const {
        auto idx = t * sourceList_.getNSrc() * network_.getNSta() + \
            src2idx_.at(src->getId()) * network_.getNSta() + sta2idx_.at(sta.getId());
        return x_[idx];
    }

    const GRBVar& SchedulerILP::getY(
        unsigned int t, 
        const std::shared_ptr<const AbstractSource> src) const {
        auto idx = t * sourceList_.getNSrc() + src2idx_.at(src->getId());
        return y_[idx];
    }

    const GRBVar& SchedulerILP::getZ(
        const Station& sta,
        const std::size_t idx) const {
        return z_[sta2idx_.at(sta.getId()) * SKY_COVERAGE_CELL_COUNT + idx];
    }

    const std::size_t SchedulerILP::calculateCell(
        unsigned int t, 
        const std::shared_ptr<const AbstractSource> src,
        Station& sta) const {
        constexpr double el_space = halfpi / 2.;

        PointingVector pv{ sta.getId(), src->getId() };
        pv.setTime(t * minScan_);
        sta.calcAzEl_rigorous(src, pv);

        std::size_t row = static_cast<std::size_t>( floorl( pv.getEl() / el_space ) );
        std::size_t idx;
        switch ( row ) {
            case 0: {
                double n = 9;
                double az_space = twopi / n;
                std::size_t col = static_cast<std::size_t>( roundl( util::wrap2twoPi( pv.getAz() ) / az_space ) );
                if ( static_cast<double>(col) > n - 1 ) col = 0;
                idx = col;
                break;
            }
            default: {
                double n = 4;
                double az_space = twopi / n;
                std::size_t col = static_cast<std::size_t>( roundl( util::wrap2twoPi( pv.getAz() ) / az_space ) );
                if ( static_cast<double>(col) > n - 1 ) col = 0;
                idx = 9 + col;
                break;
            }
        }

        return idx;
    }

    unsigned int SchedulerILP::calculateSlewTime(
        Station& sta, 
        const std::shared_ptr<const AbstractSource> currSrc, 
        const std::shared_ptr<const AbstractSource> nextSrc,
        const double currT,
        const double nextT) const {
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
    
    void SchedulerILP::start() noexcept {
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
                    if(SchedulerILP::getX(i, src, sta).get(GRB_DoubleAttr_X) < 0.5) continue;
                    
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

            std::vector<Scan> scansCurr;
            for(std::vector<PointingVector>& pvSub : pvs) {
                if(pvSub.empty()) continue;
                for(PointingVector& pv : pvSub) {
                    network_.refSkyCoverage(network_.getStaid2skyCoverageId().at(pv.getStaid())).update(pv);
                }

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

                for(unsigned long staid : scanCurr.getStationIds()) {
                    network_.refStation(staid).addObservingTime(minScan_);
                }

                std::vector<Observation> obs;

                std::vector<unsigned long> staids = scanCurr.getStationIds();
                for(std::size_t i = 0; i < staids.size(); ++i) {
                    for(std::size_t j = i + 1; j < staids.size(); ++j) {
                        auto b = network_.getBaseline(std::make_pair(staids[i], staids[j]));
                        obs.emplace_back(b.getId(), staids[i], staids[j], scanCurr.getSourceId(), i * minScan_, minScan_);
                    }
                }

                scanCurr.setObservations(obs);

                sourceList_.refSource(scanCurr.getSourceId())->update(scanCurr.getNSta(), scanCurr.getNObs(), minScan_, true);
                scansCurr.push_back(scanCurr);
            }
            
            if(scansCurr.size() == 1) {
                scans_.push_back(scansCurr[0]);
            } else for(const Scan& s : scansCurr) {
                std::vector<PointingVector> pvSub;
                for(std::size_t i = 0; i < s.getNSta(); ++i) 
                    pvSub.push_back(s.getPointingVector(i));        

                std::vector<PointingVector> pvSubEnd(pvSub);
                for(PointingVector& end : pvSubEnd) {
                    end.setTime(end.getTime() + minScan_);
                    network_.refStation(end.getStaid()).calcAzEl_rigorous(sourceList_.getSource(end.getSrcid()), end);
                }

                Scan scanSub(pvSub, s.getTimes(), s.getObservations());
                scanSub.setPointingVectorsEndtime(pvSubEnd);

                scans_.push_back(scanSub);
            }

#if 0
            // TODO: Statistics do not consider subnetting scans unless this is toggled on
            // However, if subnetting scans are accumulated then the schedule has incorrect tragets
            if(scansCurr.size() == 1) {
                scans_.push_back(scansCurr[0]);
            } else if(!scansCurr.empty()) {
                std::vector<Observation> obsMerged;
                std::vector<PointingVector> pvsMerged;
                for(const Scan& s : scansCurr) {
                    for(const Observation& o : s.getObservations()) obsMerged.push_back(o);
                    // obsMerged.insert(obsMerged.end(), s.getObservations().begin(), s.getObservations().end());
                    for(std::size_t i = 0; i < s.getNSta(); ++i) pvsMerged.push_back(s.getPointingVector(i));                                  
                }

                std::vector<PointingVector> pvsMergedEnd(pvsMerged);
                for(PointingVector& end : pvsMergedEnd) {
                    end.setTime(end.getTime() + minScan_);
                    network_.refStation(end.getStaid()).calcAzEl_rigorous(sourceList_.getSource(end.getSrcid()), end);
                }

                ScanTimes timesMerged(pvsMerged.size());
                timesMerged.setObservingStarts(i * minScan_);
                timesMerged.setObservingTimes(minScan_);

                Scan scanSub(pvsMerged, timesMerged, obsMerged);
                scanSub.setPointingVectorsEndtime(pvsMergedEnd);
                scans_.push_back(scanSub);
            }
#endif

            for(size_t j = 0; j < eolsMask.size(); ++j) if(eolsMask[j]) eols[j] = i * minScan_;
        }
        
        for(Station& sta : network_.refStations()) {
            Station::Statistics stat;
            stat.totalPreobTime = 0;
            stat.totalFieldSystemTime = 0;
            stat.totalSlewTime = 0;
            stat.totalObservingTime = sta.getTotalObservingTime();

            PointingVector prev = sta2pv0_.at(sta.getId());
            for(const Scan& scan : scans_) {
                if(scan.findIdxOfStationId(sta.getId())) {
                    stat.totalFieldSystemTime += sta.getPARA().systemDelay;
                    unsigned int slewTime = SchedulerILP::calculateSlewTime(sta,
                        sourceList_.getSource(prev.getSrcid()), 
                        sourceList_.getSource(scan.getSourceId()), 
                        prev.getTime(),
                        scan.getTimes().getScanTime());
#ifdef VIESCHEDPP_LOG
                    BOOST_LOG_TRIVIAL( info ) << sta.getName() << " is slewing from " 
                        << sourceList_.getSource(prev.getSrcid())->getName() << " to " 
                        << sourceList_.getSource(scan.getSourceId())->getName() << " over " 
                        << slewTime << " seconds";
#else
                    std::cout << sta.getName() << " is slewing from " 
                        << sourceList_.getSource(prev.getSrcid())->getName() << " to " 
                        << sourceList_.getSource(scan.getSourceId())->getName() << " over " 
                        << slewTime << " seconds";
#endif
                    stat.totalSlewTime += slewTime;
                    for(std::size_t i = 0; i < scan.getNSta(); ++i) {
                        if(scan.getPointingVector(i).getStaid() == sta.getId()) {
                            prev = scan.getPointingVector(i);
                            break;
                        }
                    }
                    stat.totalPreobTime += sta.getPARA().preob;
                }
            }

            stat.totalIdleTime = TimeSystem::duration - stat.totalObservingTime - stat.totalSlewTime;
            stat.totalObservingTime -= stat.totalPreobTime + stat.totalFieldSystemTime;
            
            sta.setStatistics(stat);
        }
    }
#else
    void SchedulerILP::initialize() noexcept { /* STUB */ }
    void SchedulerILP::start() noexcept {
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

    Subcon SchedulerILP::createSubcon( const std::shared_ptr<Subnetting> &subnetting, Scan::ScanType type,
                         const boost::optional<StationEndposition> &endposition ) noexcept {
        return Scheduler::createSubcon(subnetting, type, endposition);
    }

    Subcon SchedulerILP::allVisibleScans( Scan::ScanType type, const boost::optional<StationEndposition> &endposition,
                            bool doNotObserveSourcesWithinMinRepeat ) noexcept {
        return Scheduler::allVisibleScans(type, endposition, doNotObserveSourcesWithinMinRepeat);
    }

    void SchedulerILP::update( Scan &scan, std::ofstream &of ) noexcept {
        return Scheduler::update(scan, of);
    }

    void SchedulerILP::consideredUpdate( unsigned long n1scans, unsigned long n2scans, int depth, std::ofstream &of ) noexcept {
        return Scheduler::consideredUpdate(n1scans, n2scans, depth, of);
    }

    void SchedulerILP::statistics( std::ofstream &of ) {
        return Scheduler::statistics(of);
    }

    void SchedulerILP::highImpactScans( HighImpactScanDescriptor &himp, std::ofstream &of ) {
        return Scheduler::highImpactScans(himp, of);
    }

    void SchedulerILP::calibratorBlocks( std::ofstream &of ) {
        return Scheduler::calibratorBlocks(of);
    }

    void SchedulerILP::parallacticAngleBlocks( std::ofstream &of ) {
        return Scheduler::parallacticAngleBlocks(of);
    }

    void SchedulerILP::differentialParallacticAngleBlocks( std::ofstream &of ) {
        return Scheduler::differentialParallacticAngleBlocks(of);
    }

    bool SchedulerILP::checkAndStatistics( std::ofstream &of ) noexcept {
        return Scheduler::checkAndStatistics(of);
    }

    void SchedulerILP::checkSatelliteAvoidance() {
        return Scheduler::checkSatelliteAvoidance();
    }
} // namespace VieVS
