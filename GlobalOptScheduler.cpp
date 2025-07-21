#include "GlobalOptScheduler.h"
#include "Scan/PointingVector.h"
#include "Scan/Subcon.h"

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

        for(const auto src : sourceList_.getSources())
            src2idx_.insert(std::make_pair(src->getId(), src2idx_.size()));
        for(const auto& sta : network_.getStations())
            sta2idx_.insert(std::make_pair(sta.getId(), sta2idx_.size()));

        env_ = new GRBEnv(true);
        env_->start();

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
        for(const Station& sta : network_.getStations())
        for(const auto currSrc : sourceList_.getSources()) for(const auto nextSrc : sourceList_.getSources())
        for(unsigned int currT = 0, nextT; currT < blockCount_ - 1; ++currT) for(nextT = currT + 1; nextT < blockCount_; ++nextT) {
            unsigned int slewT = GlobalOptScheduler::calculateSlewTime(sta, currSrc, nextSrc);
            if((nextT - currT) * (minScan_ - 1) >= slewT) continue;

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
            }
        }

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

    unsigned int calculateSlewTime(
        const Station& sta, 
        const std::shared_ptr<const AbstractSource> currSrc, 
        const std::shared_ptr<const AbstractSource> nextSrc) {
        const PointingVector currVec(sta.getId(), currSrc->getId());
        const PointingVector nextVec(sta.getId(), nextSrc->getId());

        return sta.getAntenna().slewTime(currVec, nextVec);
    }
    
    void GlobalOptScheduler::start() noexcept {
        model_->optimize();

        std::vector<unsigned int> eols(network_.getNSta(), 0);

        for(unsigned int i = 0; i < blockCount_; ++i) {
            std::vector<bool> eolsMask(network_.getNSta(), false);
            std::vector<std::vector<PointingVector>> pvs;
            for(const auto src : sourceList_.getSources()) {
                std::vector<PointingVector> pvCurr;
                for(Station& sta : network_.refStations()) {
                    // skip if no observation is made
                    if(GlobalOptScheduler::getX(i, src, sta, true).get(GRB_DoubleAttr_X) < 0.5) continue;
                    
                    // station participates in scan at this timestep
                    PointingVector pv(sta.getId(), src->getId());
                    pv.setTime(i * minScan_);
                    sta.calcAzEl_rigorous(src, pv);
                    pvCurr.emplace_back(sta.getId(), src->getId());

                    eolsMask[sta2idx_[sta.getId()]] = true;
                }

                if(pvCurr.size()) pvs.push_back(std::move(pvCurr));
            }

            if(pvs.size() == 1) {
                // single scan
                std::vector<unsigned int> eolsCurr;
                for(const PointingVector& pv : pvs[0]) eolsCurr.emplace_back(eols[sta2idx_[pv.getStaid()]]);
                
                scans_.emplace_back(pvs[0], eolsCurr);
            } else if(pvs.size()) {
                // subnetting scan
                unsigned int nsta = 0;
                for(const std::vector<PointingVector>& pvSub : pvs) nsta += pvSub.size();
                ScanTimes times(nsta);

                std::vector<unsigned int> eolsCurr;
                for(const std::vector<PointingVector>& pvSub : pvs) 
                    for(const PointingVector& pv : pvSub) 
                        eolsCurr.emplace_back(eols[sta2idx_[pv.getStaid()]]);

                times.setEndOfLastScan(eolsCurr);
                times.setObservingTimes(minScan_);
                times.setObservingStarts(i * minScan_);

                std::vector<PointingVector> pv;
                pv.reserve(nsta);
                for(const std::vector<PointingVector>& pvSub : pvs) 
                    pv.insert(pv.end(), pvSub.begin(), pvSub.end());
                
                std::vector<Observation> obs;
                for(const std::vector<PointingVector>& pvSub : pvs) {
                    for(const PointingVector& pvFst : pvSub) {
                        for(const PointingVector& pvSnd : pvSub) {
                            if(pvFst.getStaid() == pvSnd.getStaid()) continue;
                            obs.emplace_back(
                                network_.getBaseline(std::make_pair(pvFst.getStaid(), pvSnd.getStaid())).getId(),
                                pvFst.getStaid(),
                                pvSnd.getStaid(),
                                pvFst.getSrcid(),
                                pvFst.getTime(),
                                minScan_
                            );
                        }
                    }
                }

                scans_.emplace_back(pv, times, obs);

            }

            for(size_t j = 0; j < eolsMask.size(); ++j) if(eolsMask[j]) eols[j] += minScan_;

        }
    }
#else
    void GlobalOptScheduler::initialize() noexcept { /* STUB */ }
    void GlobalOptScheduler::start() noexcept {
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
