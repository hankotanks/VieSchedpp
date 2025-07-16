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
                std::string yId = "y_" + std::to_string(i) + "_" + src->getName();
                y_.emplace_back(model_->addVar(0.0, 1.0, 0.0, GRB_BINARY, yId));

                for(const Station& sta : network_.getStations()) {
                    std::string xId = "x_" + std::to_string(i) + "_" + src->getName() + "_" + sta.getName();
                    x_.emplace_back(model_->addVar(0.0, 1.0, 0.0, GRB_BINARY, xId));
                }
            }
        }

        // add constraints

        // each station can only observe one source at a time
        for(unsigned int i = 0; i < blockCount_; ++i) {
            for(const Station& sta : network_.getStations()) {
                GRBLinExpr expr(0.0);
                for(const auto src : sourceList_.getSources()) expr += GlobalOptScheduler::getX(i, src, sta, true);

                std::string cId = "c0_" + std::to_string(i) + "_" + sta.getName();
                model_->addConstr(expr <= 1, cId);
            }
        }

        // must be sufficient time to slew between two targets
        for(const Station& sta : network_.getStations())
        for(const auto currSrc : sourceList_.getSources()) for(const auto nextSrc : sourceList_.getSources())
        for(unsigned int currT = 0, nextT; currT < blockCount_ - 1; ++currT) for(nextT = currT + 1; nextT < blockCount_; ++nextT) {
            unsigned int slewT = GlobalOptScheduler::calculateSlewTime(currT, nextT, sta, currSrc, nextSrc);
            if((nextT - currT) * (minScan_ - 1) >= slewT) continue;

            GRBLinExpr expr = GlobalOptScheduler::getX(currT, currSrc, sta, true) + \
                GlobalOptScheduler::getX(nextT, nextSrc, sta, true);

            std::string cId = "c1_" + std::to_string(currT) + "_" + \
                sta.getName() + "_" + std::to_string(nextT) + "_" + \
                currSrc->getName() + "_" + nextSrc->getName();
            model_->addConstr(expr <= 1, cId);
        }
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
        const unsigned int currT,
        const unsigned int nextT,
        const Station& sta, 
        const std::shared_ptr<const AbstractSource> currSrc, 
        const std::shared_ptr<const AbstractSource> nextSrc) {
        // TODO: Implement function
        return 0;
    }
    
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

        Scheduler::start();
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
