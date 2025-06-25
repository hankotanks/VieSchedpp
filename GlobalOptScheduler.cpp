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
        unsigned int minScan = getMinScan(network_);

        std::cout << "src: " << sourceList_.getNSrc() << std::endl;
        std::cout << "sta: " << network_.getNSta() << std::endl;
        std::cout << "dur: " << TimeSystem::duration / minScan << std::endl;

        for(const auto src : sourceList_.getSources()) {
            src2idx_.insert(std::make_pair(src->getId(), src2idx_.size()));
        }

        for(const auto& sta : network_.getStations()) {
            sta2idx_.insert(std::make_pair(sta.getId(), sta2idx_.size()));
        }

        env_ = new GRBEnv(true);
        env_->start();
        model_ = new GRBModel(*env_);
        
        // add variables
        for(unsigned int i = 0; i < TimeSystem::duration / minScan; ++i) {
            for(const auto src : sourceList_.getSources()) {
                std::string yId("y_");
                yId += std::to_string(i);
                yId += "_";
                yId += src->getName();
                y_.emplace_back(model_->addVar(0.0, 1.0, 0.0, GRB_BINARY, yId));

                for(const Station& sta : network_.getStations()) {
                    std::string xId("x_");
                    xId += std::to_string(i);
                    xId += "_";
                    xId += src->getName();
                    xId += "_";
                    xId += sta.getName();
                    x_.emplace_back(model_->addVar(0.0, 1.0, 0.0, GRB_BINARY, xId));
                }
            }
        }
    }

    const GRBVar& GlobalOptScheduler::getX(const unsigned int t, 
        std::shared_ptr<AbstractSource> src, const Station& sta) {

        unsigned int minScan = getMinScan(network_);

        auto idx = t / (TimeSystem::duration / minScan) * \
            sourceList_.getNSrc() * network_.getNSta() + \
            src2idx_[src->getId()] * network_.getNSta() + sta2idx_[sta.getId()];

        return x_[idx];
    }

    const GRBVar& GlobalOptScheduler::getY(const unsigned int t, std::shared_ptr<AbstractSource> src) {
        unsigned int minScan = getMinScan(network_);

        auto idx = t / (TimeSystem::duration / minScan) * \
            sourceList_.getNSrc() + src2idx_[src->getId()];

        return y_[idx];
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
