#include "Model.h"
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>
#include "gurobi_c.h"

#ifdef WITH_GUROBI
#include "gurobi_c++.h"
#endif // WITH_GUROBI

#include "Misc/TimeSystem.h"
#include "Scan/PointingVector.h"

#define UNREACHABLE std::logic_error((boost::format("unreachable: %d") % __LINE__).str())

namespace {
#ifdef WITH_GUROBI
    void initGurobi(GRBEnv*& env, GRBModel*& model) {
        try {
            env = new GRBEnv(true);
            env->start();
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Started GRB environment";
#else
            std::cout << "[info] Started GRB environment";
#endif
            model = new GRBModel(*env);
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Initialized GRB model";
#else
            std::cout << "[info] Initialized GRB model";
#endif    
        } catch (GRBException& e) {
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( error ) << "Gurobi Exception (" << e.getErrorCode() << "): " << e.getMessage();
#else
            std::cout << "[error] Gurobi Exception (" << e.getErrorCode() << "): " << e.getMessage();
#endif
            throw e;
        }
    }
#endif // WITH_GUROBI
} // private

namespace VieVS {
    Model::Model(VieVS::Network& network, VieVS::SourceList& sourceList, unsigned int blockLength, unsigned int windowLength) : 
        network_(network), sourceList_(sourceList), 
        blockLength_(blockLength), blockCount_(TimeSystem::duration / blockLength - 1),
        windowLength_(windowLength), windowBlockCount_((windowLength + blockLength - 1) / blockLength) {
        coverage_ = std::make_unique<ModelCoverage13>();
#ifdef WITH_GUROBI
        initGurobi(env_, model_);

        // build sta2idx_
        for(const Station& s : network_.getStations()) {
            sta2idx_.insert(std::make_pair(s.getId(), sta2idx_.size()));
        }

        // build bln2idx_
        for(const Baseline& b : network_.getBaselines()) {
            bln2idx_.insert(std::make_pair(b.getId(), bln2idx_.size()));
        }

        // build src2idx_
        for(const auto q : sourceList_.getSources()) {
            src2idx_.insert(std::make_pair(q->getId(), src2idx_.size()));
        }

        // StaActive
        size_t count = 0;
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const auto q : sourceList_.getSources()) {
                for(Station& s : network_.refStations()) {
                    // make sure source is visible at this time
                    if(!Model::checkStationVisibility(t, q, s)) continue;
                    // create variable
                    Model::addVar(ModelKey::StaActive(this, q, s, t), 0.0, 1.0, 0.0, GRB_BINARY);
                    count++;
                }
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " StaActive variables to model";
#else
        std::cout << "[info] Added " << count << " StaActive variables to model";
#endif

        // BlnActive
        count = 0;
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const auto q : sourceList_.getSources()) {
                for(const Baseline& b : network_.getBaselines()) {
                    Station& s1 = network_.refStation(b.getStaid1());
                    Station& s2 = network_.refStation(b.getStaid2());
                    if(!Model::checkStationVisibility(t, q, s1)) continue;
                    if(!Model::checkStationVisibility(t, q, s2)) continue;
                    Model::addVar(ModelKey::BlnActive(this, q, b, t), 0.0, 1.0, 0.0, GRB_BINARY);
                    count++;
                }
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " BlnActive variables to model";
#else
        std::cout << "[info] Added " << count << " BlnActive variables to model";
#endif

        // StaCoverage
        count = 0;
        for(const Station& s : network_.getStations()) {
            for(std::size_t c = 0; c < coverage_->cellCount(); ++c) {
                Model::addVar(ModelKey::StaCoverage(this, s, c), 0.0, 1.0, 0.0, GRB_BINARY);
                count++;
            } 
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " StaCoverage variables to model";
#else
        std::cout << "[info] Added " << count << " StaCoverage variables to model";
#endif

        // update the model to make sure variables are accessible
        model_->update();

        // s can only observe one q at time t
        count = 0;
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const Station& s : network_.getStations()) {
                GRBLinExpr lhs;
                for(const auto q : sourceList_.getSources()) {
                    if(auto var = getVar(ModelKey::StaActive(this, q, s, t))) {
                        lhs += *var;
                    }
                }
                model_->addConstr(lhs <= 1, "c0_exclusive");
                count++;
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " observation exclusivity constraints to model";
#else
        std::cout << "[info] Added " << count << " observation exclusivity constraints to model";
#endif

        // if s is observing q at t,
        // >= other station must be active for the same q, t
        count = 0;
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const auto q : sourceList_.getSources()) {
                for(const Station& s1 : network_.getStations()) {
                    GRBVar lhs;
                    GRBLinExpr rhs;
                    if(auto var = getVar(ModelKey::StaActive(this, q, s1, t))) {
                        lhs = *var;
                    } else goto next_s;
                    for(const Station& s2 : network_.getStations()) {
                        if(s1.getId() == s2.getId()) continue;
                        if(auto var = getVar(ModelKey::StaActive(this, q, s2, t))) {
                            rhs += *var;
                        }
                    }
                    model_->addConstr(lhs <= rhs, "c1_pairwise");
                    count++;
next_s:
                    (void) nullptr;
                }
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " pairwise observation constraints to model";
#else
        std::cout << "[info] Added " << count << " pairwise observation constraints to model";
#endif

        // if <s1, s2> is active at t, both must observe q at t
        count = 0;
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const auto q : sourceList_.getSources()) {
                for(const Baseline& b : network_.getBaselines()) {
                        const Station& s1 = network_.getStation(b.getStaid1());
                        const Station& s2 = network_.getStation(b.getStaid2());
                        GRBVar lhs, rhs;
                        if(auto var = getVar(ModelKey::BlnActive(this, q, b, t))) {
                            lhs = *var;
                        } else goto next_b;
                        if(auto var = getVar(ModelKey::StaActive(this, q, s1, t))) {
                            rhs = *var;
                        } else throw UNREACHABLE;
                        model_->addConstr(lhs <= rhs, "c2_baseline");
                        if(auto var = getVar(ModelKey::StaActive(this, q, s2, t))) {
                            rhs = *var;
                        } else throw UNREACHABLE;
                        model_->addConstr(lhs <= rhs, "c2_baseline");
                        count += 2;
next_b:
                        (void) nullptr;
                    }
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " baseline constraints to model";
#else
        std::cout << "[info] Added " << count << " baseline constraints to model";
#endif

        // there must be sufficient time in [t1, t2) for s to slew between q1, q2
        count = 0;
        for(Station& s : network_.refStations()) {
            for(const auto q1 : sourceList_.getSources()) {
                for(const auto q2 : sourceList_.getSources()) {
                    if(q1->getId() == q2->getId()) continue;
                    for(size_t t1 = 0; t1 < blockCount_ - 1; ++t1) {
                        for(size_t t2 = t1 + 1; t2 < blockCount_; ++t2) {
                            size_t slew = Model::calculateSlewTime(s, q1, q2, t1, t2);
                            if(t2 - t1 > slew) continue;
                            GRBLinExpr lhs;
                            if(auto var = getVar(ModelKey::StaActive(this, q1, s, t1))) {
                                lhs += *var;
                            } else goto next_t1;
                            if(auto var = getVar(ModelKey::StaActive(this, q2, s, t2))) {
                                lhs += *var;
                            } else goto next_t2;
                            model_->addConstr(lhs <= 1, "c3_slew");
                            count++;
next_t2:
                            (void) nullptr;
                        }
next_t1:
                        (void) nullptr;
                    }
                }
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " slew constraints to model";
#else
        std::cout << "[info] Added " << count << " slew constraints to model";
#endif

        // c is 'hit' if >= observations occurred over schedule duration
        count = 0;
        for(Station& s : network_.refStations()) {
            for(size_t c = 0; c < coverage_->cellCount(); ++c) {
                GRBLinExpr lhs, rhs;
                if(auto var = getVar(ModelKey::StaCoverage(this, s, c))) {
                    lhs += *var;
                } else throw UNREACHABLE;
                for(size_t t = 0; t < blockCount_; ++t) {
                    for(const auto q : sourceList_.getSources()) {
                        if(coverage_->calculateCell(this, t, q, s) != c) continue;
                        if(auto var = getVar(ModelKey::StaActive(this, q, s, t))) {
                            rhs += *var;
                        }
                    }
                }
                model_->addConstr(lhs <= rhs, "c4_coverage");
                count++;
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " sky coverage constraints to model";
#else
        std::cout << "[info] Added " << count << " sky coverage constraints to model";
#endif

        GRBLinExpr obj;

        // coverage objective
        double co = 1.0 / static_cast<double>(coverage_->cellCount()) / static_cast<double>(network_.getNSta());
        for(const Station& s : network_.refStations()) {
            for(size_t c = 0; c < coverage_->cellCount(); ++c) {
                if(auto var = getVar(ModelKey::StaCoverage(this, s, c))) {
                    obj += (*var) * co;
                } else throw UNREACHABLE;
            }
        }
        model_->setObjective(obj, GRB_MAXIMIZE);

#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Finished building ILP model";
#else
            std::cout << "[info] Finished building ILP model";
#endif
#endif // WITH_GUROBI
    }

    bool Model::optimize(void) {
        model_->update();
        model_->optimize();

        if(model_->get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "No optimal solution found";
#else
            std::cout << "[info] No optimal solution found";
#endif
            return false;
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Completed optimization";
#else
        std::cout << "[info] Completed optimization";
#endif
        return true;
    }

    std::vector<Scan> Model::optimize(std::vector<Scan>& scans) {
        Model::loadScans(scans);

        // std::cout << "Initial:" << std::endl;
        // Model::dump(GRB_DoubleAttr_Start);

        if(!Model::optimize()) return {};

        // std::cout << "Optimized" << std::endl;
        // Model::dump(GRB_DoubleAttr_X);

        return Model::readScans();
    }
}

// ModelCoverage13 implementation
namespace VieVS {
    std::size_t ModelCoverage13::cellCount(void) const noexcept { 
        return 13; 
    }

    std::size_t ModelCoverage13::calculateCell(const Model* model, size_t t, 
        const std::shared_ptr<const AbstractSource> q,
        Station& s) const noexcept {
        constexpr double el_space = halfpi / 2.;

        PointingVector pv{ s.getId(), q->getId() };
        pv.setTime(t * model->blockLength_);
        s.calcAzEl_rigorous(q, pv);

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
}

// helper implementations
namespace VieVS {
    bool Model::checkStationVisibility(size_t t, 
        std::shared_ptr<const VieVS::AbstractSource> q, Station& s) const noexcept  {
        // make sure source is visible at this time
        PointingVector pv0(s.getId(), q->getId());
        pv0.setTime(t * blockLength_);
        s.calcAzEl_rigorous( q, pv0);
        if(!s.isVisible(pv0, q->getPARA().minElevation)) return false;
        PointingVector pvf(s.getId(), q->getId());
        pvf.setTime((t + 1) * blockLength_);
        s.calcAzEl_rigorous( q, pvf);
        return s.isVisible(pvf, q->getPARA().minElevation);
    }

    size_t Model::calculateSlewTime(Station& s, 
        const std::shared_ptr<const AbstractSource> q1, 
        const std::shared_ptr<const AbstractSource> q2,
        size_t t1, size_t t2) const noexcept {
        if(q1->getId() == q2->getId()) return 0;
        
        PointingVector pv1(s.getId(), q1->getId());
        PointingVector pv2(s.getId(), q2->getId());

        pv1.setTime(t1 * blockLength_);
        pv2.setTime(t2 * blockLength_);
        s.calcAzEl_rigorous(q1, pv1);
        s.calcAzEl_rigorous(q2, pv2);

        PointingVector tempVec(pv2);
        if(!s.isVisible(tempVec, q2->getPARA().minElevation)) 
            return std::numeric_limits<size_t>::max();

        unsigned int t_slew = s.getAntenna().slewTime(pv1, pv2);
        unsigned int t_const = s.getPARA().systemDelay + s.getPARA().preob;

        return (t_slew + t_const + blockLength_ - 1) / blockLength_;
    }

#ifdef WITH_GUROBI
    boost::optional<GRBVar> Model::getVar(const ModelKey& key) const noexcept {
        try {
            return var_.at(key);
        } catch(...) {
            return boost::none;
        }
    }

    GRBVar& Model::addVar(const ModelKey& key, double lb, double ub, double obj, char vtype) {
        auto ret = var_.insert(std::make_pair(key, model_->addVar(lb, ub, obj, vtype)));
        if(!ret.second) throw UNREACHABLE;
        return ret.first->second;
    }

    bool Model::ModelKey::operator<(const Model::ModelKey& other) const {
        if(type != other.type) return type < other.type;
        switch(type) {
            case ModelKey::ModelKeyType::sta_active:
                if(key.sta_active.q != other.key.sta_active.q) return key.sta_active.q < other.key.sta_active.q;
                if(key.sta_active.s != other.key.sta_active.s) return key.sta_active.s < other.key.sta_active.s;
                return key.sta_active.t < other.key.sta_active.t;
            case ModelKey::ModelKeyType::bln_active:
                if(key.bln_active.q != other.key.bln_active.q) return key.bln_active.q < other.key.bln_active.q;
                if(key.bln_active.b != other.key.bln_active.b) return key.bln_active.b < other.key.bln_active.b;
                return key.bln_active.t < other.key.bln_active.t;
            case ModelKey::ModelKeyType::sta_coverage:
                if(key.sta_coverage.s != other.key.sta_coverage.s) return key.sta_coverage.s < other.key.sta_coverage.s;
                return key.sta_coverage.c < other.key.sta_coverage.c;
            default: throw UNREACHABLE;
        }
    }

    Model::ModelKey Model::ModelKey::StaActive(const Model* model, 
        std::shared_ptr<const VieVS::AbstractSource> const q, 
        const Station& s, size_t t) {
        ModelKey key{};
        key.type = ModelKey::ModelKeyType::sta_active;
        key.key.sta_active.q = model->src2idx_.at(q->getId());
        key.key.sta_active.s = model->sta2idx_.at(s.getId());
        key.key.sta_active.t = t;
        return key;
    }

    Model::ModelKey Model::ModelKey::BlnActive(const Model* model, 
        std::shared_ptr<const VieVS::AbstractSource> const q, 
        const Baseline& b, size_t t) {
        ModelKey key{};
        key.type = ModelKey::ModelKeyType::bln_active;
        key.key.bln_active.q = model->src2idx_.at(q->getId());
        key.key.bln_active.b = model->bln2idx_.at(b.getId());
        key.key.bln_active.t = t;
        return key;
    }

    Model::ModelKey Model::ModelKey::StaCoverage(const Model* model, const Station& s, size_t c) {
        ModelKey key{};
        key.type = ModelKey::ModelKeyType::sta_coverage;
        key.key.sta_coverage.s = model->sta2idx_.at(s.getId());
        key.key.sta_coverage.c = c;
        return key;
    }

    void Model::loadScans(const std::vector<Scan>& scans) {
        if(scans.empty()) {
#ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "Skipped loading MIP solution. No scans provided";
#else
            std::cout << "[info] Skipped loading MIP solution. No scans provided";
#endif 
            return;
        }

        // only explicitly set all variables to 0 if we have a starting point
        GRBVar* vars = model_->getVars();
        for(size_t i = 0; i < model_->get(GRB_IntAttr_NumVars); ++i) {
            vars[i].set(GRB_DoubleAttr_Start, 0.0);
        }

        // populate starting values from given scans
        for(const Scan& scan : scans) {
            std::shared_ptr<const VieVS::AbstractSource> const q = sourceList_.getSource(scan.getSourceId());
            const ScanTimes& scanTimes = scan.getTimes();
            
            // populate BlnActive variables
            for(const Observation& obs : scan.getObservations()) {
                const Baseline& b = network_.getBaseline(obs.getBlid());
                // observation start blocks
                size_t t10 = (scanTimes.getObservingTime(b.getStaid1()) + blockLength_ - 1) / blockLength_;
                size_t t20 = (scanTimes.getObservingTime(b.getStaid2()) + blockLength_ - 1) / blockLength_;
                // the number of blocks each station is observing
                size_t t1f = t10 + scanTimes.getObservingDuration(b.getStaid1()) / blockLength_;
                size_t t2f = t10 + scanTimes.getObservingDuration(b.getStaid2()) / blockLength_;
                t1f = std::min(t1f, blockCount_);
                t2f = std::min(t2f, blockCount_);
                
                for(size_t t = std::max(t10, t20); t < std::min(t1f, t2f); ++t) {
                    if(auto var = getVar(ModelKey::BlnActive(this, q, b, t))) {
                        var->set(GRB_DoubleAttr_Start, 1.0);
                    } else if(t + 1 != std::min(t1f, t2f)) {
                        throw UNREACHABLE;
                    }
                }
            }

            // populate StaActive variables
            for(unsigned long sId : scan.getStationIds()) {
                Station& s = network_.refStation(sId);
                size_t t0 = (scanTimes.getObservingTime(sId) + blockLength_ - 1) / blockLength_;
                size_t tf = t0 + scanTimes.getObservingDuration(sId) / blockLength_;
                tf = std::min(tf, blockCount_);
                for(size_t t = t0; t < tf; ++t) {
                    if(auto var = getVar(ModelKey::StaActive(this, q, s, t))) {
                        var->set(GRB_DoubleAttr_Start, 1.0);
                    } else if(t + 1 != tf) {
#ifdef VIESCHEDPP_LOG
                        BOOST_LOG_TRIVIAL( warning ) << q->getName() << "is not visible by " << s.getName() << 
                            " for the last " << (tf - 1) << " out of " << tf << " time segments";
#else
                        std::cout << "[warning] " << q->getName() << "is not visible by " << s.getName() << 
                            " for the last " << (tf - 1) << " out of " << tf << " time segments";
#endif
                        break;
                        // TODO: Consider if this is actually acceptable
                        // throw UNREACHABLE;
                    }
                }
            }
        }

        model_->update();

        // populate StaConverage variables
        for(Station& s : network_.refStations()) {
            for(size_t c = 0; c < coverage_->cellCount(); ++c) {
                for(const auto q : sourceList_.getSources()) {
                    for(size_t t = 0; t < blockCount_; ++t) {
                        if(coverage_->calculateCell(this, t, q, s) != c) continue;
                        if(auto var = getVar(ModelKey::StaActive(this, q, s, t))) {
                            if(var->get(GRB_DoubleAttr_Start) > 0.0) {
                                if(auto var = getVar(ModelKey::StaCoverage(this, s, c))) {
                                    var->set(GRB_DoubleAttr_Start, 1.0);
                                    goto next_c;
                                } else throw UNREACHABLE;
                            }
                        }
                    }
                }
next_c:
                (void) nullptr;
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Loaded preliminary result into ILP model";
#else
        std::cout << "[info] Loaded preliminary result into ILP model";
#endif
    }

#if 0
    void Model::activeScanAppend(std::vector<Scan>& scansActive, 
        std::shared_ptr<const VieVS::AbstractSource> const q, 
        const Station& s, size_t t) const noexcept {
        bool found = false;
        for(Scan& scan : scansActive) {
            std::vector<unsigned long> sIds = scan.getStationIds();
            if(std::find(sIds.begin(), sIds.end(), s.getId()) == sIds.end()) continue;
            Scan()
        }
    }   
#endif

    bool Model::ScanBuilder::append(const Model* model, std::shared_ptr<const VieVS::AbstractSource> const q, 
        Station& s, size_t t) noexcept {
        if(qId != q->getId()) return false;
        try {
            size_t& eols = sData.at(s.getId()).second;
            if(eols == t) {
                ++eols;
                return true;
            }
            return false;
        } catch(...) {
            PointingVector pv(s.getId(), qId);
            pv.setTime(t * model->blockLength_);
            s.calcAzEl_rigorous(q, pv);
            sData.insert(std::make_pair(s.getId(), std::make_pair(pv, t + 1)));
            return true;
        }
    }

    Scan Model::ScanBuilder::finish(const Model* model, const std::vector<unsigned int>& slewTime) const noexcept {
        std::vector<PointingVector> pointingVectors;
        std::vector<PointingVector> pointingVectorsEnd;
        std::vector<unsigned int> endOfLastScan;

        typedef std::pair<const unsigned long, std::pair<PointingVector, size_t>> Entry;
        std::transform(sData.begin(), sData.end(), std::back_inserter(pointingVectors),
            [](const Entry& entry) { return entry.second.first; });

        size_t blockLength = model->blockLength_;
        std::transform(sData.begin(), sData.end(), std::back_inserter(endOfLastScan),
            [blockLength](const Entry& entry) { return entry.second.second * blockLength; });

        pointingVectorsEnd.reserve(pointingVectors.size());
        for(size_t i = 0; i < pointingVectors.size(); ++i) {
            const PointingVector& pv0 = pointingVectors[i];
            Station& s = model->network_.refStation(pv0.getStaid());
            std::shared_ptr<const VieVS::AbstractSource> const q = model->sourceList_.getSource(qId);
            
            PointingVector pve(pv0.getStaid(), qId);
            pve.setTime(endOfLastScan[i]);
            s.calcAzEl_rigorous(q, pve);
            pointingVectorsEnd.push_back(pve);
        }

        unsigned int scanStart = std::numeric_limits<unsigned int>::max();

        std::vector<unsigned int> fieldSystemTime;
        std::vector<unsigned int> preob;
        std::vector<unsigned int> observingTimes;

        for(size_t i = 0; i < pointingVectors.size(); ++i) {
            const PointingVector& pv0 = pointingVectors[i];
            const PointingVector& pve = pointingVectorsEnd[i];
            const Station& s = model->network_.getStation(pv0.getStaid());

            scanStart = std::min(scanStart, pv0.getTime());

            fieldSystemTime.push_back(s.getPARA().systemDelay);
            preob.push_back(s.getPARA().preob);
            observingTimes.push_back(pve.getTime() - pv0.getTime());
        }

        Scan scan(pointingVectors, endOfLastScan, Scan::ScanType::standard);
        scan.setPointingVectorsEndtime(pointingVectorsEnd);
        scan.setScanTimes(endOfLastScan, fieldSystemTime, slewTime, preob, scanStart, observingTimes);

        return scan;
    }

    void Model::ActiveScans::append(std::shared_ptr<const VieVS::AbstractSource> const q, 
        Station& s, size_t t) noexcept {
        for(Model::ScanBuilder& scan : scans_) {
            if(scan.append(model_, q, s, t)) return;
        }

        PointingVector pv(s.getId(), q->getId());
        pv.setTime(t * model_->blockLength_);
        s.calcAzEl_rigorous(q, pv);

        Model::ScanBuilder scan;
        scan.qId = q->getId();
        scan.sData.insert(std::make_pair(s.getId(), std::make_pair(pv, t + 1)));

        scans_.emplace_back(scan);
    }

    void Model::ActiveScans::updateScans(std::vector<Scan>& scans, size_t t) {
        std::vector<size_t> remove;
        for(size_t i = 0, j = 0; i < scans_.size(); ++i, j = 0) {
            for(const auto& sEntry : scans_[i].sData) {
                const auto& sData = sEntry.second;
                if(sData.second > t) ++j;
            }

            if(j < 2) {
                std::vector<unsigned int> scanTime;
                for(const auto& sEntryInner : scans_[i].sData) {
                    const Station& s = model_->network_.getStation(sEntryInner.first);
                    if(scans.empty()) {
                        scanTime.push_back(0);
                    } else {
                        boost::optional<PointingVector> pv0 = boost::none;
                        for(size_t i = scans.size(); i-- > 0;) {
                            if(auto j = scans[i].findIdxOfStationId(sEntryInner.first)) {
                                pv0 = scans[i].getPointingVector(*j, Timestamp::end);
                                break;
                            }
                        }
                        unsigned int slew = 0;
                        if(pv0) {
                            slew = s.getAntenna().slewTime(*pv0, sEntryInner.second.first);
                        }
                        scanTime.push_back(slew);
                    }
                }
                Scan scan = scans_[i].finish(model_, scanTime);
                scans.emplace_back(scan);
                remove.emplace_back(i);
            }
        }

        std::reverse(remove.begin(), remove.end());
        for(size_t i : remove) {
            scans_.erase(scans_.begin() + i);
        }
    }

    std::vector<Scan> Model::readScans(void) const noexcept {
        std::vector<Scan> scans;

        Model::ActiveScans scansActive(this);
        for(size_t t = 0; t <= blockCount_; ++t) {
            for(Station& s : network_.refStations()) {
                for(const auto q : sourceList_.getSources()) {
                    if(auto var = getVar(ModelKey::StaActive(this, q, s, t))) {
                        if(var->get(GRB_DoubleAttr_X) < 0.5) continue;
                        scansActive.append(q, s, t);
                    }
                }
            }
            scansActive.updateScans(scans, t);
        }
        return scans;
    }

    void Model::dump(GRB_DoubleAttr attr) const noexcept {
        for(const Station& s : network_.getStations()) {
            std::cout << s.getName() << std::endl;
            for(size_t t = 0; t < blockCount_; ++t) {
                for(const auto q : sourceList_.getSources()) {
                    if(auto var = getVar(ModelKey::StaActive(this, q, s, t))) {
                        if(var->get(attr) > 0.0) {
                            std::cout << 'X';
                            goto next_t;
                        }
                    }
                }
                std::cout << " ";
next_t:
                (void) nullptr;
            }
            std::cout << std::endl;
        }
    }
#endif // WITH_GUROBI
}
