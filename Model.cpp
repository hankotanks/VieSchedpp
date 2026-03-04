#include "Model.h"
#include <stdexcept>
#include <tuple>
#include "gurobi_c++.h"

#include "Misc/TimeSystem.h"
#include "Scan/PointingVector.h"

#define SKY_COVERAGE_CELL_COUNT 13

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

        // StaActive and BlnActive
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const auto q : sourceList_.getSources()) {
                for(Station& s : network_.refStations()) {
                    // make sure source is visible at this time
                    if(!Model::checkStationVisibility(t, q, s)) continue;
                    // create variable
                    Model::addVar(ModelKey::StaActive(this, q, s, t), 0.0, 1.0, 0.0, GRB_BINARY);
                }
            }
        }

        // BlnActive
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const auto q : sourceList_.getSources()) {
                for(const Baseline& b : network_.getBaselines()) {
                    Station& s1 = network_.refStation(b.getStaid1());
                    Station& s2 = network_.refStation(b.getStaid2());
                    if(!Model::checkStationVisibility(t, q, s1)) continue;
                    if(!Model::checkStationVisibility(t, q, s2)) continue;
                    Model::addVar(ModelKey::BlnActive(this, q, b, t), 0.0, 1.0, 0.0, GRB_BINARY);
                }
            }
        }

        // Coverage
        for(const Station& s : network_.getStations()) {
            for(std::size_t c = 0; c < SKY_COVERAGE_CELL_COUNT; ++c) {
                Model::addVar(ModelKey::Coverage(this, s, c), 0.0, 1.0, 0.0, GRB_BINARY);
            } 
        }

        model_->update();

        // s can only observe one q at time t
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const Station& s : network_.getStations()) {
                GRBLinExpr lhs;
                for(const auto q : sourceList_.getSources()) {
                    if(auto var = getVar(ModelKey::StaActive(this, q, s, t))) {
                        lhs += *var;
                    }
                }
                model_->addConstr(lhs <= 1);
            }
        }

        // if s is observing q at t,
        // >= other station must be active for the same q, t
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
                    model_->addConstr(lhs <= rhs);
next_s:
                    (void) nullptr;
                }
            }
        }

        // if <s1, s2> is active at t, both must observe q at t
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
                        } else goto next_b;
                        model_->addConstr(lhs <= rhs);
                        if(auto var = getVar(ModelKey::StaActive(this, q, s2, t))) {
                            rhs = *var;
                        } else goto next_b;
                        model_->addConstr(lhs <= rhs);
next_b:
                        (void) nullptr;
                    }
            }
        }

        // there must be sufficient time in [t1, t2) for s to slew between q1, q2
        for(Station& s : network_.refStations()) {
            for(const auto q1 : sourceList_.getSources()) {
                for(const auto q2 : sourceList_.getSources()) {
                    for(size_t t1 = 0; t1 < blockCount_; ++t1) {
                        for(size_t t2 = t1 + 1; t2 < blockCount_; ++t2) {
                            size_t slew = Model::calculateSlewTime(s, q1, q2, t1, t2);
                            if((t2 - t1) >= slew) continue;
                            GRBLinExpr lhs;
                            if(auto var = getVar(ModelKey::StaActive(this, q1, s, t1))) {
                                lhs += *var;
                            } else goto next_t1;
                            if(auto var = getVar(ModelKey::StaActive(this, q2, s, t2))) {
                                lhs += *var;
                            } else goto next_t2;
                            model_->addConstr(lhs <= 1);
next_t2:
                            (void) nullptr;
                        }
next_t1:
                        (void) nullptr;
                    }
                }
            }
        }

        // c is 'hit' if >= observations occurred over schedule duration
        for(Station& s : network_.refStations()) {
            for(size_t c = 0; c < SKY_COVERAGE_CELL_COUNT; ++c) {
                GRBLinExpr lhs, rhs;
                if(auto var = getVar(ModelKey::Coverage(this, s, c))) {
                    lhs += *var;
                } else throw std::logic_error("unreachable");
                for(size_t t = 0; t < blockCount_; ++t) {
                    for(const auto q : sourceList_.getSources()) {
                        if(Model::calculateCell(t, q, s) != c) continue;
                        if(auto var = getVar(ModelKey::StaActive(this, q, s, t))) {
                            rhs += *var;
                        }
                    }
                }
                model_->addConstr(lhs <= rhs);
            }
        }

        GRBLinExpr obj;

        // coverage objective
        double co = 1.0 / static_cast<double>(SKY_COVERAGE_CELL_COUNT) / static_cast<double>(network_.getNSta());
        for(const Station& s : network_.refStations()) {
            for(size_t c = 0; c < SKY_COVERAGE_CELL_COUNT; ++c) {
                if(auto var = getVar(ModelKey::Coverage(this, s, c))) {
                    obj += *var;
                } else throw std::logic_error("unreachable");
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

    void Model::optimize(std::vector<Scan>& scans) {
        // only explicitly set all variables to 0 if we have a starting point
        if(!scans.empty()) {
            GRBVar* vars = model_->getVars();
            for(size_t i = 0; i < model_->get(GRB_IntAttr_NumVars); ++i) {
                vars[i].set(GRB_DoubleAttr_Start, 0.0);
            }
        }
        // populate starting values from given scans
        for(const Scan& scan : scans) {
            std::shared_ptr<const VieVS::AbstractSource> const q = sourceList_.getSource(scan.getSourceId());
            const ScanTimes& scanTimes = scan.getTimes();
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
                        throw std::logic_error("unreachable");
                    }
                }
            }
            for(unsigned long sId : scan.getStationIds()) {
                const Station& s = network_.getStation(sId);
                size_t t0 = (scanTimes.getObservingTime(sId) + blockLength_ - 1) / blockLength_;
                size_t tf = t0 + scanTimes.getObservingDuration(sId) / blockLength_;
                tf = std::min(tf, blockCount_);
                for(size_t t = t0; t < tf; ++t) {
                    if(auto var = getVar(ModelKey::StaActive(this, q, s, t))) {
                        var->set(GRB_DoubleAttr_Start, 1.0);
                    } else if(t + 1 != tf) {
                        throw std::logic_error("unreachable");
                    }
                }
            }
        }
        for(Station& s : network_.refStations()) {
            for(size_t c = 0; c < SKY_COVERAGE_CELL_COUNT; ++c) {
                for(const auto q : sourceList_.getSources()) {
                    for(size_t t = 0; t < blockCount_; ++t) {
                        if(calculateCell(t, q, s) != c) continue;
                        if(auto var = getVar(ModelKey::StaActive(this, q, s, t))) {
                            if(var->get(GRB_DoubleAttr_Start) > 0.0) {
                                if(auto var = getVar(ModelKey::Coverage(this, s, c))) {
                                    var->set(GRB_DoubleAttr_Start, 1.0);
                                    goto next_c;
                                } else throw std::logic_error("unreachable");
                            }
                        }
                    }
                }
next_c:
                (void) nullptr;
            }
        }
        model_->optimize();
        // TODO: update scans
    }

    void Model::optimize(void) {
        std::vector<Scan> scans;
        Model::optimize(scans);
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
        if(!s.isVisible(pv0)) return false;
        PointingVector pvf(s.getId(), q->getId());
        pvf.setTime((t + 1) * blockLength_);
        s.calcAzEl_rigorous( q, pvf);
        return s.isVisible(pvf);
    }

    std::size_t Model::calculateCell(size_t t, 
        const std::shared_ptr<const AbstractSource> q,
        Station& s) const noexcept {
        constexpr double el_space = halfpi / 2.;

        PointingVector pv{ s.getId(), q->getId() };
        pv.setTime(t * blockLength_);
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
        if(!s.isVisible(tempVec)) return std::numeric_limits<size_t>::max();

        unsigned int slew = s.getAntenna().slewTime(pv1, pv2);
        return (slew + blockLength_ - 1) / blockLength_;
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
        if(!ret.second) throw std::logic_error("unreachable");
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
            case ModelKey::ModelKeyType::coverage:
                if(key.coverage.s != other.key.coverage.s) return key.coverage.s < other.key.coverage.s;
                return key.coverage.c < other.key.coverage.c;
            default: throw std::logic_error("unreachable");
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

    Model::ModelKey Model::ModelKey::Coverage(const Model* model, const Station& s, size_t c) {
        ModelKey key{};
        key.type = ModelKey::ModelKeyType::coverage;
        key.key.coverage.s = model->sta2idx_.at(s.getId());
        key.key.coverage.c = c;
        return key;
    }
#endif // WITH_GUROBI
}
