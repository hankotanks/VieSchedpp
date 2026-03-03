#include "Model.h"
#include <tuple>
#include "gurobi_c++.h"

#ifdef WITH_GUROBI

#include "Misc/TimeSystem.h"
#include "Scan/PointingVector.h"

#define SKY_COVERAGE_CELL_COUNT 13

namespace {
    void initGurobi(GRBEnv* env, GRBModel* model) {
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
} // private

namespace VieVS {
    Model::Model(VieVS::Network& network, VieVS::SourceList& sourceList, size_t blockLength) : 
        network_(network), sourceList_(sourceList), 
        blockLength_(blockLength), blockCount_(TimeSystem::duration / blockLength - 1) {
        
        initGurobi(env_, model_);

        // StaActive and BlnActive
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const auto q : sourceList_.getSources()) {
                for(const Station& s1 : network_.getStations()) {
                    // make sure source is visible at this time
                    if(!Model::checkStationVisibility(t, q, s1)) continue;
                    // create variable
                    GRBVar var = model_->addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    var_.emplace(ModelKey::StaActive(q, s1, t), var.index());
                    // baseline variables
                    for(const Station& s2 : network_.getStations()) {
                        if(s1.getId() >= s2.getId()) continue;
                        if(!Model::checkStationVisibility(t, q, s2)) continue;
                        // create variable
                        var = model_->addVar(0.0, 1.0, 0.0, GRB_BINARY);
                        var_.emplace(ModelKey::BlnActive(q, s1, s2, t), var.index());
                    }
                }
            }
        }

        // Coverage
        for(const Station& s : network_.getStations()) {
            for(std::size_t c = 0; c < SKY_COVERAGE_CELL_COUNT; ++c) {
                GRBVar var = model_->addVar(0.0, 1.0, 0.0, GRB_BINARY);
                var_.emplace(ModelKey::Coverage(s, c), var.index());
            } 
        }

        // s can only observe one q at time t
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const Station& s : network_.getStations()) {
                GRBLinExpr lhs;
                for(const auto q : sourceList_.getSources()) {
                    if(auto var = Model::getVar(ModelKey::StaActive(q, s, t))) {
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
                    if(auto var = Model::getVar(ModelKey::StaActive(q, s1, t))) {
                        lhs = *var;
                    } else goto next_s;
                    for(const Station& s2 : network_.getStations()) {
                        if(s1.getId() == s2.getId()) continue;
                        if(auto var = Model::getVar(ModelKey::StaActive(q, s2, t))) {
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
                for(const Station& s1 : network_.getStations()) {
                    for(const Station& s2 : network_.getStations()) {
                        if(s1.getId() >= s2.getId()) continue;
                        GRBVar lhs, rhs;
                        if(auto var = Model::getVar(ModelKey::StaActive(q, s1, t))) {
                            rhs = *var;
                        } else goto next_s1;
                        if(auto var = Model::getVar(ModelKey::BlnActive(q, s1, s2, t))) {
                            lhs = *var;
                        } else goto next_s2;
                        model_->addConstr(lhs <= rhs);
                        if(auto var = Model::getVar(ModelKey::StaActive(q, s2, t))) {
                            rhs = *var;
                        } else goto next_s2;
                        model_->addConstr(lhs <= rhs);
next_s2:
                        (void) nullptr;
                    }
next_s1:
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
                            if(auto var = Model::getVar(ModelKey::StaActive(q1, s, t1))) {
                                lhs += *var;
                            } else goto next_t1;
                            if(auto var = Model::getVar(ModelKey::StaActive(q2, s, t2))) {
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
                if(auto var = Model::getVar(ModelKey::Coverage(s, c))) {
                    lhs += *var;
                } else throw std::logic_error("unreachable");
                for(size_t t = 0; t < blockCount_; ++t) {
                    for(const auto q : sourceList_.getSources()) {
                        if(Model::calculateCell(t, q, s) != c) continue;
                        if(auto var = Model::getVar(ModelKey::StaActive(q, s, t))) {
                            rhs += *var;
                        }
                    }
                }
                model_->addConstr(lhs <= rhs);
            }
        }
    }

    bool Model::checkStationVisibility(size_t t, 
        std::shared_ptr<const VieVS::AbstractSource> q, const Station& s) const noexcept  {
        // make sure source is visible at this time
        PointingVector pv0(s.getId(), q->getId());
        PointingVector pvf(s.getId(), q->getId());
        pv0.setTime(t * blockLength_);
        pvf.setTime((t + 1) * blockLength_);
        return s.isVisible(pv0) && s.isVisible(pvf);
    }

    boost::optional<GRBVar> Model::getVar(const ModelKey& key) const noexcept {
        try {
            return model_->getVar(var_.at(key));
        } catch(...) {
            return boost::none;
        }
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
}
#endif // WITH_GUROBI