/*
 *  VieSched++ Very Long Baseline Interferometry (VLBI) Scheduling Software
 *  Copyright (C) 2018  Matthias Schartner
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Model.h"

// system
#include <memory>
#include <numeric>

// gurobi
#ifdef WITH_GUROBI
#include "gurobi_c++.h"
#endif // WITH_GUROBI

namespace VieVS {
void Model::prepare(const Window& window) {
    // build source submask if windowSourceRestriction enabled
    std::set<unsigned long> sourceSubMask;
    if(windowSourceRestriction_) {
        for(size_t t : sol_.getBlocks(window.t0, window.tn)) {
            for(const Station& s : sol_.getStations()) {
                for(const auto q : sol_.getSources(t, s)) {
                    if(*sol_.getSol(Solution::Key::StaActive(&sol_, q, s, t))) {
                        sourceSubMask.insert(q->getId());
                    }
                }
            }
        }
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Restricted window source list to " << sourceSubMask.size() << " due to windowSourceRestriction = true";
#else
        std::cout << "[info] Restricted window source list to " << sourceSubMask.size() << " due to windowSourceRestriction = true";
#endif
    }

    // add constraints
    Model::constrExclusive(window, sourceSubMask);
    Model::constrPairwise(window, sourceSubMask);
    Model::constrBaseline(window, sourceSubMask);
    Model::constrBaselineExclusivity(window, sourceSubMask);
    Model::constrSlew(window, sourceSubMask);
    Model::constrSNR(window, sourceSubMask);
    Model::constrDuration(window, sourceSubMask);
    Model::constrCoverage(window, sourceSubMask);

    model_->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

    // sky coverage objective
    if(objectives_.count("coverage") == 0 || objectives_["coverage"].first) {
        model_->setObjectiveN(Model::objSkyCov(window, sourceSubMask), 0, 2);
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added coverage objective";
#else
        std::cout << "[info] Added coverage objective";
#endif
        if(objectives_["coverage"].second > 0) {
            model_->getMultiobjEnv(0).set(GRB_DoubleParam_TimeLimit, objectives_["coverage"].second);
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Set coverage objective's timeout to " << objectives_["coverage"].second << " seconds";
#else
        std::cout << "[info] Set coverage objective's timeout to " << objectives_["coverage"].second << " seconds";
#endif
        }
    }

    // baseline objective
    if(objectives_.count("baseline") == 0 || objectives_["baseline"].first) {
        model_->setObjectiveN(Model::objBaselines(window, sourceSubMask), 1, 1);
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added baseline objective";
#else
        std::cout << "[info] Added baseline objective";
#endif
        if(objectives_["baseline"].second > 0) {
            model_->getMultiobjEnv(1).set(GRB_DoubleParam_TimeLimit, objectives_["baseline"].second);
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Set baseline objective's timeout to " << objectives_["baseline"].second << " seconds";
#else
        std::cout << "[info] Set baseline objective's timeout to " << objectives_["baseline"].second << " seconds";
#endif
        }
    }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Finished building ILP model";
#else
        std::cout << "[info] Finished building ILP model";
#endif
}

#ifdef WITH_GUROBI
void Model::constrExclusive(const Window& window, const std::set<unsigned long>& sourceSubMask) {
    // s can only observe one q at time t
    size_t count = 0;
    for(size_t t : sol_.getBlocks(window.t0, window.tf)) {
        for(const Station& s : sol_.getStations()) {
            GRBLinExpr lhs;
            size_t lhsCount = 0;
            for(const auto q : sol_.getSources(t, s)) {
                if(!sourceSubMask.empty() && sourceSubMask.count(q->getId()) == 0) continue;
                lhs += *getVar(Solution::Key::StaActive(&sol_, q, s, t));
                lhsCount++;
            }
            if(lhsCount > 0) {
                model_->addConstr(lhs <= 1, "constr_exclusive[" + s.getName() + ", " + std::to_string(t) + "]");
                count++;
            }
        }
    }

#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "Added " << count << " observation exclusivity constraints to model";
#else
    std::cout << "[info] Added " << count << " observation exclusivity constraints to model";
#endif
}

void Model::constrBaseline(const Window& window, const std::set<unsigned long>& sourceSubMask) {
    // if <s1, s2> is active at t, both must observe q at t
    size_t count = 0;
    for(size_t t : sol_.getBlocks(window.t0, window.tf)) {
        for(const auto q : sol_.getSources()) {
            if(!sourceSubMask.empty() && sourceSubMask.count(q->getId()) == 0) continue;
            for(const Baseline& b : sol_.getBaselines(t, q)) {
                auto s = sol_.getStations(b);
                const Station& s1 = s.first;
                const Station& s2 = s.second;
                GRBVar lhs, rhs;
                lhs = *getVar(Solution::Key::BlnActive(&sol_, q, b, t));
                rhs = *getVar(Solution::Key::StaActive(&sol_, q, s1, t));
                model_->addConstr(lhs <= rhs, "constr_baseline[<" + s1.getName() + ", " + s2.getName() + ">, " + q->getName() + ", " + std::to_string(t) + "]");
                rhs = *getVar(Solution::Key::StaActive(&sol_, q, s2, t));
                model_->addConstr(lhs <= rhs, "constr_baseline[<" + s2.getName() + ", " + s1.getName() + ">, " + q->getName() + ", " + std::to_string(t) + "]");
                count += 2;
            }
        }
    }

#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "Added " << count << " baseline constraints to model";
#else
    std::cout << "[info] Added " << count << " baseline constraints to model";
#endif
}

void Model::constrBaselineExclusivity(const Window& window, const std::set<unsigned long>& sourceSubMask) {
    // b can only observe one q at time t
    size_t count = 0;
    for(size_t t : sol_.getBlocks(window.t0, window.tf)) {
        for(const Baseline& b : sol_.getBaselines()) {
            GRBLinExpr lhs;
            size_t lhsCount = 0;
            for(const auto q : sol_.getSources(t, b)) {
                if(!sourceSubMask.empty() && sourceSubMask.count(q->getId()) == 0) continue;
                lhs += *getVar(Solution::Key::BlnActive(&sol_, q, b, t));
                lhsCount++;
            }
            if(lhsCount > 0) {
                model_->addConstr(lhs <= 1, "constr_baseline_exclusive[" + b.getName() + ", " + std::to_string(t) + "]");
                count++;
            }
        }
    }

#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "Added " << count << " baseline observation exclusivity constraints to model";
#else
    std::cout << "[info] Added " << count << " baseline observation exclusivity constraints to model";
#endif
}

void Model::constrPairwise(const Window& window, const std::set<unsigned long>& sourceSubMask) {
    // if s is observing q at t,
    // >= N-1 other station must be active for the same q, t
    size_t count = 0;
    for(size_t t : sol_.getBlocks(window.t0, window.tf)) {
        for(const auto q : sol_.getSources()) {
            if(!sourceSubMask.empty() && sourceSubMask.count(q->getId()) == 0) continue;
            for(const Station& s1 : sol_.getStations(t, q)) {
                GRBVar lhs = *getVar(Solution::Key::StaActive(&sol_, q, s1, t));
                GRBLinExpr rhs;
                size_t rhsCount = 0;
                for(const Station& s2 : sol_.getStations(t, q)) {
                    if(s1.getId() == s2.getId()) continue;
                    rhs += *getVar(Solution::Key::StaActive(&sol_, q, s2, t));
                    rhsCount++;
                }
                if(rhsCount > 0) {
                    model_->addConstr(rhs >= (q->getPARA().minNumberOfSites - 1) * lhs, "constr_pairwise[" + s1.getName() + ", " + q->getName() + ", " + std::to_string(t) + "]");
                    count++;
                }
            }
        }
    }

#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "Added " << count << " pairwise observation constraints to model";
#else
    std::cout << "[info] Added " << count << " pairwise observation constraints to model";
#endif
}

void Model::constrDuration(const Window& window, const std::set<unsigned long>& sourceSubMask) {
    size_t count = 0;
    for(const Station& s : sol_.getStations()) {
        for(const auto q : sol_.getSources()) {
            if(!sourceSubMask.empty() && sourceSubMask.count(q->getId()) == 0) continue;
            for(size_t t2 : sol_.getBlocks(window.t0, window.tn)) {
                size_t maxScan = sol_.getBlocks(std::min(q->getPARA().maxScan, s.getPARA().maxScan));
                if(t2 < maxScan) continue;
                // look backwards by minScan segments and forbid
                GRBLinExpr lhs;
                for(size_t k = 0; k <= maxScan; ++k) {
                    size_t t1 = t2 - k;
                    if(t1 < window.t0) {
                        if(auto result = sol_.getSol(Solution::Key::StaActive(&sol_, q, s, t1))) {
                            if(*result) {
                                --maxScan;
                            } else break;
                        } else break;
                    } else if(auto var = getVar(Solution::Key::StaActive(&sol_, q, s, t1))) {
                        lhs += *var;
                    }
                }
                model_->addConstr(lhs <= maxScan, "constr_duration[" + s.getName() + ", " + q->getName() + " , " + std::to_string(t2) + "]");
                ++count;
            }
        }
    }
#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "Added " << count << " max scan duration constraints to model";
#else
    std::cout << "[info] Added " << count << " max scan duration constraints to model";
#endif
}

void Model::constrSNR(const Window& window, const std::set<unsigned long>& sourceSubMask) {
    size_t count = 0;
    auto add = [this](GRBLinExpr& expr, const Solution::Key& key) {
        if(auto var = this->getVar(key)) {
            expr += *var;
        } else if(auto result = this->sol_.getSol(key)) {
            expr += (*result) ? 1 : 0;
        } else return false;
        return true;
    };

    for(const auto q : sol_.getSources()) {
        if(!sourceSubMask.empty() && sourceSubMask.count(q->getId()) == 0) continue;
        for(const Baseline& b : sol_.getBaselines()) {
            auto s = sol_.getStations(b);
            const Station& s1 = s.first;
            const Station& s2 = s.second;
            for(size_t t1 : sol_.getBlocks(window.tp, window.tn, q, b)) {
                auto key = Solution::Key::BlnActive(&sol_, q, b, t1);
                size_t dur = sol_.getMinObs(key);
                GRBVar lhs = *getVar(key);
                GRBLinExpr rhs;
                if(!add(rhs, Solution::Key::BlnActive(&sol_, q, b, t1 + 1))) continue;
                for(size_t i = 1; i < dur; ++i) {
                    GRBLinExpr inner{rhs};
                    add(inner, Solution::Key::BlnActive(&sol_, q, b, t1 - i));
                    model_->addConstr(lhs <= inner, "constr_snr[<" + s1.getName() + ", " + s2.getName() + ">, " + q->getName() + " , " + std::to_string(t1) + ", " + std::to_string(dur) + ", " + std::to_string(i) + "]");
                    count++;
                }   
            }
        }
    }
#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "Added " << count << " SNR constraints to model";
#else
    std::cout << "[info] Added " << count << " SNR constraints to model";
#endif
}

void Model::constrSlew(const Window& window, const std::set<unsigned long>& sourceSubMask) {
    // there must be sufficient time in [t1, t2) for s to slew between q1, q2
    size_t count = 0;
    for(const Station& s : sol_.getStations()) {
        for(const auto q1 : sol_.getSources()) {
            if(!sourceSubMask.empty() && sourceSubMask.count(q1->getId()) == 0) continue;
            for(size_t t1 : sol_.getBlocks(window.tp, window.tf, q1, s)) {
                auto key_from = Solution::Key::StaActive(&sol_, q1, s, t1);
                auto lhs = *getVar(key_from);
                GRBLinExpr rhs;
                for(const auto q2 : sol_.getSources()) {
                    if(q1->getId() == q2->getId()) continue;
                    if(!sourceSubMask.empty() && sourceSubMask.count(q2->getId()) == 0) continue;
                    for(size_t t2 : sol_.getBlocks(t1, window.tn, q2, s)) {
                        auto key_to = Solution::Key::StaActive(&sol_, q2, s, t2);
                        size_t t_slew = sol_.getSlew(key_from, key_to);
                        if(t1 + t_slew < t2) continue;
                        rhs += *getVar(key_to);
                    }
                }
                model_->addGenConstrIndicator(lhs, 1, rhs == 0, "constr_slew[" + s.getName() + ", " + q1->getName() + " , " + std::to_string(t1) + "]");
                count++;
            }
        }
    }
#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "Added " << count << " slew constraints to model";
#else
    std::cout << "[info] Added " << count << " slew constraints to model";
#endif
}

void Model::constrCoverage(const Window& window, const std::set<unsigned long>& sourceSubMask) {
    // c is 'hit' if >= observations occurred over schedule duration
    size_t count = 0;
    for(const Station& s : sol_.getStations()) {
        for(size_t c = 0; c < coverage_->cellCount(); ++c) {
            GRBLinExpr lhs = *getVar(Solution::Key::StaCoverage(&sol_, s, c));
            GRBLinExpr rhs;
            size_t rhsCount = 0;
            for(size_t t : sol_.getBlocks(window.t0, window.tf)) {
                for(const auto q : sol_.getSources(t, s)) {
                    if(!sourceSubMask.empty() && sourceSubMask.count(q->getId()) == 0) continue;
                    auto pv = sol_.getPointingVector(Solution::Key::StaActive(&sol_, q, s, t));
                    if(coverage_->calculateCell(*pv) != c) continue;
                    rhs += *getVar(Solution::Key::StaActive(&sol_, q, s, t));
                    rhsCount++;
                }
            }
            if(rhsCount > 0) {
                model_->addConstr(lhs <= rhs, "constr_coverage[" + s.getName() + ", " + std::to_string(c) + "]");
                count++;
            }
        }
    }

#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "Added " << count << " sky coverage constraints to model";
#else
    std::cout << "[info] Added " << count << " sky coverage constraints to model";
#endif
}

GRBLinExpr Model::objSkyCov(const Window& window, const std::set<unsigned long>& sourceSubMask) {
    auto sta = sol_.getStations();
    // coverage objective
    GRBLinExpr obj;
    double co = 1.0 / static_cast<double>(coverage_->cellCount()) / static_cast<double>(sta.size());
    for(const Station& s : sta) {
        for(size_t c = 0; c < coverage_->cellCount(); ++c) {
            obj += *getVar(Solution::Key::StaCoverage(&sol_, s, c)) * co;
        }
    }

    return obj;
}

GRBLinExpr Model::objBaselines(const Window& window, const std::set<unsigned long>& sourceSubMask) {
    // baseline occurrence
    std::map<unsigned long, double> bLength;
    for(const Baseline& b : sol_.getBaselines()) {
        auto s = sol_.getStations(b);
        const Station& s1 = s.first;
        const Station& s2 = s.second;
        double length = s1.getPosition()->getDistance(*s2.getPosition());
        bLength.insert(std::make_pair(b.getId(), length));
    }

    auto it = std::max_element(bLength.begin(), bLength.end(),
        [](const auto& l1, const auto& l2) { return l1.second < l2.second; });
    double bLengthMax = it->second;

    std::map<unsigned long, double> bCo;
    std::transform(bLength.begin(), bLength.end(), std::inserter(bCo, bCo.end()), 
        [bLengthMax](const auto& entry) { return std::make_pair(entry.first, std::exp(entry.second / bLengthMax)); });

    double bSum = std::accumulate(bCo.begin(), bCo.end(), 0.0, 
        [](double acc, const auto& entry) { return acc + entry.second; });

    std::for_each(bCo.begin(), bCo.end(), 
        [bSum](auto& entry) { entry.second /= bSum; });

    GRBLinExpr obj;
    for(const Baseline& b : sol_.getBaselines()) {
        double co = bCo.at(b.getId());
        co /= static_cast<double>(window.tf - window.t0);
        for(size_t t : sol_.getBlocks(window.t0, window.tf)) {
            for(const auto q : sol_.getSources(t, b)) {
                if(!sourceSubMask.empty() && sourceSubMask.count(q->getId()) == 0) continue;
                obj += *getVar(Solution::Key::BlnActive(&sol_, q, b, t)) * co;
            }
        }
    }

    return obj;
}
#endif // WITH_GUROBI
}
