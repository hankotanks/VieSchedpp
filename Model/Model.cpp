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

#include <limits>
#include <memory>
#include <numeric>

#ifdef WITH_GUROBI
#include "gurobi_c++.h"
#endif // WITH_GUROBI

namespace VieVS {
void Model::prepare(size_t tp, size_t t0, size_t tf, size_t tn) {
    Model::constrPairwise(t0, tf);
    Model::constrBaseline(t0, tf);
    Model::constrSlew(tp, t0, tf, tn);
    Model::constrDuration(tp, t0, tf, tn);
    Model::constrCoverage(t0, tf);

#if 0
    // NOTE: should be included in Model::constrSlew now
    Model::constExclusive(t0, tf);
#endif

#if 1
    // TODO: unfinished
    Model::constrSNR(tp, t0, tf, tn);
#endif

    model_->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

    model_->setObjectiveN(Model::objSkyCov(), 0, 2);
    model_->getMultiobjEnv(0).set(GRB_DoubleParam_TimeLimit, 600.0);

    model_->setObjectiveN(Model::objBaselines(t0, tf), 1, 1);
    model_->getMultiobjEnv(1).set(GRB_DoubleParam_TimeLimit, 150.0);

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Finished building ILP model";
#else
        std::cout << "[info] Finished building ILP model";
#endif
}

#ifdef WITH_GUROBI
void Model::constrExclusive(size_t t0, size_t tf) {
    // s can only observe one q at time t
    size_t count = 0;
    for(size_t t : ModelBase::getBlocks(t0, tf)) {
        for(Station& s : ModelBase::getStations()) {
            GRBLinExpr lhs;
            size_t lhsCount = 0;
            for(const auto q : ModelBase::getSources(t, s)) {
                lhs += *getVar(ModelKey::StaActive(this, q, s, t));
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

void Model::constrBaseline(size_t t0, size_t tf) {
    // if <s1, s2> is active at t, both must observe q at t
    size_t count = 0;
    for(size_t t : ModelBase::getBlocks(t0, tf)) {
        for(const auto q : ModelBase::getSources()) {
            for(const Baseline& b : ModelBase::getBaselines(t, q)) {
                const Station& s1 = network_.getStation(b.getStaid1());
                const Station& s2 = network_.getStation(b.getStaid2());
                GRBVar lhs, rhs;
                lhs = *getVar(ModelKey::BlnActive(this, q, b, t));
                rhs = *getVar(ModelKey::StaActive(this, q, s1, t));
                model_->addConstr(lhs <= rhs, "constr_baseline[<" + s1.getName() + ", " + s2.getName() + ">, " + q->getName() + ", " + std::to_string(t) + "]");
                rhs = *getVar(ModelKey::StaActive(this, q, s2, t));
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

void Model::constrPairwise(size_t t0, size_t tf) {
    // if s is observing q at t,
    // >= other station must be active for the same q, t
    size_t count = 0;
    for(size_t t : ModelBase::getBlocks(t0, tf)) {
        for(const auto q : ModelBase::getSources()) {
            for(const Station& s1 : ModelBase::getStations(t, q)) {
                GRBVar lhs = *getVar(ModelKey::StaActive(this, q, s1, t));
                GRBLinExpr rhs;
                size_t rhsCount = 0;
                for(const Station& s2 : ModelBase::getStations(t, q)) {
                    if(s1.getId() == s2.getId()) continue;
                    rhs += *getVar(ModelKey::StaActive(this, q, s2, t));
                    rhsCount++;
                }
                if(rhsCount > 0) {
                    model_->addConstr(rhs >= std::max(static_cast<unsigned int>(2), q->getPARA().minNumberOfSites - 1) * lhs, "constr_pairwise[" + s1.getName() + ", " + q->getName() + ", " + std::to_string(t) + "]");
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

void Model::constrDuration(size_t tp, size_t t0, size_t tf, size_t tn) {
    size_t count = 0;
    for(Station& s : ModelBase::getStations()) {
        for(const auto q : ModelBase::getSources()) {
            for(size_t t2 : ModelBase::getBlocks(tp, tn)) {
                size_t maxScan = (std::min(q->getPARA().maxScan, s.getPARA().maxScan) + blockLength_ - 1) / blockLength_;
                if(t2 < maxScan) continue;
                // look backwards by minScan segments and forbid
                GRBLinExpr lhs;
                for(size_t k = 0; k <= maxScan; ++k) {
                    size_t t1 = t2 - k;
                    if(t1 < tp) {
                        if(auto sol = getSol(ModelKey::StaActive(this, q, s, t1))) {
                            if(*sol) {
                                --maxScan;
                            } else break;
                        } else break;
                    } else if(auto var = getVar(ModelKey::StaActive(this, q, s, t1))) {
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

void Model::constrSNR(size_t tp, size_t t0, size_t tf, size_t tn) {
    size_t count = 0;
    for(const auto q : ModelBase::getSources()) {
        for(const Baseline& b : ModelBase::getBaselines()) {
            for(size_t t1 : ModelBase::getBlocks(tp, tn, q, b)) {
                GRBVar lhs = *getVar(ModelKey::BlnActive(this, q, b, t1));
                size_t dur = std::numeric_limits<size_t>::max();
                for(auto& mode : modes_->getModes()) {
                    dur = std::min(dur, Model::calculateMinObs(t1, q, b, mode));
                }
                // we want to get the # of blocks before and after to extend it, clipped to the complete schedule bounds
                // if the duration would extend before the beginning of the entire schedule
                size_t dur_prior = 0;
                size_t dur_after = 0;
                if(dur > 0) {
                    dur_prior = std::min(dur - 1, t1);
                    if(t1 + 1 >= blockCount_) {
                        dur_after = 0;
                    } else {
                        dur_after = std::min(dur, blockCount_ - t1);
                    }
                }
                if(t1 - dur_prior >= tf || t1 + dur_after <= t0) continue;
                // next, compute all active baselines outside the observation window
                size_t active = 0;
                GRBLinExpr rhs;
                for(size_t t2 : ModelBase::getBlocks(t1 - dur_prior, t1 + dur_after, q, b)) {
                    if(t2 >= t0 && t2 < tf) {
                        rhs += *getVar(ModelKey::BlnActive(this, q, b, t2));
                    } else if(*getSol(ModelKey::BlnActive(this, q, b, t2))) {
                        active++;
                    }
                }
                if(rhs.size() > 0) {
                    auto s1 = network_.getStation(b.getStaid1());
                    auto s2 = network_.getStation(b.getStaid2());
                    model_->addConstr(rhs + active >= lhs * (dur - 1), "constr_snr[<" + s1.getName() + ", " + s2.getName() + ">, " + q->getName() + " , " + std::to_string(t1) + ", " + std::to_string(dur) + "]");
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

void Model::constrSlew(size_t tp, size_t t0, size_t tf, size_t tn) {
    // there must be sufficient time in [t1, t2) for s to slew between q1, q2
    size_t count = 0;
    for(Station& s : ModelBase::getStations()) {
        for(const auto q1 : ModelBase::getSources()) {
            for(size_t t1 : ModelBase::getBlocks(t0, tf, q1, s)) {
                auto lhs = *getVar(ModelKey::StaActive(this, q1, s, t1));
                GRBLinExpr rhs;
                for(const auto q2 : ModelBase::getSources()) {
                    if(q1->getId() == q2->getId()) continue;
                    for(size_t t2 : ModelBase::getBlocks(t1, tf, q2, s)) {
                        size_t t_slew = Model::calculateSlewTime(s, q1, q2, t1, t2);
                        if(t1 + t_slew < t2) continue;
                        rhs += *getVar(ModelKey::StaActive(this, q2, s, t2));
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

    if(t0 > tp) {
        count = 0;
        for(Station& s : ModelBase::getStations()) {
            for(const auto q1 : ModelBase::getSources()) { // starting
                for(const auto q2 : ModelBase::getSources()) { // ending
                    if(q1->getId() == q2->getId()) continue;
                    for(size_t t1 : ModelBase::getBlocks(tp, t0, q1, s)) { // starting
                        if(*getSol(ModelKey::StaActive(this, q1, s, t1))) {
                            // check if any slew windows extend into the active optimization window
                            size_t t_slew = Model::calculateSlewTime(s, q1, q2, t1, t0);
                            if(t1 + t_slew >= t0) {
                                for(size_t t2 : ModelBase::getBlocks(t0, t1 + t_slew + 1, q2, s)) { // ending
                                    auto var = *getVar(ModelKey::StaActive(this, q2, s, t2));
                                    var.set(GRB_DoubleAttr_LB, 0.0);
                                    var.set(GRB_DoubleAttr_UB, 0.0);
                                    ++count;
                                }
                                goto next_backward;
                            }
                        }
                    }
next_backward:;
                }
            }
        }
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Forbade " << count << " potential observations due to backward slew violations";
#else
        std::cout << "[info] Forbade " << count << " potential observations due to backward slew violations";
#endif
    }
    
    if(tf < tn) {
        count = 0;
        for(Station& s : ModelBase::getStations()) {
            for(const auto q1 : ModelBase::getSources()) { // ending
                for(const auto q2 : ModelBase::getSources()) { // starting
                    if(q1->getId() == q2->getId()) continue;
                    for(size_t t1 : ModelBase::getBlocks(tf, tn, q1, s)) { // ending
                        if(*getSol(ModelKey::StaActive(this, q1, s, t1))) {
                            // check if any slew windows extend into the active optimization window
                            size_t t_slew = Model::calculateSlewTime(s, q2, q1, tf, t1);
                            if(t1 - t_slew < tf) {
                                // force these variables to 0
                                for(size_t t2 : ModelBase::getBlocks(t1 - t_slew, tf, q2, s)) {
                                    auto var = *getVar(ModelKey::StaActive(this, q2, s, t2));
                                    var.set(GRB_DoubleAttr_LB, 0.0);
                                    var.set(GRB_DoubleAttr_UB, 0.0);
                                    count++;
                                }
                                goto next_forward;
                            }
                        }
                    }
next_forward:;
                }
            }
        }
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Forbade " << count << " potential observations due to forward slew violations";
#else
        std::cout << "[info] Forbade " << count << " potential observations due to forward slew violations";
#endif
    }
}

void Model::constrCoverage(size_t t0, size_t tf) {
    // c is 'hit' if >= observations occurred over schedule duration
    size_t count = 0;
    for(Station& s : ModelBase::getStations()) {
        for(size_t c = 0; c < coverage_->cellCount(); ++c) {
            GRBLinExpr lhs = *getVar(ModelKey::StaCoverage(this, s, c));
            GRBLinExpr rhs;
            size_t rhsCount = 0;
            for(size_t t : ModelBase::getBlocks(t0, tf)) {
                for(const auto q : ModelBase::getSources(t, s)) {
                    if(coverage_->calculateCell(this, t, q, s) != c) continue;
                    rhs += *getVar(ModelKey::StaActive(this, q, s, t));
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

GRBLinExpr Model::objSkyCov() {
    // coverage objective
    GRBLinExpr obj;
    double co = 1.0 / static_cast<double>(coverage_->cellCount()) / static_cast<double>(network_.getNSta());
    for(const Station& s : ModelBase::getStations()) {
        for(size_t c = 0; c < coverage_->cellCount(); ++c) {
            obj += *getVar(ModelKey::StaCoverage(this, s, c)) * co;
        }
    }

    return obj;
}

GRBLinExpr Model::objBaselines(size_t t0, size_t tf) {
    // baseline occurrence
    std::map<unsigned long, double> bLength;
    for(const Baseline& b : ModelBase::getBaselines()) {
        const Station& s1 = network_.getStation(b.getStaid1());
        const Station& s2 = network_.getStation(b.getStaid2());
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
    for(const Baseline& b : ModelBase::getBaselines()) {
        double co = bCo.at(b.getId());
#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << network_.getStation(b.getStaid1()).getName() << 
        "-" << network_.getStation(b.getStaid2()).getName() << " weighting: " << co;
#else
    std::cout << "[info] " << network_.getStation(b.getStaid1()).getName() << 
        "-" << network_.getStation(b.getStaid2()).getName() << " weighting: " << co;
#endif
        co /= static_cast<double>(tf - t0);
        for(size_t t : ModelBase::getBlocks(t0, tf)) {
            for(const auto q : ModelBase::getSources(t, b)) {
                obj += *getVar(ModelKey::BlnActive(this, q, b, t)) * co;
            }
        }
    }

    return obj;
}
#endif // WITH_GUROBI
}
