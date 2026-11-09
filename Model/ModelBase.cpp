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

#include "ModelBase.h"

// system
#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

// gurobi
#ifdef WITH_GUROBI
#include "gurobi_c++.h"
#endif // WITH_GUROBI

// VieSchedpp
#include "../Misc/TimeSystem.h"
#include "../Scan/PointingVector.h"

#define UNREACHABLE std::logic_error((boost::format("unreachable: %d") % __LINE__).str())

namespace {
#ifdef WITH_GUROBI
void initGurobi(std::unique_ptr<GRBEnv>& env) {
    try {
        env = std::make_unique<GRBEnv>(true);
        env->start();
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
        throw e;
    }
}

void initModel(std::unique_ptr<GRBEnv>& env, std::unique_ptr<GRBModel>& model) {
    try {
        model = std::make_unique<GRBModel>(*env);
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
ModelBase::ModelBase(VieVS::Network& network, VieVS::SourceList& sourceList, 
    const std::set<unsigned long>& sourceMask, 
    const std::shared_ptr<const ObservingMode>& modes,
    unsigned int blockLength, unsigned int windowLength) :
    sol_(network, sourceList, sourceMask, modes, blockLength, windowLength) {
    coverage_ = std::make_unique<ModelCoverage13>();
#ifdef WITH_GUROBI
    initGurobi(env_);
#endif // WITH_GUROBI
}

bool ModelBase::optimize(void) {
#ifdef WITH_GUROBI
    for (size_t i = 0; i < sol_.windowCount_; ++i) {
        size_t t0 = i * sol_.windowBlockCount_;
        size_t tp = (i == 0) ? 0 : (t0 - sol_.windowBlockCount_);
        size_t tf = std::min(t0 + sol_.windowBlockCount_, sol_.blockCount_);
        size_t tn = std::min(tf + sol_.windowBlockCount_, sol_.blockCount_);

        // reinitialize the model
        initModel(env_, model_);

        // clear the variable map
        var_.clear();

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Complete window spans " << tp * sol_.blockLength_ << " to " << tn * sol_.blockLength_;
#else
        std::cout << "[info] Optimizing window spans " << tp * blockLength_ << " to " << tn * blockLength_;
#endif

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Optimizing between " << t0 * sol_.blockLength_ << " and " << tf * sol_.blockLength_;
#else
        std::cout << "[info] Optimizing between " << t0 * blockLength_ << " and " << tf * blockLength_;
#endif

        // StaActive
        size_t count = 0;
        for(size_t t : sol_.getBlocks(tp, tn)) {
            for(const auto q : sol_.getSources()) {
                for(const Station& s : sol_.getStations(t, q)) {
                    // create variable
                    Solution::Key key = Solution::Key::StaActive(&sol_, q, s, t);
                    GRBVar& var = ModelBase::addVar(key, 0.0, 1.0, 0.0, GRB_BINARY);
                    var.set(GRB_DoubleAttr_Start, 0.0);
                    auto sol = sol_.getSol(key);
                    if(*sol) var.set(GRB_DoubleAttr_Start, 1.0);
                    count++;
                }
            }
        }
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " station activity variables to model";
#else
        std::cout << "[info] Added " << count << " StaActive variables to model";
#endif

        // BlnActive
        count = 0;
        for(size_t t : sol_.getBlocks(tp, tn)) {
            for(const auto q : sol_.getSources()) {
                for(const Baseline& b : sol_.getBaselines(t, q)) {
                    Solution::Key key = Solution::Key::BlnActive(&sol_, q, b, t);
                    GRBVar& var = ModelBase::addVar(key, 0.0, 1.0, 0.0, GRB_BINARY);
                    var.set(GRB_DoubleAttr_Start, 0.0);
                    auto sol = sol_.getSol(key);
                    if(*sol) var.set(GRB_DoubleAttr_Start, 1.0);
                    count++;
                }
            }
        }
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " baseline activity variables to model";
#else
        std::cout << "[info] Added " << count << " BlnActive variables to model";
#endif

        if(t0 > tp && false) {
            for(const Station& s : sol_.getStations()) {
                for(const auto q1 : sol_.getSources()) { // starting
                    for(const auto q2 : sol_.getSources()) { // ending
                        if(q1->getId() == q2->getId()) continue;
                        auto key_to = Solution::Key::StaActive(&sol_, q2, s, t0);
                        for(size_t t1 : sol_.getBlocks(tp, t0, q1, s)) { // starting
                            auto key_from = Solution::Key::StaActive(&sol_, q1, s, t1 + 1);
                            if(*sol_.getSol(key_from)) {
                                // check if any slew windows extend into the active optimization window
                                size_t t_slew = sol_.getSlew(key_from, key_to);
                                if(t1 + t_slew > t0) {
                                    for(size_t t2 : sol_.getBlocks(t0, t1 + t_slew + 1, q2, s)) { // ending
                                        Solution::Key key = Solution::Key::StaActive(&sol_, q2, s, t2);
                                        if(auto sol = sol_.getSol(key)) *sol = false;
                                        auto var = *getVar(key);
                                        var.set(GRB_DoubleAttr_LB, 0.0);
                                        var.set(GRB_DoubleAttr_UB, 0.0);
                                    }
                                    goto next_backward;
                                }
                            }
                        }
next_backward:;
                    }
                }
            }
        }
        
        if(tf < tn && false) {
            for(const Station& s : sol_.getStations()) {
                for(const auto q1 : sol_.getSources()) { // ending
                    for(const auto q2 : sol_.getSources()) { // starting
                        if(q1->getId() == q2->getId()) continue;
                        auto key_from = Solution::Key::StaActive(&sol_, q2, s, tf + 1);
                        for(size_t t1 : sol_.getBlocks(tf, tn, q1, s)) { // ending
                            auto key_to = Solution::Key::StaActive(&sol_, q1, s, t1);
                            if(*sol_.getSol(key_to)) {
                                // check if any slew windows extend into the active optimization window
                                size_t t_slew = sol_.getSlew(key_from, key_to);
                                if(t1 - t_slew < tf) {
                                    // force these variables to 0
                                    for(size_t t2 : sol_.getBlocks(t1 - t_slew, tf, q2, s)) {
                                        Solution::Key key = Solution::Key::StaActive(&sol_, q2, s, t2);
                                        if(auto sol = sol_.getSol(key)) *sol = false;
                                        auto var = *getVar(key);
                                        var.set(GRB_DoubleAttr_LB, 0.0);
                                        var.set(GRB_DoubleAttr_UB, 0.0);
                                    }
                                    goto next_forward;
                                }
                            }
                        }
next_forward:;
                    }
                }
            }
        }

        // StaCoverage
        count = 0;
        for(const Station& s : sol_.getStations()) {
            for(std::size_t c = 0; c < coverage_->cellCount(); ++c) {
                ModelBase::addVar(Solution::Key::StaCoverage(&sol_, s, c), 0.0, 1.0, 0.0, GRB_BINARY);
                count++;
            } 
        }
        // populate StaConverage variables
        for(const Station& s : sol_.getStations()) {
            for(size_t c = 0; c < coverage_->cellCount(); ++c) {
                for(const auto q : sol_.getSources()) {
                    for(size_t t : sol_.getBlocks(t0, tf, q, s)) {
                        auto key = Solution::Key::StaActive(&sol_, q, s, t);
                        auto pv = sol_.getPointingVector(key);
                        if(coverage_->calculateCell(*pv) != c) continue;
                        if(*sol_.getSol(key)) {
                            auto var = *getVar(Solution::Key::StaCoverage(&sol_, s, c));
                            var.set(GRB_DoubleAttr_Start, 1.0);
                            goto next;
                        }
                    }
                }
next:;
            }
        }
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " sky coverage variables to model";
#else
        std::cout << "[info] Added " << count << " StaCoverage variables to model";
#endif

        // update the model to make sure variables are accessible
        model_->update();

        if(i == 0) {
            std::string dump = ModelBase::dump(0, sol_.blockCount_);
            std::cout << dump;
            dumps_.emplace_back(dump);
        }

        // freeze the previous optimization window...
        ModelBase::apply(tp, t0, [](GRBVar& var) {
            double val = var.get(GRB_DoubleAttr_Start);
            var.set(GRB_DoubleAttr_LB, val);
            var.set(GRB_DoubleAttr_UB, val);
        });

        // ... and the next optimization window
        ModelBase::apply(tf, tn, [](GRBVar& var) {
            double val = var.get(GRB_DoubleAttr_Start);
            var.set(GRB_DoubleAttr_LB, val);
            var.set(GRB_DoubleAttr_UB, val);
        });

        // optimize the window
        std::vector<std::tuple<std::string, GRBLinExpr, int, unsigned int>> obj = this->prepare({ .tp = tp, .tn = tn, .t0 = t0, .tf = tf });
        for(size_t i = 0; i < obj.size(); ++i) {
            if(i > 0) {
                if(std::get<2>(obj[i - 1]) == GRB_MAXIMIZE) {
                    model_->addConstr(std::get<1>(obj[i - 1]) + 1e-6 >= model_->get(GRB_DoubleAttr_ObjVal));
                } else {
                    model_->addConstr(std::get<1>(obj[i - 1]) <= model_->get(GRB_DoubleAttr_ObjVal) + 1e-6);
                }
            }
            model_->update();
            // model_->reset(0);
            model_->setObjective(std::get<1>(obj[i]), std::get<2>(obj[i]));
            model_->update();
            if(std::get<3>(obj[i]) > 0) {
                model_->getEnv().set(GRB_DoubleParam_TimeLimit, std::get<3>(obj[i]));
            }
#ifdef VIESCHEDPP_LOG
                BOOST_LOG_TRIVIAL( info ) << "Optimizing " << std::get<0>(obj[i]);
#else
                std::cout << "[info] Optimizing " << std::get<0>(obj[i]);
#endif
            model_->optimize();
            int status = model_->get(GRB_IntAttr_Status);
            if(status == GRB_INFEASIBLE) {
                model_->computeIIS();
                model_->write("/tmp/iis.ilp");
            }
            // error checking
            bool success = true;
            if(status == GRB_TIME_LIMIT) {
                if(model_->get(GRB_IntAttr_SolCount) == 0) {
                    success = false;
                }
            } else if(status != GRB_OPTIMAL && status != GRB_SUBOPTIMAL) {
                success = false;
            }

            if(!success) {
#ifdef VIESCHEDPP_LOG
                BOOST_LOG_TRIVIAL( info ) << "No optimal solution found between " << t0 * sol_.blockLength_ << " and " << tf * sol_.blockLength_;
#else
                std::cout << "[info] No optimal solution found between " << t0 * sol_.blockLength_ << " and " << tf * sol_.blockLength_;
#endif
            }
            this->post({ .tp = tp, .tn = tn, .t0 = t0, .tf = tf }, std::get<0>(obj[i]));
        }
        
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Completed optimization between " << t0 * sol_.blockLength_ << " and " << tf * sol_.blockLength_;
#else
        std::cout << "[info] Completed optimization between " << t0 * blockLength_ << " and " << tf * blockLength_;
#endif

        // copy results back into solution
        for(size_t t : sol_.getBlocks(t0, tf)) {
            for(const auto q : sol_.getSources()) {
                for(const Station& s : sol_.getStations(t, q)) {
                    Solution::Key key = Solution::Key::StaActive(&sol_, q, s, t);
                    auto sol = sol_.getSol(key);
                    auto var = ModelBase::getVar(key);
                    (*sol) = (var->get(GRB_DoubleAttr_X) > 0.5);
                }
                for(const Baseline& b : sol_.getBaselines(t, q)) {
                    Solution::Key key = Solution::Key::BlnActive(&sol_, q, b, t);
                    auto sol = sol_.getSol(key);
                    auto var = ModelBase::getVar(key);
                    (*sol) = (var->get(GRB_DoubleAttr_X) > 0.5);
                }
            }
        }

        std::string dump = ModelBase::dump(t0, tf);
        if(tf < tn) {
            const auto& prev = dumps_.back();
            std::cout << prev;
            std::cout << dump;
        }
        dumps_.emplace_back(dump);
    }

    for(const auto& dump : dumps_) {
        std::cout << dump;
    }
    
    return true;
#else // WITH_GUROBI
    return false;
#endif // WITH_GUROBI   
}

std::vector<Scan> ModelBase::optimize(std::vector<Scan>& scans) {
#ifdef WITH_GUROBI
    ModelBase::loadScans(scans);
    if(!ModelBase::optimize()) return {};
    std::vector<Scan> newScans = ModelBase::readScans();
    ModelBase::packScans(newScans);
    return newScans;
#else // WITH_GUROBI
    return {};
#endif // WITH_GUROBI
}
}

// ModelCoverage13 implementation
namespace VieVS {
std::size_t ModelCoverage13::cellCount(void) const noexcept { 
    return 13; 
}

std::size_t ModelCoverage13::calculateCell(const PointingVector& pv) const noexcept {
    constexpr double el_space = halfpi / 2.;

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
#ifdef WITH_GUROBI
boost::optional<GRBVar> ModelBase::getVar(const Solution::Key& key) const noexcept {
    try {
        return var_.at(key);
    } catch(...) {
        return boost::none;
    }
}

GRBVar& ModelBase::addVar(const Solution::Key& key, double lb, double ub, double obj, char vtype) {
    auto ret = var_.insert(std::make_pair(key, model_->addVar(lb, ub, obj, vtype, key.name)));
    if(!ret.second) throw UNREACHABLE;
    return ret.first->second;
}

std::set<std::tuple<const Observation*, size_t, size_t>> ModelBase::validateScan(const std::set<std::tuple<const Observation*, size_t, size_t>>& obs) {
    if(obs.empty()) return {};
    const auto q = sol_.sourceList_.getSource(std::get<0>(*obs.begin())->getSrcid());
    size_t minNumberOfSites = static_cast<size_t>(q->getPARA().minNumberOfSites);
    // find the combined largest span
    size_t t_start_comb = std::numeric_limits<size_t>::max();
    size_t t_end_comb = 0;
    for(const auto& data : obs) {
        const Observation* obs = std::get<0>(data);
        t_start_comb = std::min(t_start_comb, std::get<1>(data));
        t_end_comb = std::max(t_end_comb, std::get<2>(data));
    }   
    if(t_start_comb >= t_end_comb) return {};

    // construct a set of participating stations over each time segment
    std::vector<std::set<unsigned long>> sitesPerSegment{t_end_comb - t_start_comb};
    for(const auto& data : obs) {
        const Observation* obs = std::get<0>(data);
        size_t t_obs_start = std::get<1>(data);
        size_t t_obs_end = std::get<2>(data);
        for(size_t t_obs = t_obs_start - t_start_comb; t_obs < t_obs_end - t_start_comb; ++t_obs) {
            sitesPerSegment[t_obs].emplace(obs->getStaid1());
            sitesPerSegment[t_obs].emplace(obs->getStaid2());
        }
    }   

#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( warning ) << "Finished building sitesPerSegment!";
#else
    std::cout << "[warning] Finished building sitesPerSegment!";
#endif

    // find the span where the minNumberOfSites is respected
    size_t first = 0;
    while(first < sitesPerSegment.size() && sitesPerSegment[first].size() < minNumberOfSites) ++first;
    if(first == sitesPerSegment.size()) return {}; // skip to next scan, nothing here can be added to the warm start
    size_t last = first;
    while(last + 1 < sitesPerSegment.size() && sitesPerSegment[last + 1].size() >= minNumberOfSites) ++last;
    size_t new_start = t_start_comb + first;
    size_t new_end   = t_start_comb + last + 1;
    t_start_comb = new_start;
    t_end_comb   = new_end;
    
#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( warning ) << "Finished finding final span!";
#else
    std::cout << "[warning] Finished finding final span!";
#endif  

    // one more check for SNR
    std::set<std::tuple<const Observation*, size_t, size_t>> obsValidInner;
    for(const auto& data : obs) {
        const Observation* obs = std::get<0>(data);
        const Baseline& b = sol_.network_.getBaseline(obs->getBlid());
        const Station& s1 = sol_.network_.getStation(b.getStaid1());
        const Station& s2 = sol_.network_.getStation(b.getStaid2());
        size_t t_start_curr = std::max(t_start_comb, std::get<1>(data));
        size_t t_end_curr = std::min(t_end_comb, std::get<2>(data));
        bool viable = true;
        for(size_t t = t_start_curr; t < t_end_curr; ++t) {
            size_t dur = sol_.getMinObs(Solution::Key::BlnActive(&sol_, q, b, t_start_curr));
            if(t_end_curr - t_start_curr < dur) {
                viable = false;
                break;
            }
        }
        if(viable) {
            obsValidInner.emplace(obs, std::max(t_start_comb, std::get<1>(data)), std::min(t_end_comb, std::get<2>(data)));
        }
    }        
    return obsValidInner;
}

void ModelBase::loadScans(std::vector<Scan>& scans) {
    if(scans.empty()) {
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Skipped loading MIP solution. No scans provided";
#else
        std::cout << "[info] Skipped loading MIP solution. No scans provided";
#endif 
        return;
    }

    // populate starting values from given scans
    for(const Scan& scan : scans) {
        const ScanTimes& scanTimes = scan.getTimes();

        // populate BlnActive variables
        std::map<unsigned long, std::set<std::tuple<const Observation*, size_t, size_t>>> obsValid;
        for(const Observation& obs : scan.getObservations()) {
            const auto q = sol_.sourceList_.getSource(obs.getSrcid());
            if(sol_.sourceMask_.count(q->getId()) == 0) continue;
            const Baseline& b = sol_.network_.getBaseline(obs.getBlid());
            Station& s1 = sol_.network_.refStation(b.getStaid1());
            Station& s2 = sol_.network_.refStation(b.getStaid2());

            // observation start blocks
            size_t t10 = scanTimes.getObservingTime(Timestamp::start) / sol_.blockLength_; // (scanTimes.getObservingTime(s1.getId())) / sol_.blockLength_;
            size_t t20 = scanTimes.getObservingTime(Timestamp::start) / sol_.blockLength_; // (scanTimes.getObservingTime(s2.getId())) / sol_.blockLength_;

            if(t10 > sol_.blockCount_ || t20 > sol_.blockCount_) continue;

            // the number of blocks each station is observing
            size_t t1f = scanTimes.getObservingTime(s1.getId(), Timestamp::end) / sol_.blockLength_;
            size_t t2f = scanTimes.getObservingTime(s2.getId(), Timestamp::end) / sol_.blockLength_;
            t1f = std::min(t1f, sol_.blockCount_);
            t2f = std::min(t2f, sol_.blockCount_);
            while(t1f > 0 && !sol_.getSol(Solution::Key::StaActive(&sol_, q, s1, t1f))) --t1f;
            while(t2f > 0 && !sol_.getSol(Solution::Key::StaActive(&sol_, q, s2, t2f))) --t2f;
            // check that these are possible, if they aren't then advance by one segment
            if(sol_.getSol(Solution::Key::StaActive(&sol_, q, s1, t10))) {
                for(const auto q2 : sol_.getSources()) {
                    if(q->getId() == q2->getId()) continue;
                    if(auto sol = sol_.getSol(Solution::Key::StaActive(&sol_, q2, s1, t10))) {
                        if(*sol) {
                            ++t10;
                            break;
                        }
                    }
                }
            }
            for(size_t t = t10; t-- > 0;) {
                for(const auto q2 : sol_.getSources()) {
                    if(q->getId() == q2->getId()) continue;
                    if(auto sol = sol_.getSol(Solution::Key::StaActive(&sol_, q2, s1, t))) {
                        if(*sol) {
                            // then we need to check slew time
                            size_t t_slew = sol_.getSlew(Solution::Key::StaActive(&sol_, q2, s1, t + 1), 
                                Solution::Key::StaActive(&sol_, q, s1, t10));
                            if(t10 == 30) {
                                std::cout << s1.getName() << " " << s2.getName() << " t: " << t << " t10: " << t10 << " t_slew: " << t_slew << std::endl;
                            }
                            if(t10 - t <= t_slew) t10 += t_slew + 1 - (t10 - t);
                            goto terminate_s1;
                        }
                    }
                }
            }
terminate_s1:;
            while(t10 < sol_.blockCount_ && !sol_.getSol(Solution::Key::StaActive(&sol_, q, s1, t10))) ++t10;

            // start by flooring the scan starts, then check if there is enough slew time
            if(sol_.getSol(Solution::Key::StaActive(&sol_, q, s2, t20))) {
                for(const auto q2 : sol_.getSources()) {
                    if(q->getId() == q2->getId()) continue;
                    if(auto sol = sol_.getSol(Solution::Key::StaActive(&sol_, q2, s2, t20))) {
                        if(*sol) {
                            ++t20;
                            break;
                        }
                    }
                }
            }
            for(size_t t = t20; t-- > 0;) {
                for(const auto q2 : sol_.getSources()) {
                    if(q->getId() == q2->getId()) continue;
                    if(auto sol = sol_.getSol(Solution::Key::StaActive(&sol_, q2, s2, t))) {
                        if(*sol) {
                            // then we need to check slew time
                            size_t t_slew = sol_.getSlew(Solution::Key::StaActive(&sol_, q2, s2, t + 1), 
                                Solution::Key::StaActive(&sol_, q, s2, t20));
                            if(t20 == 30) {
                                std::cout << s1.getName() << " " << s2.getName() << " t: " << t << " t20: " << t20 << " t_slew: " << t_slew << std::endl;
                            }
                            if(t20 - t <= t_slew) t20 += t_slew + 1 - (t20 - t);
                            goto terminate_s2;
                        }
                    }
                }
            }
terminate_s2:;
            while(t20 < sol_.blockCount_ && !sol_.getSol(Solution::Key::StaActive(&sol_, q, s2, t20))) ++t20;
            
            // check if any starts exceed the ends
            if(t10 >= t1f) continue;
            if(t20 >= t2f) continue;
            
            size_t t_start = std::max(t10, t20);
            size_t t_end = std::min(t1f, t2f);
            if(t_start >= t_end) continue;

            size_t maxScan = (std::min(std::min(s1.getPARA().maxScan, s2.getPARA().maxScan), q->getPARA().maxScan) + sol_.blockLength_ - 1) / sol_.blockLength_;
            if(t_end - t_start > maxScan) continue;

            bool viable = true;
            for(size_t t = t_start; t < t_end; ++t) {
                size_t dur = sol_.getMinObs(Solution::Key::BlnActive(&sol_, q, b, t));
                if(t_end - t_start < dur) {
                    viable = false;
                    break;
                }
            }
            
            if(viable) {
                if(obsValid.count(q->getId()) == 0) {
                    obsValid.emplace(q->getId(), std::set<std::tuple<const Observation*, size_t, size_t>>{});
                }
                obsValid[q->getId()].emplace(&obs, t_start, t_end);
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( warning ) << "Finished finding valid observations!";
#else
        std::cout << "[warning] Finished finding valid observations!";
#endif

        for(auto obsValidQ : obsValid) {
            auto q = sol_.sourceList_.getSource(std::get<0>(obsValidQ));
            std::set<std::tuple<const Observation*, size_t, size_t>> obsValidCurr;
            std::set<std::tuple<const Observation*, size_t, size_t>> obsValidNext = std::get<1>(obsValidQ);
            do {
                obsValidCurr = obsValidNext;
                obsValidNext = ModelBase::validateScan(obsValidCurr);
            } while(obsValidNext.size() < obsValidCurr.size() && !obsValidNext.empty());

            if(obsValidNext.empty()) continue;

            // finally, we have the combined start and end
            // now we have to add each observation to the warm start
            for(const auto& data : obsValidNext) {
                const Observation* obs = std::get<0>(data);
                const Baseline& b = sol_.network_.getBaseline(obs->getBlid());
                const Station& s1 = sol_.network_.getStation(b.getStaid1());
                const Station& s2 = sol_.network_.getStation(b.getStaid2());
                for(size_t t = std::get<1>(data); t < std::get<2>(data); ++t) {
                    *sol_.getSol(Solution::Key::StaActive(&sol_, q, s1, t)) = true;
                    *sol_.getSol(Solution::Key::StaActive(&sol_, q, s2, t)) = true;
                    *sol_.getSol(Solution::Key::BlnActive(&sol_, q, b, t)) = true;
                }
            }
        }
    }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( warning ) << "Finished adding observations!";
#else
        std::cout << "[warning] Finished adding observations!";
#endif

#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "Loaded preliminary result into ILP model";
#else
    std::cout << "[info] Loaded preliminary result into ILP model";
#endif
}

std::vector<Scan> ModelBase::readScans(void) const noexcept {
    std::map<unsigned long, char> qId;
    for(const unsigned long q : sol_.sourceMask_) {
        qId.insert(std::make_pair(q, static_cast<char>(qId.size() + '!')));
    }

    std::map<unsigned int, std::vector<Interval>> intervals;
    for(Station& s : sol_.network_.refStations()) intervals[s.getId()] = {};

    for(size_t t = 0; t <= sol_.blockCount_; ++t) {
        for(Station& s : sol_.network_.refStations()) {
            std::vector<Interval>& s_intervals = intervals[s.getId()];
            for(const auto q : sol_.getSources()) {
                if(auto sol = sol_.getSol(Solution::Key::StaActive(&sol_, q, s, t))) {
                    if(*sol) {
                        if(s_intervals.empty()) {
                            s_intervals.emplace_back(t, t, s.getId(), q->getId());
                        } else {
                            Interval& last = s_intervals.back();
                            if(last.q == q->getId() && last.end + 1 == t) {
                                last.end++;
                            } else {
                                s_intervals.emplace_back(t, t, s.getId(), q->getId());
                            }
                        }
                    }
                }
            }
        }
    }

    std::map<unsigned int, std::vector<Interval>> bySource;
    for(const auto& entry : intervals) {
        for(const Interval& interval : entry.second) {
            bySource[interval.q].push_back(interval);
        }
    }

    std::vector<std::pair<unsigned int, ActiveGroup>> groups;
    for(auto& sourceEntry : bySource) {
        const unsigned int q = sourceEntry.first;
        auto& sourceIntervals = sourceEntry.second;
        std::sort(sourceIntervals.begin(), sourceIntervals.end(), [](const Interval& a, const Interval& b) {
            if(a.start != b.start) return a.start < b.start;
            if(a.end != b.end) return a.end < b.end;
            return a.s < b.s;
        });

        ActiveGroup current;
        bool haveGroup = false;

        size_t i = 0;
        while (i < sourceIntervals.size()) {
            const size_t epoch = sourceIntervals[i].start;
            size_t j = i;
            while(j < sourceIntervals.size() &&
                sourceIntervals[j].start == epoch) {
                ++j;
            }

            // does a station join a pre-existing scan
            bool lateJoin = false;
            if(haveGroup) {
                bool observationActive = false;
                for(const auto& entry : current.stations) {
                    if(entry.second.end >= epoch) {
                        observationActive = true;
                        break;
                    }
                }
                if(observationActive) {
                    for(size_t k = i; k < j; ++k) {
                        if(current.stations.find(sourceIntervals[k].s) ==
                            current.stations.end()) {
                            lateJoin = true;
                            break;
                        }
                    }
                }
            }

            if(haveGroup && !lateJoin) {
                bool observationActive = false;
                for(const auto& entry : current.stations) {
                    if(entry.second.end >= epoch) {
                        observationActive = true;
                        break;
                    }
                }
                if(!observationActive) {
                    groups.emplace_back(q, current);
                    current = ActiveGroup{};
                    haveGroup = false;
                }
            }

            if(lateJoin) {
                ActiveGroup continuation;
                // close ongoing scan at this point, then start a new one with this new participant
                for(auto it = current.stations.begin();
                    it != current.stations.end();) {
                    Interval& oldInterval = it->second;
                    if(oldInterval.end >= epoch) {
                        Interval continued = oldInterval;
                        continued.start = epoch;
                        oldInterval.end = epoch - 1;
                        continuation.stations.emplace(it->first, continued);
                        continuation.continuations.insert(it->first);
                    }
                    if (oldInterval.start > oldInterval.end) {
                        it = current.stations.erase(it);
                    } else {
                        ++it;
                    }
                }

                groups.emplace_back(q, current);
                current = std::move(continuation);
            }

            // add (and readd) stations
            for(size_t k = i; k < j; ++k) {
                const Interval& interval = sourceIntervals[k];
                auto it = current.stations.find(interval.s);
                if(it == current.stations.end()) {
                    current.stations.emplace(interval.s, interval);
                } else {
                    // Merge adjacent or overlapping intervals for
                    // the same station, if encountered.
                    it->second.start = std::min(it->second.start, interval.start);
                    it->second.end = std::max(it->second.end, interval.end);
                }
            }

            haveGroup = true;
            i = j;
        }

        // handle final group
        if (haveGroup) {
            groups.emplace_back(q, current);
        }
    }

    std::sort(groups.begin(), groups.end(), [](const auto& a, const auto& b) {
        size_t a_first = std::numeric_limits<size_t>::max();
        for (const auto& entry : std::get<1>(a).stations) a_first = std::min(a_first, entry.second.start);
        size_t b_first = std::numeric_limits<size_t>::max();
        for (const auto& entry : std::get<1>(b).stations) b_first = std::min(b_first, entry.second.start);
        return a_first < b_first;
    });

    std::vector<Scan> scans;

    std::map<unsigned int, std::pair<unsigned int, size_t>> eols;
    for(const auto& entry : groups) {
        // printGroup(entry.first, entry.second);
        const auto& q_id = std::get<0>(entry);
        const auto q = sol_.sourceList_.getSource(q_id);
        const auto& group = std::get<1>(entry);

        std::vector<PointingVector> pointingVectors;
        std::vector<PointingVector> pointingVectorsEnd;
        for(const auto& group_entry : group.stations) {
            const auto& s_id = std::get<0>(group_entry);
            const auto& s_interval = std::get<1>(group_entry);
            PointingVector s_pv_0{s_id, q_id};
            s_pv_0.setTime(s_interval.start * sol_.blockLength_);
            sol_.network_.refStation(s_id).calcAzEl_rigorous(q, s_pv_0);
            pointingVectors.emplace_back(s_pv_0);
            PointingVector s_pv_f{s_id, q_id};
            s_pv_f.setTime((s_interval.end + 1) * sol_.blockLength_);
            sol_.network_.refStation(s_id).calcAzEl_rigorous(q, s_pv_f);
            pointingVectorsEnd.emplace_back(s_pv_f);
        }
        std::vector<unsigned int> endOfLastScan;
        for(const auto& group_entry : group.stations) {
            const auto& s_id = std::get<0>(group_entry);
            if(eols.count(s_id) > 0 && group.continuations.count(s_id) == 0) {
                endOfLastScan.emplace_back((eols.at(s_id).second + 1) * sol_.blockLength_);
            } else {
                const auto& s_interval = std::get<1>(group_entry);
                endOfLastScan.emplace_back(s_interval.start * sol_.blockLength_);
            }
        }
        Scan scan{pointingVectors, endOfLastScan, Scan::ScanType::standard};
        scan.setPointingVectorsEndtime(pointingVectorsEnd);
        size_t idx = 0;
        for(const auto& group_entry : group.stations) {
            const auto& s_id = std::get<0>(group_entry);
            auto& s = sol_.network_.refStation(s_id);
            const auto& s_interval = std::get<1>(group_entry);
            if(eols.count(s_id) > 0 && group.continuations.count(s_id) == 0) {
                const auto& q_id_0 = eols.at(s_id).first;
                const auto q_0 = sol_.sourceList_.getSource(q_id_0);
                PointingVector s_pv_0{s_id, q_id_0};
                s_pv_0.setTime((eols.at(s_id).second + 1) * sol_.blockLength_);
                s.calcAzEl_rigorous(q_0, s_pv_0);
                PointingVector s_pv_f{s_id, q_id};
                s_pv_f.setTime(s_interval.start * sol_.blockLength_);
                s.calcAzEl_rigorous(q, s_pv_f);
                unsigned int slew = s.getAntenna().slewTime(s_pv_0, s_pv_f);
                unsigned int idle = (s_pv_f.getTime() - s_pv_0.getTime()) - slew - s.getPARA().preob - s.getPARA().systemDelay;
                // slew time is from (eols + systemDelay) to (s_interval.start * blockLength - preob)
                scan.addTimes(idx, s.getPARA().systemDelay, slew, idle, s.getPARA().preob);
                scan.referenceTime().setPreobTime(idx, s.getPARA().preob);
            } else {
                scan.addTimes(idx, 0, 0, 0);
            }
            idx++;
        }

        std::vector<unsigned int> scanTimes;
        for(const auto& group_entry : group.stations) {
            const auto& s_id = std::get<0>(group_entry);
            auto& s = sol_.network_.refStation(s_id);
            const auto& s_interval = std::get<1>(group_entry);
            scanTimes.emplace_back((s_interval.end - s_interval.start + 1) * sol_.blockLength_);
        }
        scan.referenceTime().setObservingTimes(scanTimes, false);
        scan.constructObservations(sol_.network_, q);

        scans.emplace_back(scan);

        for(const auto& group_entry : group.stations) {
            const auto& s_id = group_entry.first;
            const auto& s_interval = group_entry.second;
            eols[s_id] = {q_id, s_interval.end};
        }
    }

    return scans;
}

void dumpTimes(const ScanTimes& times, unsigned int idx) {
    unsigned int start = times.getFieldSystemTime(idx);
    unsigned int slew = times.getSlewTime(idx);
    unsigned int idle = times.getIdleTime(idx);
    unsigned int preob = times.getPreobTime(idx);
    unsigned int obs = times.getObservingTime(idx, Timestamp::start);
    unsigned int end = times.getObservingTime(idx, Timestamp::end);
    std::cout << start << ", " << slew << ", " << idle << ", " << preob << ", " << obs << ", " << end << std::endl;
}

void ModelBase::packScans(std::vector<Scan>& scans) const noexcept {
    std::map<unsigned int, std::vector<size_t>> scansByStation;
    for(const Station& station : sol_.network_.getStations()) {
        scansByStation.emplace(station.getId(), std::vector<size_t>{});
    }

    for(size_t i = 0; i < scans.size(); ++i) {
        for(unsigned int s_id : scans[i].getStationIds()) {
            auto it = scansByStation.find(s_id);
            if(it != scansByStation.end()) it->second.push_back(i);
        }
    }

    // sort each station's scan by time
    for(auto& entry : scansByStation) {
        const unsigned int s_id = entry.first;
        auto& stationScans = entry.second;
        std::sort(stationScans.begin(), stationScans.end(), [&scans, s_id](size_t a, size_t b) {
            const auto s_idx_a = *scans[a].findIdxOfStationId(s_id);
            const auto s_idx_b = *scans[b].findIdxOfStationId(s_id);
            const unsigned int start_a = scans[a].getTimes().getObservingTime(s_idx_a);
            const unsigned int start_b = scans[b].getTimes().getObservingTime(s_idx_b);
            return start_a != start_b ? start_a < start_b : a < b;
        });
    }

    // build downtime map, same as in Solution.cpp
    // there's a better way to handle downtime, but at least it works
    std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int>>> downtime;
    for(Station& station : sol_.network_.refStations()) {
        const unsigned int s_id = station.getId();

        std::map<size_t, bool> events;
        for(const auto& event : station.refParaForMultiScheduling()) {
            if(event.PARA.available) {
                events.emplace((event.time + sol_.blockLength_ - 1) / sol_.blockLength_, true);
            } else {
                events.emplace(event.time / sol_.blockLength_, false);
            }
        }

        bool available = true;
        size_t prev = 0;
        bool has_prev = false;
        for(const auto& ev : events) {
            const auto& block = std::get<0>(ev);
            const auto& available_new = std::get<1>(ev);
            if(has_prev && !available) {
                const auto block_start = static_cast<unsigned long long>(prev) * sol_.blockLength_;
                const auto block_end = static_cast<unsigned long long>(block) * sol_.blockLength_;
                downtime[s_id].emplace_back(static_cast<unsigned int>(block_start), static_cast<unsigned int>(block_end));
            }
            available = available_new;
            prev = block;
            has_prev = true;
        }

        if(has_prev && !available) {
            const auto block_start = static_cast<unsigned long long>(prev) * sol_.blockLength_;
            downtime[s_id].emplace_back(static_cast<unsigned int>(block_start), std::numeric_limits<unsigned int>::max());
        }
    }

    // helper to check if a scan intersects a down period
    auto crossesDowntime = [&downtime](unsigned int s_id, unsigned int start, unsigned int end) -> bool {
        const auto it = downtime.find(s_id);
        if(it == downtime.end() || start >= end) return false;
        for(const auto& entry : it->second) {
            if(start < std::get<1>(entry) && end > std::get<0>(entry)) return true;
        }
        return false;
    };

    // update the start or end pointing vector of a station in a scan
    auto refreshPointingVector = [this](Scan& scan, unsigned int s_id, Timestamp ts) {
        const unsigned long s_idx = *scan.findIdxOfStationId(s_id);
        PointingVector& pv = scan.referencePointingVector(s_idx, ts);
        pv.setTime(scan.getTimes().getObservingTime(s_idx, ts));
        sol_.network_.refStation(s_id).calcAzEl_rigorous(sol_.sourceList_.getSource(scan.getSourceId()), pv);
    };

    // Record same-source extensions for collective validation in the
    // clipping phase. No subsequent-scan overhead is removed during growth.
    // because VieVS::Scan objects assume a unified start time,
    // its not possible to represent a single station exiting and entering a scan at multiple points
    // instead, we eliminate all the overhead time from a 'continuation' (i.e. another scan of the same source directly after the first)
    // and treat it as an extension
    struct ScanExtension {
        size_t currentScan;
        size_t nextScan;
        unsigned int stationId;
        unsigned int originalEnd;
        unsigned int nextStart;
        PointingVector originalEndPV;
    };
    std::vector<ScanExtension> extensions;
    using ExtensionKey = std::pair<size_t, unsigned int>;
    std::map<ExtensionKey, size_t> extensionByScanStation;
    
    // grow observations
    for(const auto& entry : scansByStation) {
        const unsigned int s_id = entry.first;
        const auto& stationScans = entry.second;

        Station& station = sol_.network_.refStation(s_id);

        bool isFirstScan = true;

        const unsigned int constTimes = station.getPARA().systemDelay + station.getPARA().preob;
        for(size_t i = 0; i + 1 < stationScans.size(); ++i) {
            const size_t currentScanIdx = stationScans[i];
            const size_t nextScanIdx = stationScans[i + 1];
            Scan& current = scans[currentScanIdx];
            Scan& next = scans[nextScanIdx];
            ScanTimes& currentTimes = current.referenceTime();
            ScanTimes& nextTimes = next.referenceTime();

            const unsigned long currentIdx = *current.findIdxOfStationId(s_id);
            const unsigned long nextIdx = *next.findIdxOfStationId(s_id);

            const unsigned int originalEnd = currentTimes.getObservingTime(currentIdx, Timestamp::end);
            const unsigned int nextStart = nextTimes.getObservingTime(nextIdx, Timestamp::start);
            if(nextStart <= originalEnd) continue;

            const unsigned int nextIdle = nextTimes.getIdleDuration(nextIdx);

            const bool sameSource = current.getSourceId() == next.getSourceId();
            if(sameSource) {
                if(crossesDowntime(s_id, originalEnd, nextStart)) continue;
                extensions.push_back({
                    currentScanIdx,
                    nextScanIdx,
                    s_id,
                    originalEnd,
                    nextStart,
                    current.getPointingVector(static_cast<int>(currentIdx), Timestamp::end)
                });

                extensionByScanStation.emplace(ExtensionKey{currentScanIdx, s_id}, extensions.size() - 1);
                currentTimes.setObservingTime(currentIdx, nextStart, Timestamp::end);
                refreshPointingVector(current, s_id, Timestamp::end);

                continue;
            }

            // if the sources are different, do not consume actually-needed overhead like slew or preob
            if(nextIdle == 0) continue;

            const PointingVector originalEndPV = current.getPointingVector(static_cast<int>(currentIdx), Timestamp::end);

            const bool skipSlew = isFirstScan;
            auto slewIsFeasible = [&](unsigned int candidateEnd) -> bool {
                if(candidateEnd >= nextStart) return false;
                currentTimes.setObservingTime(currentIdx, candidateEnd, Timestamp::end);

                refreshPointingVector(current, s_id, Timestamp::end);
                refreshPointingVector(next, s_id, Timestamp::start);

                const PointingVector& thisEnd = current.getPointingVector(static_cast<int>(currentIdx), Timestamp::end);
                const PointingVector& nextStartPV = next.getPointingVector(static_cast<int>(nextIdx), Timestamp::start);
                const auto oslewtime = skipSlew ? boost::optional<unsigned int>(0) : station.slewTime(
                    thisEnd, nextStartPV, currentTimes.getObservingDuration(currentIdx), true);
                if(!oslewtime.is_initialized()) return false;
                
                // this mirrors the validation logic from checkAndStatistics
                const unsigned long long required = static_cast<unsigned long long>(*oslewtime) + constTimes;
                const unsigned long long available = static_cast<unsigned long long>(nextStart) - candidateEnd;
                return available >= required;
            };

            unsigned int maxAllowedEnd = std::numeric_limits<unsigned int>::max();

            const auto downtimeIt = downtime.find(s_id);
            if(downtimeIt != downtime.end()) {
                for(const auto& entry : downtimeIt->second) {
                    if(std::get<0>(entry) >= originalEnd) {
                        maxAllowedEnd = std::get<0>(entry);
                        break;
                    }
                }
            }

            unsigned int shift = nextIdle;

            // check to make sure we aren't extending into a downtime period
            if(maxAllowedEnd <= originalEnd) continue;

            shift = std::min(shift, maxAllowedEnd - originalEnd);

            // iteratively try to improve observation time while keeping slew valid
            bool converged = false;
            constexpr unsigned int maxIterations = 32;
            for(unsigned int iteration = 0; iteration < maxIterations; ++iteration) {
                const unsigned long long candidateEndTemp = static_cast<unsigned long long>(originalEnd) + shift;
                if(candidateEndTemp >= nextStart || candidateEndTemp > std::numeric_limits<unsigned int>::max()) shift = 0;

                const unsigned int candidateEnd = originalEnd + shift;
                if(!slewIsFeasible(candidateEnd)) {
                    currentTimes.setObservingTime(currentIdx, originalEnd, Timestamp::end);
                    refreshPointingVector(current, s_id, Timestamp::end);
                    if(shift == 0) break;
                    --shift;
                    continue;
                }

                const PointingVector& thisEnd = current.getPointingVector(static_cast<int>(currentIdx), Timestamp::end);
                const PointingVector& nextStartPV = next.getPointingVector(static_cast<int>(nextIdx), Timestamp::start);

                const auto oslewtime = station.slewTime(thisEnd, nextStartPV,
                    currentTimes.getObservingDuration(currentIdx));
                if(!oslewtime.is_initialized()) {
                    currentTimes.setObservingTime(currentIdx, originalEnd, Timestamp::end);
                    refreshPointingVector(current, s_id, Timestamp::end);
                    shift = 0;
                    break;
                }

                const unsigned long long required = static_cast<unsigned long long>(*oslewtime) + constTimes;
                const unsigned long long available = static_cast<unsigned long long>(nextStart) - candidateEnd;

                const unsigned int safeShift = available >= required ? static_cast<unsigned int>(std::min<unsigned long long>(nextIdle, static_cast<unsigned long long>(shift) + available - required)) : 0;
                if(safeShift == shift) {
                    converged = true;
                    break;
                }
                shift = std::min(shift, safeShift);

                currentTimes.setObservingTime(currentIdx, originalEnd, Timestamp::end);
                refreshPointingVector(current, s_id, Timestamp::end);
            }

            if(!converged) {
                currentTimes.setObservingTime(currentIdx, originalEnd, Timestamp::end);
                PointingVector& pv = current.referencePointingVector(static_cast<int>(currentIdx), Timestamp::end);
                pv = originalEndPV;
                continue;
            }

            const unsigned int finalEnd = originalEnd + shift;
            if(!slewIsFeasible(finalEnd)) {
                currentTimes.setObservingTime(currentIdx, originalEnd, Timestamp::end);
                PointingVector& pv = current.referencePointingVector(static_cast<int>(currentIdx), Timestamp::end);
                pv = originalEndPV;
                continue;
            }

            nextTimes.setIdleTime(nextIdx, nextIdle - shift);
            isFirstScan = false;
        }
    }

    // given potential ScanExtensions, 
    // we have to accept or reject them depending on whether or not minNumberOfSites is respected
    std::vector<bool> accepted(extensions.size(), true);

    bool changed = true;
    while(changed) {
        changed = false;
        for(size_t e = 0; e < extensions.size(); ++e) {
            if(!accepted[e]) continue;

            const ScanExtension& extension = extensions[e];
            const Scan& current = scans[extension.currentScan];

            const unsigned int minSites = sol_.sourceList_.getSource(current.getSourceId())->getPARA().minNumberOfSites;

            unsigned int participatingStations = 0;
            for(unsigned int otherId : current.getStationIds()) {
                const unsigned long otherIdx = *current.findIdxOfStationId(otherId);
                const unsigned int otherStart = current.getTimes().getObservingTime(otherIdx, Timestamp::start);
                unsigned int otherEnd = current.getTimes().getObservingTime(otherIdx, Timestamp::end);

                const auto proposalIt = extensionByScanStation.find(ExtensionKey{extension.currentScan, otherId});
                if(proposalIt != extensionByScanStation.end()) {
                    const size_t otherProposal = proposalIt->second;
                    if(accepted[otherProposal]) {
                        otherEnd = extensions[otherProposal].nextStart;
                    } else {
                        otherEnd = extensions[otherProposal].originalEnd;
                    }
                }

                if(otherStart < extension.nextStart && otherEnd >= extension.nextStart) ++participatingStations;
            }

            if(participatingStations < minSites) {
                accepted[e] = false;
                changed = true;
            }
        }
    }

    // incorporate accepted extensions into the schedule
    for(size_t e = 0; e < extensions.size(); ++e) {
        const ScanExtension& extension = extensions[e];

        Scan& current = scans[extension.currentScan];
        Scan& next = scans[extension.nextScan];
        const unsigned long currentIdx = *current.findIdxOfStationId(extension.stationId);
        const unsigned long nextIdx = *next.findIdxOfStationId(extension.stationId);

        if(!accepted[e]) {
            current.referenceTime().setObservingTime(currentIdx, extension.originalEnd, Timestamp::end);
            PointingVector& pv = current.referencePointingVector(currentIdx, Timestamp::end);
            pv = extension.originalEndPV;
            continue;
        }

        // the extension is valid, so we can remove overhead from the subsequent observation
        next.referenceTime().removeOverheadBeforeObservation(static_cast<int>(nextIdx), extension.nextStart);
    }

    // clip lengths such that no station observes alone
    for(Scan& scan : scans) {
        ScanTimes& times = scan.referenceTime();
        if(scan.getStationIds().size() < 2) continue;

        std::vector<std::pair<unsigned int, unsigned long>> observingEnds;
        observingEnds.reserve(scan.getStationIds().size());
        for(unsigned int s_id : scan.getStationIds()) {
            const unsigned long s_idx = *scan.findIdxOfStationId(s_id);
            observingEnds.emplace_back(times.getObservingTime(s_idx, Timestamp::end), s_idx);
        }

        std::sort(observingEnds.begin(), observingEnds.end(), [](const auto& a, const auto& b) {
            return a.first > b.first;
        });

        const unsigned int clipEnd = observingEnds[1].first;
        for(const auto& item : observingEnds) {
            const unsigned int end = item.first;
            const unsigned long s_idx = item.second;
            if(end <= clipEnd) continue;

            const unsigned int start = times.getObservingTime(s_idx, Timestamp::start);
            if(clipEnd > start) times.setObservingTime(s_idx, clipEnd, Timestamp::end);
        }
    }

    for(Scan& scan : scans) {
        for(unsigned int s_id : scan.getStationIds()) {
            refreshPointingVector(scan, s_id, Timestamp::start);
            refreshPointingVector(scan, s_id, Timestamp::end);
        }
        scan.updateObservingTime();
    }

    // final validation check for overlaps
    std::map<size_t, std::set<unsigned int>> stationsToRemove;
    for(const auto& entry : scansByStation) {
        const unsigned int s_id = entry.first;
        std::vector<size_t> stationScans;
        for(size_t i = 0; i < scans.size(); ++i) {
            if(scans[i].findIdxOfStationId(s_id)) stationScans.push_back(i);
        }
        // sort by scan start
        std::sort(stationScans.begin(), stationScans.end(), [&scans, s_id](size_t a, size_t b) {
            const unsigned int a_start = scans[a].getTimes().getFieldSystemTime(
                *scans[a].findIdxOfStationId(s_id), Timestamp::start);
            const unsigned int b_start = scans[b].getTimes().getFieldSystemTime(
                *scans[b].findIdxOfStationId(s_id), Timestamp::start);
            return a_start != b_start ? a_start < b_start : a < b;
        });

        bool havePreviousObservation = false;
        unsigned int lastRetainedObservingEnd = 0;
        for(size_t scanIdx : stationScans) {
            const Scan& scan = scans[scanIdx];
            const unsigned long s_idx = *scan.findIdxOfStationId(s_id);
            const unsigned int fieldSystemStart = scan.getTimes().getFieldSystemTime(s_idx, Timestamp::start);
            const unsigned int observingEnd = scan.getTimes().getObservingTime(s_idx, Timestamp::end);
            // reject participation if the scan start is before previous scan (for this station) finished
            if(havePreviousObservation && fieldSystemStart < lastRetainedObservingEnd) {
                stationsToRemove[scanIdx].insert(s_id);
                continue;
            }
            // store as this station's last scan
            lastRetainedObservingEnd = observingEnd;
            havePreviousObservation = true;
        }
    }

    // remove stations from observations in reverse order if some overlap has been produced (which shouldn't happen),
    // but I haven't had time to validate every edge case, so I'm keeping it here
    for(const auto& entry : stationsToRemove) {
        Scan& scan = scans[entry.first];
        const auto source = sol_.sourceList_.getSource(scan.getSourceId());
        std::vector<int> indicesToRemove;
        for(unsigned int s_id : entry.second) {
            const auto s_idx = scan.findIdxOfStationId(s_id);
            if(s_idx) indicesToRemove.push_back(static_cast<int>(*s_idx));
        }
        std::sort(indicesToRemove.begin(), indicesToRemove.end(), std::greater<int>());
        for(int s_idx : indicesToRemove) scan.removeStation(s_idx, source);
    }

    // now, remove any scans that violate minNumberOfSites
    std::vector<Scan> validScans;
    validScans.reserve(scans.size());
    for(Scan& scan : scans) {
        const unsigned int minSites = sol_.sourceList_.getSource(scan.getSourceId())->getPARA().minNumberOfSites;
        if(scan.getStationIds().size() < minSites) continue;
        validScans.push_back(std::move(scan));
    }
    scans = std::move(validScans);

    // once again, clip the longest observing station to the second longest, 
    // in case a removed station was the co-observer of the previous longest observing station
    for(Scan& scan : scans) {
        const unsigned int minSites = sol_.sourceList_.getSource(scan.getSourceId())->getPARA().minNumberOfSites;
        const size_t stationCount = scan.getStationIds().size();
        if(stationCount < 2) continue;

        std::vector<std::pair<unsigned int, unsigned long>> observingEnds;
        observingEnds.reserve(stationCount);
        for(unsigned int s_id : scan.getStationIds()) {
            const unsigned long s_idx = *scan.findIdxOfStationId(s_id);
            observingEnds.emplace_back(scan.getTimes().getObservingTime(s_idx, Timestamp::end), s_idx);
        }
        std::sort(observingEnds.begin(), observingEnds.end(), [](const auto& a, const auto& b) {
            return a.first > b.first;
        });

        const unsigned int clipEnd = minSites <= stationCount ? observingEnds[minSites - 1].first : 0;

        ScanTimes& times = scan.referenceTime();
        for(const auto& item : observingEnds) {
            const unsigned int end = item.first;
            const unsigned long s_idx = item.second;
            if(end <= clipEnd) continue;

            const unsigned int start = times.getObservingTime(s_idx, Timestamp::start);
            if(clipEnd > start) times.setObservingTime(s_idx, clipEnd, Timestamp::end);
        }
    }

    // update pointing vectors according to the new scan times
    for(Scan& scan : scans) {
        for(unsigned int s_id : scan.getStationIds()) {
            refreshPointingVector(scan, s_id, Timestamp::start);
            refreshPointingVector(scan, s_id, Timestamp::end);
        }
        scan.updateObservingTime();
    }
}


std::string ModelBase::dump(size_t t0, size_t tf) const noexcept {
    std::ostringstream output;
    std::map<unsigned long, char> qId;
    for(const unsigned long q : sol_.sourceMask_) {
        qId.insert(std::make_pair(q, static_cast<char>(qId.size() + '!')));
    }
    output << "[" << t0 << ", " << tf << "]" << std::endl;
    for(const Station& s : sol_.network_.getStations()) {
        output << s.getName() << std::endl;
        for(size_t t = 0; t < sol_.blockCount_; ++t) {
            for(const auto q : sol_.getSources()) {
                if(auto sol = sol_.getSol(Solution::Key::StaActive(&sol_, q, s, t))) {
                    if(*sol) {
                        output << qId.at(q->getId());
                        goto next;
                    }
                }
            }
            output << " ";
next:;
        }
        output << std::endl;
    }
    return output.str();
}
#endif // WITH_GUROBI
}
