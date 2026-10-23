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
            model_->reset(0);
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
            // std::map<int, double> val;
            // for(const auto q : sol_.getSources()) {
            //     for(size_t t : sol_.getBlocks(t0, tf)) {
            //         for(const Station& s : sol_.getStations(t, q)) {
            //             auto var = getVar(Solution::Key::StaActive(&sol_, q, s, t));
            //             val.emplace(var->index(), var->get(GRB_DoubleAttr_X));
            //         }
            //         for(const Baseline& b : sol_.getBaselines(t, q)) {
            //             auto var = getVar(Solution::Key::BlnActive(&sol_, q, b, t));
            //             val.emplace(var->index(), var->get(GRB_DoubleAttr_X));
            //         }
            //     }
            // }
            // for(const Station& s : sol_.getStations()) {
            //     for(size_t c = 0; c < coverage_->cellCount(); ++c) {
            //         auto var = getVar(Solution::Key::StaCoverage(&sol_, s, c));
            //         val.emplace(var->index(), var->get(GRB_DoubleAttr_X));
            //     }
            // }
            // model_->reset(1);
            // for(const auto q : sol_.getSources()) {
            //     for(size_t t : sol_.getBlocks(t0, tf)) {
            //         for(const Station& s : sol_.getStations(t, q)) {
            //             auto var = getVar(Solution::Key::StaActive(&sol_, q, s, t));
            //             var->set(GRB_DoubleAttr_X, val[var->index()]);
            //         }
            //         for(const Baseline& b : sol_.getBaselines(t, q)) {
            //             auto var = getVar(Solution::Key::BlnActive(&sol_, q, b, t));
            //             var->set(GRB_DoubleAttr_X, val[var->index()]);
            //         }
            //     }
            // }
            // for(const Station& s : sol_.getStations()) {
            //     for(size_t c = 0; c < coverage_->cellCount(); ++c) {
            //         auto var = getVar(Solution::Key::StaCoverage(&sol_, s, c));
            //         var->set(GRB_DoubleAttr_X, val[var->index()]);
            //     }
            // }
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
                return false;
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
    return ModelBase::readScans();
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
#if 0
while (t10 < sol_.blockCount_) {
    // Find the next block where q is active.
    while (t10 < sol_.blockCount_ &&
           !sol_.getSol(Solution::Key::StaActive(&sol_, q, s1, t10)))
        ++t10;

    if (t10 >= sol_.blockCount_)
        break;

    size_t requiredStart = t10;

    // Check every preceding active source.
    for (size_t t = t10; t-- > 0;) {
        for (const auto q2 : sol_.getSources()) {
            if (q2->getId() == q->getId())
                continue;

            auto active = sol_.getSol(
                Solution::Key::StaActive(&sol_, q2, s1, t));
          
            if (!active || !*active)
                continue;

            const size_t t_slew = sol_.getSlew(
                Solution::Key::StaActive(&sol_, q2, s1, t + 1),
                Solution::Key::StaActive(&sol_, q, s1, t10));

            if(t10 == 30) {
                // std::cout << s1.getName() << ","<< s2.getName()<< "," << t_slew <<","<<q2->getName()<<","<<t<<","<<q->getName()<<","<<t10<<std::endl;
            }
            if ((s1.getName() == "ZELENCHK" || s1.getAlternativeName() == "ZELENCHK") &&
                (q2->getName() == "1144+402" || q2->getAlternativeName()  == "1144+402") &&
                (q->getName() == "3C371" || q->getAlternativeName() == "3C371") &&
                t == 23 && t10 == 30) {
                std::cout << "WARM START slew=" << t_slew << '\n';
            }

            // constrSlew requires t10 - t > t_slew.
            if (t10 - t <= t_slew) {
                requiredStart = std::max(requiredStart, t + t_slew + 1);
            }
        }
    }

    if (requiredStart == t10)
        break; // No slew conflict.

    t10 = requiredStart;
}
while (t20 < sol_.blockCount_) {
    // Find the next block where q is active.
    while (t20 < sol_.blockCount_ &&
           !sol_.getSol(Solution::Key::StaActive(&sol_, q, s2, t20)))
        ++t20;

    if (t20 >= sol_.blockCount_)
        break;

    size_t requiredStart = t20;

    // Check every preceding active source.
    for (size_t t = t20; t-- > 0;) {
        for (const auto q2 : sol_.getSources()) {
            if (q2->getId() == q->getId())
                continue;

            auto active = sol_.getSol(
                Solution::Key::StaActive(&sol_, q2, s2, t));
             
            if (!active || !*active)
                continue;

            const size_t t_slew = sol_.getSlew(
                Solution::Key::StaActive(&sol_, q2, s2, t + 1),
                Solution::Key::StaActive(&sol_, q, s2, t20));
            if(t20 == 30) {
                // std::cout << s1.getName() <<","<< s2.getName() << "," << t_slew <<","<<q2->getName()<<","<<t<<","<<q->getName()<<","<<t20<<std::endl;
            }
            if ((s2.getName() == "ZELENCHK" || s2.getAlternativeName() == "ZELENCHK") &&
                (q2->getName() == "1144+402" || q2->getAlternativeName()  == "1144+402") &&
                (q->getName() == "3C371" || q->getAlternativeName() == "3C371") &&
                t == 23 && t20 == 30) {
                std::cout << "WARM START slew=" << t_slew << '\n';
            }

            // constrSlew requires t20 - t > t_slew.
            if (t20 - t <= t_slew) {
                requiredStart = std::max(requiredStart, t + t_slew + 1);
            }
        }
    }

    if (requiredStart == t20)
        break; // No slew conflict.

    t20 = requiredStart;
}
#else
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
#endif
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

#if 0
boost::optional<bool> ModelBase::ScanBuilder::append(const ModelBase* model, std::shared_ptr<const VieVS::AbstractSource> const q, 
    Station& s, size_t t) noexcept {
    if(qId != q->getId()) return boost::none;
    try {
        size_t& end = sData.at(s.getId()).second;
        if(end == t) {
            ++end;
            return true;
        }
        return false;
    } catch(...) {
        auto pv = model->sol_.getPointingVector(Solution::Key::StaActive(&model->sol_, q, s, t));
        sData.insert(std::make_pair(s.getId(), std::make_pair(*pv, t + 1)));
        if(t == 0) sContinuation.insert(s.getId());
        return true;
    }
}

Scan ModelBase::ScanBuilder::finish(const ModelBase* model, 
    const std::vector<unsigned int>& slewTime, 
    std::vector<unsigned int>& endOfLastScan) const noexcept {
    std::vector<PointingVector> pointingVectors;
    std::vector<PointingVector> pointingVectorsEnd;
    // std::vector<unsigned int> endOfLastScan;

    typedef std::pair<const unsigned long, std::pair<PointingVector, size_t>> Entry;
    std::transform(sData.begin(), sData.end(), std::back_inserter(pointingVectors),
        [](const Entry& entry) { return entry.second.first; });

    size_t blockLength = model->sol_.blockLength_;
    // std::transform(sData.begin(), sData.end(), std::back_inserter(endOfLastScan),
    //     [blockLength](const Entry& entry) { return entry.second.second * blockLength; });

    pointingVectorsEnd.reserve(pointingVectors.size());
    for(size_t i = 0; i < pointingVectors.size(); ++i) {
        const PointingVector& pv0 = pointingVectors[i];

        PointingVector pve(pv0.getStaid(), qId);
        pve.setTime(sData.at(pv0.getStaid()).second * blockLength);

        Station& s = model->sol_.network_.refStation(pv0.getStaid());
        std::shared_ptr<const VieVS::AbstractSource> const q = model->sol_.sourceList_.getSource(qId);
        s.calcAzEl_rigorous(q, pve);

        pointingVectorsEnd.push_back(pve);
    }

    std::vector<unsigned int> fieldSystemTime;
    std::vector<unsigned int> preob;
    std::vector<unsigned int> scanStart;
    std::vector<unsigned int> observingTimes;

    for(size_t i = 0; i < pointingVectors.size(); ++i) {
        const PointingVector& pv0 = pointingVectors[i];
        const PointingVector& pve = pointingVectorsEnd[i];

        const Station& s = model->sol_.network_.getStation(pv0.getStaid());
        if(pv0.getTime() == 0 || sContinuation.count(pv0.getStaid()) > 0) {
            fieldSystemTime.push_back(0);
            preob.push_back(0);
        } else {
            fieldSystemTime.push_back(s.getPARA().systemDelay);
            preob.push_back(s.getPARA().preob);
        }
        
        scanStart.push_back(pv0.getTime());
        observingTimes.push_back(pve.getTime());
    }

    std::vector<Observation> obs;
    for(size_t i = 0; i < pointingVectors.size(); ++i) {
        const PointingVector& pvi = pointingVectors[i];

        for(size_t j = i + 1; j < pointingVectors.size(); ++j) {
            const PointingVector& pvj = pointingVectors[j];
            
            const std::pair<unsigned long, unsigned long> staids(pvi.getStaid(), pvj.getStaid());
            unsigned long blid = model->sol_.network_.getBaseline(staids).getId();

            unsigned int startTime = std::max(pvi.getTime(), pvj.getTime());
            unsigned int endTime = std::min(pointingVectorsEnd[i].getTime(), pointingVectorsEnd[j].getTime());
            if(startTime >= endTime) continue;
            unsigned int observingTime = endTime - startTime;

            obs.emplace_back(blid, pvi.getStaid(), pvj.getStaid(), qId, startTime, observingTime);
        }
    }

    Scan scan(pointingVectors, endOfLastScan, Scan::ScanType::standard);
    scan.setPointingVectorsEndtime(pointingVectorsEnd);

    ScanTimes& times = scan.referenceTime();
    times.setEndOfLastScan( endOfLastScan );
    for ( int i = 0; i < slewTime.size(); ++i ) {
        times.addTimes( i, fieldSystemTime[i], slewTime.at( static_cast<unsigned long>( i ) ), 0 );
        times.setObservingTime( i, observingTimes[i], Timestamp::end );
    }
    times.setObservingStarts( scanStart );

    std::vector<bool> excl;
    for(const auto& sEntry : sData) {
        excl.push_back(sContinuation.count(sEntry.first) > 0);
    }

    bool valid = times.setPreobTime( preob, excl );
    if(!valid) {
#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( error ) << "Failed to set ScanTimes";
#else
        std::cout << "[error] Failed to set ScanTimes";
#endif
    }

    scan.setObservations(obs);
    
    return scan;
}

void ModelBase::ActiveScans::append(std::shared_ptr<const VieVS::AbstractSource> const q, 
    Station& s, size_t t) noexcept {
    boost::optional<const ModelBase::ScanBuilder&> scanContinuation = boost::none;
    for(ModelBase::ScanBuilder& scan : scans_) {
        if(boost::optional<bool> ret = scan.append(model_, q, s, t)) {
            if(*ret) return;
            scanContinuation = scan;
        }
        // if(scan.append(model_, q, s, t)) return;
    }

    PointingVector pv(s.getId(), q->getId());
    pv.setTime(t * model_->sol_.blockLength_);
    s.calcAzEl_rigorous(q, pv);

    ModelBase::ScanBuilder scan;
    scan.qId = q->getId();
    scan.sData.insert(std::make_pair(s.getId(), std::make_pair(pv, t + 1)));

    if(scanContinuation) {
        for(const auto& sEntry : scanContinuation->sData) {
            if(sEntry.first == s.getId() || sEntry.second.second <= t) continue;
            scan.sContinuation.insert(sEntry.first);

            const auto& sData = sEntry.second;

            PointingVector pv(sEntry.first, q->getId());
            pv.setTime(t * model_->sol_.blockLength_);
            s.calcAzEl_rigorous(q, pv);

            scan.sData.insert(std::make_pair(sEntry.first, std::make_pair(pv, t + 1)));
        }

        auto it = std::find_if(scans_.begin(), scans_.end(),
            [&](const ModelBase::ScanBuilder& sb) { return sb.qId == scanContinuation->qId; });
        if(it != scans_.end()) {
            scansToAppend_.push_back(std::move(*it));
            scans_.erase(it);
        }
    }

    scans_.emplace_back(scan);
}

void ModelBase::ActiveScans::updateScans(std::vector<Scan>& scans, size_t t) {
    for(auto& scan : scansToAppend_) {
        std::vector<unsigned int> slewTime = ModelBase::ActiveScans::computeSlewTime(scan, scans);
        std::vector<unsigned int> endOfLastScan = ModelBase::ActiveScans::computeEndOfLastScan(scan, scans);
        scans.push_back(scan.finish(model_, slewTime, endOfLastScan));
    }
    scansToAppend_.clear();

    std::vector<size_t> remove;
    for(size_t i = 0, j = 0; i < scans_.size(); ++i, j = 0) {
        for(const auto& sEntry : scans_[i].sData) {
            const auto& sData = sEntry.second;
            if(sData.second > t) ++j;
        }

        if(j < 2) {
            remove.emplace_back(i);

            std::vector<unsigned int> slewTime = ModelBase::ActiveScans::computeSlewTime(scans_[i], scans);
            std::vector<unsigned int> endOfLastScan = ModelBase::ActiveScans::computeEndOfLastScan(scans_[i], scans);
            Scan scan = scans_[i].finish(model_, slewTime, endOfLastScan);
            scans.emplace_back(scan);
        }
    }

    std::reverse(remove.begin(), remove.end());
    for(size_t i : remove) {
        scans_.erase(scans_.begin() + i);
    }
}

std::vector<unsigned int> ModelBase::ActiveScans::computeSlewTime(const ModelBase::ScanBuilder& scan, 
    const std::vector<Scan>& scans) noexcept {
    std::vector<unsigned int> slewTime;
    for(const auto& sEntryInner : scan.sData) {
        const Station& s = model_->sol_.network_.getStation(sEntryInner.first);
        if(scan.sContinuation.count(sEntryInner.first) > 0) {
            slewTime.push_back(0);
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
            slewTime.push_back(slew);
        }
    }

    return slewTime;
}

std::vector<unsigned int> ModelBase::ActiveScans::computeEndOfLastScan(const ModelBase::ScanBuilder& scan, 
    const std::vector<Scan>& scans) noexcept {
    std::vector<unsigned int> endOfLastScan;
    for(const auto& sEntryInner : scan.sData) {
        const Station& s = model_->sol_.network_.getStation(sEntryInner.first);
        if(scan.sContinuation.count(sEntryInner.first) > 0) {
            endOfLastScan.push_back(0);
        } else {
            unsigned int eols = 0;
            for(size_t i = scans.size(); i-- > 0;) {
                if(auto j = scans[i].findIdxOfStationId(sEntryInner.first)) {
                    eols = scans[i].getPointingVector(*j, Timestamp::end).getTime();
                    break;
                }
            }
            
            endOfLastScan.push_back(eols);
        }
    }

    return endOfLastScan;
}

std::vector<Scan> ModelBase::readScans(void) const noexcept {
    std::vector<Scan> scans;

    ModelBase::ActiveScans scansActive(this);
    for(size_t t = 0; t <= sol_.blockCount_; ++t) {
        for(Station& s : sol_.network_.refStations()) {
            for(const auto q : sol_.getSources()) {
                if(auto sol = sol_.getSol(Solution::Key::StaActive(&sol_, q, s, t))) {
                    if(*sol) scansActive.append(q, s, t);
                }
            }
        }

        scansActive.updateScans(scans, t);
    }

    return scans;
}
#endif

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

    for (const auto& entry : intervals) {
        for (const Interval& interval : entry.second) {
            bySource[interval.q].push_back(interval);
        }
    }

    auto printGroup = [&](unsigned int q, const ActiveGroup& group) {
        if (group.stations.empty()) {
            return;
        }

        size_t first = std::numeric_limits<size_t>::max();
        size_t last = 0;

        for (const auto& entry : group.stations) {
            first = std::min(first, entry.second.start);
            last = std::max(last, entry.second.end);
        }

        std::cout << "["
                << sol_.sourceList_.getSource(q)->getName()
                << ", " << first
                << ", " << last
                << "]\n";

        for (const Station& station : sol_.network_.getStations()) {
            auto it = group.stations.find(station.getId());
            if (it == group.stations.end()) {
                continue;
            }

            const Interval& interval = it->second;
            std::cout << station.getName() << (group.continuations.count(station.getId()) > 0 ? "(continuation)" : "(initial)") << '\n';

            for (size_t t = first; t <= last; ++t) {
                std::cout << ((t >= interval.start && t <= interval.end)
                                ? qId.at(q)
                                : ' ');
            }

            std::cout << '\n';
        }
    };

    std::vector<std::pair<unsigned int, ActiveGroup>> groups;
    for (auto& sourceEntry : bySource) {
        const unsigned int q = sourceEntry.first;
        auto& sourceIntervals = sourceEntry.second;

        std::sort(sourceIntervals.begin(), sourceIntervals.end(),
            [](const Interval& a, const Interval& b) {
                if (a.start != b.start) return a.start < b.start;
                if (a.end != b.end) return a.end < b.end;
                return a.s < b.s;
            });

        ActiveGroup current;
        bool haveGroup = false;

        size_t i = 0;
        while (i < sourceIntervals.size()) {
            // Process all station starts at this epoch together.
            const size_t epoch = sourceIntervals[i].start;
            size_t j = i;

            while (j < sourceIntervals.size() &&
                sourceIntervals[j].start == epoch) {
                ++j;
            }

            // Determine whether this epoch introduces a new station
            // into an observation that is already active.
            bool lateJoin = false;

            if (haveGroup) {
                bool observationActive = false;

                for (const auto& entry : current.stations) {
                    if (entry.second.end >= epoch) {
                        observationActive = true;
                        break;
                    }
                }

                if (observationActive) {
                    for (size_t k = i; k < j; ++k) {
                        if (current.stations.find(sourceIntervals[k].s) ==
                            current.stations.end()) {
                            lateJoin = true;
                            break;
                        }
                    }
                }
            }

            if (haveGroup && !lateJoin) {
                bool observationActive = false;

                for (const auto& entry : current.stations) {
                    if (entry.second.end >= epoch) {
                        observationActive = true;
                        break;
                    }
                }

                if (!observationActive) {
                    groups.emplace_back(q, current);
                    current = ActiveGroup{};
                    haveGroup = false;
                }
            }

            if (lateJoin) {
                ActiveGroup continuation;

                // Close the initial group immediately before the join.
                for (auto it = current.stations.begin();
                    it != current.stations.end();) {
                    Interval& oldInterval = it->second;

                    if (oldInterval.end >= epoch) {
                        Interval continued = oldInterval;
                        continued.start = epoch;

                        oldInterval.end = epoch - 1;

                        continuation.stations.emplace(
                            it->first, continued);
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

            // Add the stations starting at this epoch. This is done
            // after splitting so simultaneous joins stay in one group.
            for (size_t k = i; k < j; ++k) {
                const Interval& interval = sourceIntervals[k];
                auto it = current.stations.find(interval.s);

                if (it == current.stations.end()) {
                    current.stations.emplace(interval.s, interval);
                } else {
                    // Merge adjacent or overlapping intervals for
                    // the same station, if encountered.
                    it->second.start = std::min(
                        it->second.start, interval.start);
                    it->second.end = std::max(
                        it->second.end, interval.end);
                }
            }

            haveGroup = true;
            i = j;
        }

        // Do not forget the last group: it may still be active.
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
                s_pv_0.setTime((eols.at(s_id).second + 1) * sol_.blockLength_ + s.getPARA().systemDelay);
                s.calcAzEl_rigorous(q_0, s_pv_0);
                PointingVector s_pv_f{s_id, q_id};
                s_pv_f.setTime(s_interval.start * sol_.blockLength_ - s.getPARA().preob);
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
