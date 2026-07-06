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

#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#ifdef WITH_GUROBI
#include "gurobi_c++.h"
#endif // WITH_GUROBI

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
    network_(network), sourceList_(sourceList), sourceMask_(sourceMask), modes_(modes),
    blockLength_(blockLength), 
    blockCount_(TimeSystem::duration / blockLength),
    windowLength_(windowLength), 
    windowBlockCount_((windowLength + blockLength - 1) / blockLength),
    // TODO: This may need a +1
    windowCount_(((TimeSystem::duration / blockLength) - 3) / ((windowLength + blockLength - 1) / blockLength - 2) + 1) {
    coverage_ = std::make_unique<ModelCoverage13>();
#ifdef WITH_GUROBI
    initGurobi(env_);

#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "time segments = " << blockCount_;
    BOOST_LOG_TRIVIAL( info ) << "scan length [sec] = " << blockLength_;
    BOOST_LOG_TRIVIAL( info ) << "optimization steps = " << windowCount_;
    BOOST_LOG_TRIVIAL( info ) << "optimization window length [sec] = " << windowLength_;
    BOOST_LOG_TRIVIAL( info ) << "time segments per window = " << windowBlockCount_;
#else
    std::cout << "[info] time segments = " << blockCount_;
    std::cout << "[info] scan length [sec] = " << blockLength_;
    std::cout << "[info] optimization steps = " << windowBlockCount_;
    std::cout << "[info] optimization window length [sec] = " << windowLength_;
    std::cout << "[info] time segments per window = " << windowBlockCount_;
#endif  

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
    for(size_t t : ModelBase::getBlocks(0, blockCount_)) {
        for(const auto q : ModelBase::getSources()) {
            for(const Station& s : ModelBase::getStations(t, q)) {
                ModelBase::addSol(ModelKey::StaActive(this, q, s, t));
            }
        }
    }

    // BlnActive
    for(size_t t : ModelBase::getBlocks(0, blockCount_)) {
        for(const auto q : ModelBase::getSources()) {
            for(const Baseline& b : ModelBase::getBaselines(t, q)) {
                ModelBase::addSol(ModelKey::BlnActive(this, q, b, t));
            }
        }
    }
#endif // WITH_GUROBI
}

bool ModelBase::optimize(void) {
#ifdef WITH_GUROBI
    for (size_t i = 0; i < windowCount_; ++i) {
        size_t t0 = i * (windowBlockCount_ - 2);
        size_t tf = std::min(t0 + windowBlockCount_, blockCount_);

        // reinitialize the model
        initModel(env_, model_);

        // clear the variable map
        var_.clear();

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Optimizing between " << t0 * blockLength_ << " and " << tf * blockLength_;
#else
        std::cout << "[info] Optimizing between " << t0 * blockLength_ << " and " << tf * blockLength_;
#endif

        // StaActive
        size_t count = 0;
        for(size_t t : ModelBase::getBlocks(t0, tf)) {
            for(const auto q : ModelBase::getSources()) {
                for(Station& s : ModelBase::getStations(t, q)) {
                    // create variable
                    ModelBase::ModelKey key = ModelKey::StaActive(this, q, s, t);
                    GRBVar& var = ModelBase::addVar(key, 0.0, 1.0, 0.0, GRB_BINARY);
                    var.set(GRB_DoubleAttr_Start, 0.0);
                    if(auto sol = ModelBase::getSol(key)) {
                        if(*sol) var.set(GRB_DoubleAttr_Start, 1.0);
                    }
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
        for(size_t t : ModelBase::getBlocks(t0, tf)) {
            for(const auto q : ModelBase::getSources()) {
                for(const Baseline& b : ModelBase::getBaselines(t, q)) {
                    ModelBase::ModelKey key = ModelKey::BlnActive(this, q, b, t);
                    GRBVar& var = ModelBase::addVar(key, 0.0, 1.0, 0.0, GRB_BINARY);
                    var.set(GRB_DoubleAttr_Start, 0.0);
                    if(auto sol = ModelBase::getSol(key)) {
                        if(*sol) var.set(GRB_DoubleAttr_Start, 1.0);
                    }
                    count++;
                }
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " baseline activity variables to model";
#else
        std::cout << "[info] Added " << count << " BlnActive variables to model";
#endif

        // StaCoverage
        count = 0;
        for(const Station& s : ModelBase::getStations()) {
            for(std::size_t c = 0; c < coverage_->cellCount(); ++c) {
                ModelBase::addVar(ModelKey::StaCoverage(this, s, c), 0.0, 1.0, 0.0, GRB_BINARY);
                count++;
            } 
        }

        // populate StaConverage variables
        for(Station& s : ModelBase::getStations()) {
            for(size_t c = 0; c < coverage_->cellCount(); ++c) {
                for(const auto q : ModelBase::getSources()) {
                    for(size_t t : ModelBase::getBlocks(t0, tf, q, s)) {
                        if(coverage_->calculateCell(this, t, q, s) != c) continue;
                        if(*getSol(ModelKey::StaActive(this, q, s, t))) {
                            auto var = *getVar(ModelKey::StaCoverage(this, s, c));
                            var.set(GRB_DoubleAttr_Start, 1.0);
                            goto next_c;
                        }
                    }
                }
next_c:
                (void) nullptr;
            }
        }

#ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Added " << count << " sky coverage variables to model";
#else
        std::cout << "[info] Added " << count << " StaCoverage variables to model";
#endif

        // update the model to make sure variables are accessible
        model_->update();
        
        // optimize the window
        this->prepare(t0, tf);

        // optimize
        model_->optimize();

        int status = model_->get(GRB_IntAttr_Status);
        if(status == GRB_INFEASIBLE) {
            model_->computeIIS();
            model_->write("/tmp/iis.ilp");
        }

        // error checking
        if(status != GRB_OPTIMAL && status != GRB_SUBOPTIMAL) {
    #ifdef VIESCHEDPP_LOG
            BOOST_LOG_TRIVIAL( info ) << "No optimal solution found between " << t0 * blockLength_ << " and " << tf * blockLength_;
    #else
            std::cout << "[info] No optimal solution found between " << t0 * blockLength_ << " and " << tf * blockLength_;
    #endif
            return false;
        }

    #ifdef VIESCHEDPP_LOG
        BOOST_LOG_TRIVIAL( info ) << "Completed optimization between " << t0 * blockLength_ << " and " << tf * blockLength_;
    #else
        std::cout << "[info] Completed optimization between " << t0 * blockLength_ << " and " << tf * blockLength_;
    #endif

        // copy results back into solution
        for(size_t t = t0; t < tf; ++t) {
            for(const auto q : ModelBase::getSources()) {
                for(Station& s : ModelBase::getStations(t, q)) {
                    ModelBase::ModelKey key = ModelKey::StaActive(this, q, s, t);
                    auto sol = ModelBase::getSol(key);
                    auto var = ModelBase::getVar(key);
                    (*sol) = (var->get(GRB_DoubleAttr_X) > 0.5);
                }
                for(const Baseline& b : ModelBase::getBaselines(t, q)) {
                    ModelBase::ModelKey key = ModelKey::BlnActive(this, q, b, t);
                    auto sol = ModelBase::getSol(key);
                    auto var = ModelBase::getVar(key);
                    (*sol) = (var->get(GRB_DoubleAttr_X) > 0.5);
                }
            }
        }
    }
    
    return true;
#else // WITH_GUROBI
    return false;
#endif // WITH_GUROBI   
}

std::vector<Scan> ModelBase::optimize(std::vector<Scan>& scans) {
#ifdef WITH_GUROBI
    ModelBase::loadScans(scans);

    std::string output_greedy = ModelBase::dump();

    if(!ModelBase::optimize()) return {};

    std::cout << output_greedy;
    std::cout << "===========================================" << std::endl;
    std::cout << ModelBase::dump();

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

std::size_t ModelCoverage13::calculateCell(const ModelBase* model, size_t t, 
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
bool ModelBase::checkStationVisibility(size_t t, 
    std::shared_ptr<const VieVS::AbstractSource> q, Station& s) const noexcept {
    // make sure source is visible at this time
    PointingVector pv0(s.getId(), q->getId());
    pv0.setTime(t * blockLength_);
    s.calcAzEl_rigorous( q, pv0);
    if(!s.isVisible(pv0, q->getPARA().minElevation)) return false;
    PointingVector pvf(s.getId(), q->getId());
    pvf.setTime((t + 1) * blockLength_ - 1);
    s.calcAzEl_rigorous( q, pvf);
    return s.isVisible(pvf, q->getPARA().minElevation);
}

unsigned int ModelBase::calculateMinObsExact(unsigned int t,
    const std::shared_ptr<const AbstractSource>& q,
    const Baseline& b,
    const std::shared_ptr<const Mode> &mode) {
    boost::optional<unsigned int> fixedScanDuration = q->getPARA().fixedScanDuration;
    if(auto fixedScanDuration = q->getPARA().fixedScanDuration) {
        return *fixedScanDuration;
    }

    unsigned long staid1 = b.getStaid1();
    Station &sta1 = network_.refStation(staid1);
    unsigned long staid2 = b.getStaid2();
    Station &sta2 = network_.refStation( staid2 );

    // calculate greenwhich meridian sedirial time
    double date1 = 2400000.5;
    double date2 = TimeSystem::mjdStart + static_cast<double>( t ) / 86400.0;
    double gmst = iauGmst82( date1, date2 );

    unsigned int maxDuration = 0;
    for(auto &band : mode->getAllBands()) {
        double SEFD_src;
        if ( q->hasFluxInformation( band ) ) {
            // calculate observed flux density for each band
            SEFD_src = q->observedFlux( band, t, gmst, network_.getDxyz( staid1, staid2 ) );
        } else if ( ObservingMode::sourceBackup[band] == ObservingMode::Backup::internalModel ) {
            // calculate observed flux density based on model
            double wavelength = ObservingMode::wavelengths[band];
            SEFD_src = q->observedFlux_model( wavelength, t, gmst, network_.getDxyz( staid1, staid2 ) );
        } else {
            SEFD_src = 1e-3;
        }

        if ( SEFD_src == 0 ) {
            SEFD_src = 1e-3;
        }

        PointingVector pv1(staid1, q->getId());
        PointingVector pv2(staid2, q->getId());

        pv1.setTime(t);
        pv2.setTime(t);
        sta1.calcAzEl_rigorous(q, pv1);
        sta2.calcAzEl_rigorous(q, pv2);

        // calculate system equivalent flux density for each station
        double el1 = pv1.getEl();
        double SEFD_sta1 = sta1.getEquip().getSEFD( band, el1 );
        double el2 = pv2.getEl();
        double SEFD_sta2 = sta2.getEquip().getSEFD( band, el2 );

        // get minimum required SNR for each station, baseline and source
        double minSNR_sta1 = sta1.getPARA().minSNR.at( band );
        double minSNR_sta2 = sta2.getPARA().minSNR.at( band );
        double minSNR_bl = b.getParameters().minSNR.at( band );
        double minSNR_src = q->getPARA().minSNR.at( band );

        // maximum required minSNR
        double maxminSNR = std::max( { minSNR_src, minSNR_bl, minSNR_sta1, minSNR_sta2 } );

        // get maximum correlator synchronization time for
        double maxCorSynch1 = sta1.getPARA().midob;
        double maxCorSynch2 = sta2.getPARA().midob;
        double maxCorSynch = std::max( { maxCorSynch1, maxCorSynch2 } );

        // calc required baseline scan duration
        double efficiency = mode->efficiency( sta1.getId(), sta2.getId() );
        double anum = ( maxminSNR / ( SEFD_src * efficiency ) );
        double anu1 = SEFD_sta1 * SEFD_sta2;
        double anu2 = mode->recordingRate( staid1, staid2, band );
        if ( anu2 == 0 ) {
            return std::numeric_limits<unsigned int>::max();
        }
        double new_duration = anum * anum * anu1 / anu2 + maxCorSynch;
        new_duration = ceil( new_duration );
        auto new_duration_uint = static_cast<unsigned int>( new_duration );

        // check if required baseline scan duration is within min and max scan times of baselines
        unsigned int minScanBl = b.getParameters().minScan;
        if ( new_duration_uint < minScanBl ) {
            new_duration_uint = minScanBl;
        }
        unsigned int maxScanBl = b.getParameters().maxScan;
        if(new_duration_uint > maxDuration) {
            maxDuration = new_duration_uint;
        }
    }
    return maxDuration;
}

size_t ModelBase::calculateMinObs(size_t t,
    const std::shared_ptr<const AbstractSource>& q,
    const Baseline& b,
    const std::shared_ptr<const Mode> &mode) {
    return (ModelBase::calculateMinObsExact(t * blockLength_, q, b, mode) + blockLength_ - 1) / blockLength_;
}

unsigned int ModelBase::calculateSlewTimeExact(Station& s, 
    const std::shared_ptr<const AbstractSource> q1, 
    const std::shared_ptr<const AbstractSource> q2,
    unsigned int t1, unsigned int t2) const noexcept {
    if(q1->getId() == q2->getId()) return 0;
    
    PointingVector pv1(s.getId(), q1->getId());
    PointingVector pv2(s.getId(), q2->getId());

    pv1.setTime(t1);
    pv2.setTime(t2);
    s.calcAzEl_rigorous(q1, pv1);
    s.calcAzEl_rigorous(q2, pv2);

    PointingVector tempVec(pv2);
    if(!s.isVisible(tempVec, q2->getPARA().minElevation)) 
        return std::numeric_limits<unsigned int>::max();

    unsigned int t_slew = s.getAntenna().slewTime(pv1, pv2);
    unsigned int t_const = s.getPARA().systemDelay + s.getPARA().preob;

    return t_slew + t_const;
}

size_t ModelBase::calculateSlewTime(Station& s, 
    const std::shared_ptr<const AbstractSource> q1, 
    const std::shared_ptr<const AbstractSource> q2,
    size_t t1, size_t t2) const noexcept {
    unsigned int t = ModelBase::calculateSlewTimeExact(s, q1, q2, t1 * blockLength_ + blockLength_, t2 * blockLength_);
    return (t + blockLength_ - 1) / blockLength_;
}

std::vector<size_t> ModelBase::getBlocks(size_t t0, size_t tf) const noexcept {
    std::vector<size_t> blocks(tf - t0);
    std::iota(blocks.begin(), blocks.end(), t0);
    return blocks;
}

std::vector<size_t> ModelBase::getBlocks(size_t t0, size_t tf, const std::shared_ptr<const AbstractSource>& q, Station& s) const noexcept {
    std::vector<size_t> blocks = ModelBase::getBlocks(t0, tf);
    std::vector<size_t> blocksMasked;
    std::copy_if(blocks.begin(), blocks.end(), std::back_inserter(blocksMasked), 
        [this, q, &s](size_t obj) { return this->checkStationVisibility(obj, q, s); });
    return blocksMasked;
}

std::vector<size_t> ModelBase::getBlocks(size_t t0, size_t tf, const std::shared_ptr<const AbstractSource>& q, const Baseline& b) noexcept {
    std::vector<size_t> blocks = ModelBase::getBlocks(t0, tf);
    std::vector<size_t> blocksMasked;
    std::copy_if(blocks.begin(), blocks.end(), std::back_inserter(blocksMasked), 
        [this, q, &b](size_t obj) { 
            Station& s1 = this->network_.refStation(b.getStaid1());
            Station& s2 = this->network_.refStation(b.getStaid2());
            return this->checkStationVisibility(obj, q, s1) && this->checkStationVisibility(obj, q, s2); 
        });
    return blocksMasked;
}

std::vector<std::reference_wrapper<Station>> ModelBase::getStations() noexcept {
    std::vector<Station>& stations = network_.refStations();
    std::vector<std::reference_wrapper<Station>> refs;
    refs.reserve(stations.size());
    std::transform(stations.begin(), stations.end(), std::back_inserter(refs), 
        [](Station& obj) { return std::ref(obj); });
    return refs;
}

std::vector<std::reference_wrapper<Station>> ModelBase::getStations(size_t t, const std::shared_ptr<const AbstractSource>& q) noexcept {
    std::vector<std::reference_wrapper<Station>> refs = ModelBase::getStations();
    std::vector<std::reference_wrapper<Station>> refsMasked;
    std::copy_if(refs.begin(), refs.end(), std::back_inserter(refsMasked), 
        [this, t, q](Station& obj) { return this->checkStationVisibility(t, q, obj); });
    return refsMasked;
}

std::vector<std::reference_wrapper<const Baseline>> ModelBase::getBaselines() const noexcept {
    const std::vector<Baseline>& baselines = network_.getBaselines();
    std::vector<std::reference_wrapper<const Baseline>> refs;
    refs.reserve(baselines.size());
    std::transform(baselines.begin(), baselines.end(), std::back_inserter(refs), 
        [](const Baseline& obj) { return std::cref(obj); });
    return refs;
}

std::vector<std::reference_wrapper<const Baseline>> ModelBase::getBaselines(size_t t, const std::shared_ptr<const AbstractSource>& q) noexcept {
    std::vector<std::reference_wrapper<const Baseline>> refs = ModelBase::getBaselines();
    std::vector<std::reference_wrapper<const Baseline>> refsMasked;
    std::copy_if(refs.begin(), refs.end(), std::back_inserter(refsMasked), 
        [this, t, q](const Baseline& obj) { 
            Station& s1 = this->network_.refStation(obj.getStaid1());
            Station& s2 = this->network_.refStation(obj.getStaid2());
            return this->checkStationVisibility(t, q, s1) && this->checkStationVisibility(t, q, s2); 
        });
    return refsMasked;
}

std::vector<std::shared_ptr<const AbstractSource>> ModelBase::getSources() const noexcept {
    std::vector<std::shared_ptr<const AbstractSource>> sources = sourceList_.getSources();
    std::vector<std::shared_ptr<const AbstractSource>> sourcesMasked;
    std::copy_if(sources.begin(), sources.end(), std::back_inserter(sourcesMasked), 
        [this](const std::shared_ptr<const AbstractSource>& obj) { 
            return this->sourceMask_.count(obj->getId()) > 0; 
        });
    return sourcesMasked;
}

std::vector<std::shared_ptr<const AbstractSource>> ModelBase::getSources(size_t t, Station& s) const noexcept {
    std::vector<std::shared_ptr<const AbstractSource>> refs = ModelBase::getSources();
    std::vector<std::shared_ptr<const AbstractSource>> refsMasked;
    std::copy_if(refs.begin(), refs.end(), std::back_inserter(refsMasked), 
        [this, t, &s](const std::shared_ptr<const AbstractSource>& obj) { 
            return this->checkStationVisibility(t, obj, s);
        });
    return refsMasked;
}

std::vector<std::shared_ptr<const AbstractSource>> ModelBase::getSources(size_t t, const Baseline& b) noexcept {
    std::vector<std::shared_ptr<const AbstractSource>> refs = ModelBase::getSources();
    std::vector<std::shared_ptr<const AbstractSource>> refsMasked;
    std::copy_if(refs.begin(), refs.end(), std::back_inserter(refsMasked), 
        [this, t, &b](const std::shared_ptr<const AbstractSource>& obj) { 
            Station& s1 = this->network_.refStation(b.getStaid1());
            Station& s2 = this->network_.refStation(b.getStaid2());
            return this->checkStationVisibility(t, obj, s1) && this->checkStationVisibility(t, obj, s2); 
        });
    return refsMasked;
}

#ifdef WITH_GUROBI
boost::optional<GRBVar> ModelBase::getVar(const ModelKey& key) const noexcept {
    try {
        return var_.at(key);
    } catch(...) {
        return boost::none;
    }
}

bool* ModelBase::getSol(const ModelKey& key) noexcept {
    auto it = sol_.find(key);
    if(it == sol_.end()) return nullptr;
    return &it->second;
}

const bool* ModelBase::getSol(const ModelKey& key) const noexcept {
    auto it = sol_.find(key);
    if(it == sol_.end()) return nullptr;
    return &it->second;
}

GRBVar& ModelBase::addVar(const ModelKey& key, double lb, double ub, double obj, char vtype) {
    auto ret = var_.insert(std::make_pair(key, model_->addVar(lb, ub, obj, vtype, key.getName())));
    if(!ret.second) throw UNREACHABLE;
    return ret.first->second;
}

bool& ModelBase::addSol(const ModelKey& key) {
    auto ret = sol_.insert(std::make_pair(key, false));
    if(!ret.second) throw UNREACHABLE;
    return ret.first->second;
}

bool ModelBase::ModelKey::operator<(const ModelBase::ModelKey& other) const {
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

ModelBase::ModelKey ModelBase::ModelKey::StaActive(const ModelBase* model, 
    std::shared_ptr<const VieVS::AbstractSource> const q, 
    const Station& s, size_t t) {
    ModelKey key{};
    key.type = ModelKey::ModelKeyType::sta_active;
    key.key.sta_active.q = model->src2idx_.at(q->getId());
    key.key.sta_active.s = model->sta2idx_.at(s.getId());
    key.key.sta_active.t = t;
    key.name += "sta_active_" + s.getName() + "_" + q->getName() + "_" + std::to_string(t);
    return key;
}

ModelBase::ModelKey ModelBase::ModelKey::BlnActive(const ModelBase* model, 
    std::shared_ptr<const VieVS::AbstractSource> const q, 
    const Baseline& b, size_t t) {
    ModelKey key{};
    key.type = ModelKey::ModelKeyType::bln_active;
    key.key.bln_active.q = model->src2idx_.at(q->getId());
    key.key.bln_active.b = model->bln2idx_.at(b.getId());
    key.key.bln_active.t = t;
    key.name += "bln_active__" + model->network_.getStation(b.getStaid1()).getName() + "_" + model->network_.getStation(b.getStaid2()).getName() + "_" + q->getName() + "_" + std::to_string(t);
    return key;
}

ModelBase::ModelKey ModelBase::ModelKey::StaCoverage(const ModelBase* model, const Station& s, size_t c) {
    ModelKey key{};
    key.type = ModelKey::ModelKeyType::sta_coverage;
    key.key.sta_coverage.s = model->sta2idx_.at(s.getId());
    key.key.sta_coverage.c = c;
    key.name += "sta_coverage_" + s.getName() + "_" + std::to_string(c);
    return key;
}

void ModelBase::loadScans(const std::vector<Scan>& scans) {
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
        std::shared_ptr<const VieVS::AbstractSource> const q = sourceList_.getSource(scan.getSourceId());
        if(sourceMask_.count(q->getId()) == 0) continue;
        const ScanTimes& scanTimes = scan.getTimes();
        
        // populate BlnActive variables
        for(const Observation& obs : scan.getObservations()) {
            const Baseline& b = network_.getBaseline(obs.getBlid());
            const Station& s1 = network_.getStation(b.getStaid1());
            const Station& s2 = network_.getStation(b.getStaid2());
            // observation start blocks
            size_t t10 = (scanTimes.getObservingTime(s1.getId()) + blockLength_ - 1) / blockLength_ + 1;
            size_t t20 = (scanTimes.getObservingTime(s2.getId()) + blockLength_ - 1) / blockLength_ + 1;
            // the number of blocks each station is observing
            size_t t1f = t10 + scanTimes.getObservingDuration(s1.getId()) / blockLength_ - 1;
            size_t t2f = t20 + scanTimes.getObservingDuration(s2.getId()) / blockLength_ - 1;
            t1f = std::min(t1f, blockCount_);
            t2f = std::min(t2f, blockCount_);
            
            bool started = false;
            size_t t_start = std::max(t10, t20);
            size_t t_end = std::min(t1f, t2f);
            size_t t_delayed = 0;
            size_t t_premature = 0;
            for(size_t t = t_start; t < t_end; ++t) {
                auto solS2 = getSol(ModelKey::StaActive(this, q, s1, t));
                auto solS1 = getSol(ModelKey::StaActive(this, q, s2, t));
                auto solBL = getSol(ModelKey::BlnActive(this, q, b, t));
                if(!solS1 || !solS2 || !solBL) {
                    if(started) {
                        t_premature = t_end - t;
                        break;
                    }
                    continue;
                }
                bool s1_available = true;
                bool s2_available = true;
                for(const auto q2 : ModelBase::getSources()) {
                    if(q->getId() == q2->getId()) continue;
                    if(auto sol1 = getSol(ModelKey::StaActive(this, q2, s1, t))) {
                        if(*sol1) s1_available = false;
                    }
                    if(auto sol2 = getSol(ModelKey::StaActive(this, q2, s2, t))) {
                        if(*sol2) s2_available = false;
                    }
                }
                if(s1_available && s2_available) {
                    (*solS1) = true;
                    (*solS2) = true;
                    (*solBL) = true;
                    started = true;
                }
                
                if(!started) {
                    t_delayed = t - t_start;
                }
            }
            if(t_delayed > 0 || t_premature > 0) {
                // TODO: proper logging
                std::cout << "t_delayed: " << t_delayed * blockLength_ << "s, t_premature: " << t_premature * blockLength_ << std::endl;
            }
        }
    }

    // disable observations that don't respect the minNumberOfSites parameter due to discretization
    for(const auto q : ModelBase::getSources()) {
        unsigned int minNumberOfSites = q->getPARA().minNumberOfSites;
        for(size_t t : ModelBase::getBlocks(0, blockCount_)) {
            size_t active = 0;
            for(const Station& s : ModelBase::getStations(t, q)) {
                if(*getSol(ModelKey::StaActive(this, q, s, t))) ++active;
            }
            if(active < minNumberOfSites) {
                for(const Station& s : ModelBase::getStations(t, q)) {
                    *getSol(ModelKey::StaActive(this, q, s, t)) = false;
                }
                for(const Baseline& b : ModelBase::getBaselines(t, q)) {
                    *getSol(ModelKey::BlnActive(this, q, b, t)) = false;
                }
            }
        } 
    }

#ifdef VIESCHEDPP_LOG
    BOOST_LOG_TRIVIAL( info ) << "Loaded preliminary result into ILP model";
#else
    std::cout << "[info] Loaded preliminary result into ILP model";
#endif
}

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
        PointingVector pv(s.getId(), qId);
        pv.setTime(t * model->blockLength_);
        s.calcAzEl_rigorous(q, pv);

        sData.insert(std::make_pair(s.getId(), std::make_pair(pv, t + 1)));
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

    size_t blockLength = model->blockLength_;
    // std::transform(sData.begin(), sData.end(), std::back_inserter(endOfLastScan),
    //     [blockLength](const Entry& entry) { return entry.second.second * blockLength; });

    pointingVectorsEnd.reserve(pointingVectors.size());
    for(size_t i = 0; i < pointingVectors.size(); ++i) {
        const PointingVector& pv0 = pointingVectors[i];

        PointingVector pve(pv0.getStaid(), qId);
        pve.setTime(sData.at(pv0.getStaid()).second * blockLength);

        Station& s = model->network_.refStation(pv0.getStaid());
        std::shared_ptr<const VieVS::AbstractSource> const q = model->sourceList_.getSource(qId);
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

        const Station& s = model->network_.getStation(pv0.getStaid());
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
            unsigned long blid = model->network_.getBaseline(staids).getId();

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
    pv.setTime(t * model_->blockLength_);
    s.calcAzEl_rigorous(q, pv);

    ModelBase::ScanBuilder scan;
    scan.qId = q->getId();
    scan.sData.insert(std::make_pair(s.getId(), std::make_pair(pv, t + 1)));
    
    if(t == 0) {

    }

    if(scanContinuation) {
        for(const auto& sEntry : scanContinuation->sData) {
            if(sEntry.first == s.getId() || sEntry.second.second <= t) continue;
            scan.sContinuation.insert(sEntry.first);

            const auto& sData = sEntry.second;

            PointingVector pv(sEntry.first, q->getId());
            pv.setTime(t * model_->blockLength_);
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
        const Station& s = model_->network_.getStation(sEntryInner.first);
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
        const Station& s = model_->network_.getStation(sEntryInner.first);
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
    for(size_t t = 0; t <= blockCount_; ++t) {
        for(Station& s : network_.refStations()) {
            for(const auto q : sourceList_.getSources()) {
                if(auto sol = getSol(ModelKey::StaActive(this, q, s, t))) {
                    if(*sol) scansActive.append(q, s, t);
                }
            }
        }

        scansActive.updateScans(scans, t);
    }

    return scans;
}

std::string ModelBase::dump() const noexcept {
    std::ostringstream output;
    std::map<unsigned long, char> qId;
    for(const unsigned long q : sourceMask_) {
        qId.insert(std::make_pair(q, static_cast<char>(qId.size() + 65)));
    }
    for(const Station& s : network_.getStations()) {
        output << s.getName() << std::endl;
        for(size_t t = 0; t < blockCount_; ++t) {
            for(const auto q : sourceList_.getSources()) {
                if(auto sol = getSol(ModelKey::StaActive(this, q, s, t))) {
                    if(*sol) {
                        output << qId.at(q->getId());
                        goto next_t;
                    }
                }
            }
            output << " ";
next_t:
            (void) nullptr;
        }
        output << std::endl;
    }
    return output.str();
}
#endif // WITH_GUROBI
}
