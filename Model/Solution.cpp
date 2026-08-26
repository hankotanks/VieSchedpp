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

#include "Solution.h"

// system
#include <limits>

#define UNREACHABLE std::logic_error((boost::format("unreachable: %d") % __LINE__).str())

namespace VieVS {
Solution::Solution(VieVS::Network& network, VieVS::SourceList& sourceList, 
    const std::set<unsigned long>& sourceMask, 
    const std::shared_ptr<const ObservingMode>& modes,
    unsigned int blockLength, unsigned int windowLength) :
    network_(network), sourceList_(sourceList), sourceMask_(sourceMask), modes_(modes),
    blockLength_(blockLength), 
    blockCount_(TimeSystem::duration / blockLength),
    windowLength_(windowLength), 
    windowCount_(TimeSystem::duration / windowLength),
    windowBlockCount_((windowLength + blockLength - 1) / blockLength) {
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

    // populate sta_, bln_, src_
    std::vector<Station>& sta = network_.refStations();
    sta_.reserve(sta.size());
    std::transform(sta.begin(), sta.end(), std::back_inserter(sta_), 
        [](Station& obj) { return std::cref(obj); });
    // bln_
    std::vector<Baseline>& bln = network_.refBaselines();
    bln_.reserve(bln.size());
    std::transform(bln.begin(), bln.end(), std::back_inserter(bln_), 
        [](Baseline& obj) { return std::cref(obj); });
    // src_
    src_.reserve(sourceMask_.size());
    std::transform(sourceMask_.begin(), sourceMask_.end(), std::back_inserter(src_), 
        [this](unsigned long obj) { return this->sourceList_.getSource(obj); });
    // build sta2idx_
    for(const Station& s : sta_) {
        sta2idx_.insert(std::make_pair(s.getId(), sta2idx_.size()));
        idx2sta_.emplace_back(s.getId());
    }
    // build bln2idx_
    for(const Baseline& b : bln_) {
        bln2idx_.insert(std::make_pair(b.getId(), bln2idx_.size()));
        idx2bln_.emplace_back(b.getId());
    }
    // build src2idx_
    for(const auto q : src_) {
        src2idx_.insert(std::make_pair(q->getId(), src2idx_.size()));
        idx2src_.emplace_back(q->getId());
    }
    // populate pvs_
    for(size_t t = 0; t < blockCount_; ++t) {
        for(const auto q : src_) {
            for(Station& s : sta) {
                auto key = Key::StaActive(this, q, s, t);
                PointingVector pv0(s.getId(), q->getId());
                pv0.setTime(t * blockLength_);
                s.calcAzEl_rigorous( q, pv0);
                pvs_.emplace(key, pv0);
            }
        }
    }
    // populate vis_
    for(size_t t = 0; t < blockCount_; ++t) {
        for(const auto q : src_) {
            for(const Station& s : sta_) {
                auto key = Key::StaActive(this, q, s, t);
                bool vis;
                if(!s.isVisible(pvs_.at(key), q->getPARA().minElevation)) {
                    vis = false;
                } else {
                    if(t + 1 < blockCount_) {
                        vis = s.isVisible(pvs_.at(Key::StaActive(this, q, s, t + 1)), q->getPARA().minElevation);
                    } else {
                        vis = true;
                    }
                }
                if(!vis) continue; 
                vis_.emplace(key);
            }
        }
    }
    // populate snr_
    for(size_t t = 0; t < blockCount_; ++t) {
        for(const auto q : src_) {
            for(const Baseline& b : bln_) {
                const auto& s1 = network_.getStation(b.getStaid1());
                Key s1key = Key::StaActive(this, q, s1, t);
                if(vis_.count(s1key) == 0) continue;
                const auto& s2 = network_.getStation(b.getStaid2());
                Key s2key = Key::StaActive(this, q, s2, t);
                if(vis_.count(s2key) == 0) continue;
                size_t minObs = Solution::calculateMinObs(q, b, t);
                size_t maxObs = (std::min(q->getPARA().maxScan, std::min(s1.getPARA().maxScan, s2.getPARA().maxScan)) + blockLength_ - 1) / blockLength_;
                if(minObs > maxObs) continue;
                snr_.emplace(Key::BlnActive(this, q, b, t), minObs);
            }
        }
    }
    // populate solution space
    for(size_t t = 0; t < blockCount_; ++t) {
        for(const auto q : src_) {
            for(const Station& s : sta_) {
                auto key = Key::StaActive(this, q, s, t);
                if(vis_.count(key) == 0) continue;
                sol_.emplace(key, false);
            }
        }
    }

    for(size_t t = 0; t < blockCount_; ++t) {
        for(const auto q : src_) {
            for(const Baseline& b : bln_) {
                auto key = Key::BlnActive(this, q, b, t);
                if(snr_.count(key) == 0) continue;
                sol_.emplace(key, false);
            }
        }
    }
}

bool Solution::Key::operator<(const Solution::Key& other) const {
    if(type != other.type) return type < other.type;
    switch(type) {
        case Key::Type::sta_active:
            if(key.sta_active.q != other.key.sta_active.q) return key.sta_active.q < other.key.sta_active.q;
            if(key.sta_active.s != other.key.sta_active.s) return key.sta_active.s < other.key.sta_active.s;
            return key.sta_active.t < other.key.sta_active.t;
        case Key::Type::bln_active:
            if(key.bln_active.q != other.key.bln_active.q) return key.bln_active.q < other.key.bln_active.q;
            if(key.bln_active.b != other.key.bln_active.b) return key.bln_active.b < other.key.bln_active.b;
            return key.bln_active.t < other.key.bln_active.t;
        case Key::Type::sta_coverage:
            if(key.sta_coverage.s != other.key.sta_coverage.s) return key.sta_coverage.s < other.key.sta_coverage.s;
            return key.sta_coverage.c < other.key.sta_coverage.c;
        default: throw UNREACHABLE;
    }
}

bool Solution::Key::operator==(const Solution::Key& other) const noexcept {
    if(type != other.type) return false;
    switch(type) {
    case sta_active:
        return key.sta_active.q == other.key.sta_active.q &&
            key.sta_active.s == other.key.sta_active.s &&
            key.sta_active.t == other.key.sta_active.t;
    case bln_active:
        return key.bln_active.q == other.key.bln_active.q &&
            key.bln_active.b == other.key.bln_active.b &&
            key.bln_active.t == other.key.bln_active.t;
    case sta_coverage:
        return key.sta_coverage.s == other.key.sta_coverage.s &&
            key.sta_coverage.c == other.key.sta_coverage.c;
    }
    return false;
}

Solution::Key Solution::Key::StaActive(const Solution* model, 
    std::shared_ptr<const VieVS::AbstractSource> const q, 
    const Station& s, size_t t) {
    Key key{};
    key.type = Key::Type::sta_active;
    key.key.sta_active.q = model->src2idx_.at(q->getId());
    key.key.sta_active.s = model->sta2idx_.at(s.getId());
    key.key.sta_active.t = t;
    key.name += "sta_active_" + s.getName() + "_" + q->getName() + "_" + std::to_string(t);
    return key;
}

Solution::Key Solution::Key::BlnActive(const Solution* model, 
    std::shared_ptr<const VieVS::AbstractSource> const q, 
    const Baseline& b, size_t t) {
    Key key{};
    key.type = Key::Type::bln_active;
    key.key.bln_active.q = model->src2idx_.at(q->getId());
    key.key.bln_active.b = model->bln2idx_.at(b.getId());
    key.key.bln_active.t = t;
    key.name += "bln_active__" + model->network_.getStation(b.getStaid1()).getName() + "_" + model->network_.getStation(b.getStaid2()).getName() + "_" + q->getName() + "_" + std::to_string(t);
    return key;
}

Solution::Key Solution::Key::StaCoverage(const Solution* model, const Station& s, size_t c) {
    Key key{};
    key.type = Key::Type::sta_coverage;
    key.key.sta_coverage.s = model->sta2idx_.at(s.getId());
    key.key.sta_coverage.c = c;
    key.name += "sta_coverage_" + s.getName() + "_" + std::to_string(c);
    return key;
}

boost::optional<bool&> Solution::getSol(const Solution::Key& key) noexcept {
    auto it = sol_.find(key);
    return it == sol_.end() ? boost::none : boost::optional<bool&>(it->second);
}

boost::optional<const bool&> Solution::getSol(const Solution::Key& key) const noexcept {
    auto it = sol_.find(key);
    return it == sol_.end() ? boost::none : boost::optional<const bool&>(it->second);
}

size_t Solution::getMinObs(const Key& key) const noexcept {
    assert(key.type == Key::Type::bln_active);
    auto it = snr_.find(key);
    return it == snr_.end() ? std::numeric_limits<size_t>::max() : snr_.at(key);
}

size_t Solution::getSlew(const Key& from, const Key& to) const noexcept {
    assert(from.type == Key::Type::sta_active);
    assert(to.type == Key::Type::sta_active);
    assert(from.key.sta_active.s == to.key.sta_active.s);
    if(from.key.sta_active.q == to.key.sta_active.q) return 0;
    auto it1 = pvs_.find(from);
    if(it1 == pvs_.end()) return std::numeric_limits<size_t>::max();
    const PointingVector& pv1 = it1->second;
    auto it2 = pvs_.find(to);
    if(it2 == pvs_.end()) return std::numeric_limits<size_t>::max();
    const PointingVector& pv2 = it2->second;
    const Station& s = network_.getStation(idx2sta_.at(from.key.sta_active.s));
    unsigned int t_slew = s.getAntenna().slewTime(pv1, pv2);
    unsigned int t_const = s.getPARA().systemDelay + s.getPARA().preob;
    return (t_slew + t_const + blockLength_ - 1) / blockLength_;
}

boost::optional<PointingVector> Solution::getPointingVector(const Key& key) const noexcept {
    assert(key.type == Key::Type::sta_active);
    auto it1 = pvs_.find(key);
    if(it1 == pvs_.end()) return boost::none;
    return boost::optional<PointingVector>(it1->second);
}

std::vector<size_t> Solution::getBlocks(size_t t0, size_t tf) const noexcept {
    tf = std::min(tf, blockCount_);
    if(t0 == tf || t0 > tf) return {};
    std::vector<size_t> blocks(tf - t0);
    std::iota(blocks.begin(), blocks.end(), t0);
    return blocks;
}

std::vector<size_t> Solution::getBlocks(size_t t0, size_t tf, const std::shared_ptr<const AbstractSource>& q, const Station& s) const noexcept {
    std::vector<size_t> blocks = Solution::getBlocks(t0, tf);
    std::vector<size_t> mask;
    std::copy_if(blocks.begin(), blocks.end(), std::back_inserter(mask), 
        [this, q, &s](size_t t) { return this->getSol(Key::StaActive(this, q, s, t)); });
    return mask;
}

std::vector<size_t> Solution::getBlocks(size_t t0, size_t tf, const std::shared_ptr<const AbstractSource>& q, const Baseline& b) const noexcept {
    std::vector<size_t> blocks = Solution::getBlocks(t0, tf);
    std::vector<size_t> mask;
    std::copy_if(blocks.begin(), blocks.end(), std::back_inserter(mask), 
        [this, q, &b](size_t obj) { return this->getSol(Key::BlnActive(this, q, b, obj)); });
    return mask;
}

size_t Solution::getBlocks(unsigned int seconds) const noexcept {
    return (seconds + blockLength_ - 1) / blockLength_;
}

std::vector<std::reference_wrapper<const Station>> Solution::getStations() const noexcept {
    std::vector<std::reference_wrapper<const Station>> mask;
    mask.reserve(sta_.size());
    std::transform(sta_.begin(), sta_.end(), std::back_inserter(mask), 
        [](const Station& obj) { return std::ref(obj); });
    return mask;
}

std::pair<std::reference_wrapper<const Station>, std::reference_wrapper<const Station>> Solution::getStations(const Baseline& b) const noexcept {
    const Station& s1 = network_.getStation(b.getStaid1());
    const Station& s2 = network_.getStation(b.getStaid2());
    return std::make_pair(std::cref(s1), std::cref(s2));
}

std::vector<std::reference_wrapper<const Station>> Solution::getStations(size_t t, const std::shared_ptr<const AbstractSource>& q) const noexcept {
    std::vector<std::reference_wrapper<const Station>> mask;
    mask.reserve(sta_.size());
    std::copy_if(sta_.begin(), sta_.end(), std::back_inserter(mask), 
        [this, t, q](const Station& obj) { return this->getSol(Key::StaActive(this, q, obj, t)); });
    return mask;
}

std::vector<std::reference_wrapper<const Baseline>> Solution::getBaselines() const noexcept {
    std::vector<std::reference_wrapper<const Baseline>> mask;
    mask.reserve(bln_.size());
    std::transform(bln_.begin(), bln_.end(), std::back_inserter(mask), 
        [](const Baseline& obj) { return std::cref(obj); });
    return mask;
}

std::vector<std::reference_wrapper<const Baseline>> Solution::getBaselines(size_t t, const std::shared_ptr<const AbstractSource>& q) const noexcept {
    std::vector<std::reference_wrapper<const Baseline>> mask;
    mask.reserve(bln_.size());
    std::copy_if(bln_.begin(), bln_.end(), std::back_inserter(mask), 
        [this, t, q](const Baseline& obj) { return this->getSol(Key::BlnActive(this, q, obj, t)); });
    return mask;
}

std::vector<std::shared_ptr<const AbstractSource>> Solution::getSources() const noexcept {
    std::vector<std::shared_ptr<const AbstractSource>> mask;
    std::copy_if(src_.begin(), src_.end(), std::back_inserter(mask), 
        [this](const std::shared_ptr<const AbstractSource>& obj) { return this->sourceMask_.count(obj->getId()) > 0; });
    return mask;
}

std::vector<std::shared_ptr<const AbstractSource>> Solution::getSources(size_t t, const Station& s) const noexcept {
    std::vector<std::shared_ptr<const AbstractSource>> refsMasked;
    std::copy_if(src_.begin(), src_.end(), std::back_inserter(refsMasked), 
        [this, t, &s](const std::shared_ptr<const AbstractSource>& obj) { return this->getSol(Key::StaActive(this, obj, s, t)); });
    return refsMasked;
}

std::vector<std::shared_ptr<const AbstractSource>> Solution::getSources(size_t t, const Baseline& b) const noexcept {
    std::vector<std::shared_ptr<const AbstractSource>> mask;
    std::copy_if(src_.begin(), src_.end(), std::back_inserter(mask), 
        [this, t, &b](const std::shared_ptr<const AbstractSource>& obj) { 
            return this->getSol(Key::BlnActive(this, obj, b, t));
        });
    return mask;
}

size_t Solution::calculateMinObs(std::shared_ptr<const VieVS::AbstractSource> const q, const Baseline& b, size_t t) {
    size_t dur = std::numeric_limits<size_t>::max();
    for(auto& mode : modes_->getModes()) {
        size_t maxDuration = Solution::calculateMinObs(q, b, t, mode);
        dur = std::min(dur, maxDuration);
    }
    return dur;
}

size_t Solution::calculateMinObs(std::shared_ptr<const VieVS::AbstractSource> const q, const Baseline& b, size_t t, std::shared_ptr<const VieVS::Mode> const mode) {
    boost::optional<unsigned int> fixedScanDuration = q->getPARA().fixedScanDuration;
    if(auto fixedScanDuration = q->getPARA().fixedScanDuration) {
        return (*fixedScanDuration + blockLength_ - 1) / blockLength_ - 1;
    }

    // get stations
    unsigned long staid1 = b.getStaid1();
    const Station& sta1 = network_.refStation(staid1);
    unsigned long staid2 = b.getStaid2();
    const Station& sta2 = network_.refStation( staid2 );

    // calculate greenwhich meridian sedirial time
    double date1 = 2400000.5;
    double date2 = TimeSystem::mjdStart + static_cast<double>( t * blockLength_ ) / 86400.0;
    double gmst = iauGmst82( date1, date2 );

    // calculate system equivalent flux density for each station
    double el1 = pvs_.at(Key::StaActive(this, q, sta1, t)).getEl();
    double el2 = pvs_.at(Key::StaActive(this, q, sta2, t)).getEl();

    const auto& sta1Para = sta1.getPARA();
    const auto& sta2Para = sta2.getPARA();
    const auto& srcPara  = q->getPARA();
    const auto& blPara   = b.getParameters();

    double maxCorSynch1 = sta1Para.midob;
    double maxCorSynch2 = sta2Para.midob;
    // get maximum correlator synchronization time for
    double maxCorSynch = std::max( { maxCorSynch1, maxCorSynch2 } );

    double efficiency = mode->efficiency( sta1.getId(), sta2.getId() );

    unsigned int maxDuration = 0;
    for(auto &band : mode->getAllBands()) {
        double SEFD_src;
        if ( q->hasFluxInformation( band ) ) {
            // calculate observed flux density for each band
            SEFD_src = q->observedFlux( band, t * blockLength_, gmst, network_.getDxyz( staid1, staid2 ) );
        } else if ( ObservingMode::sourceBackup[band] == ObservingMode::Backup::internalModel ) {
            // calculate observed flux density based on model
            double wavelength = ObservingMode::wavelengths[band];
            SEFD_src = q->observedFlux_model( wavelength, t * blockLength_, gmst, network_.getDxyz( staid1, staid2 ) );
        } else {
            SEFD_src = 1e-3;
        }

        if ( SEFD_src == 0 ) {
            SEFD_src = 1e-3;
        }
        
        double SEFD_sta1 = sta1.getEquip().getSEFD( band, el1 );
        double SEFD_sta2 = sta2.getEquip().getSEFD( band, el2 );

        // get minimum required SNR for each station, baseline and source
        double minSNR_sta1 = sta1Para.minSNR.at( band );
        double minSNR_sta2 = sta2Para.minSNR.at( band );
        double minSNR_bl = blPara.minSNR.at( band );
        double minSNR_src = srcPara.minSNR.at( band );

        // maximum required minSNR
        double maxminSNR = std::max( { minSNR_src, minSNR_bl, minSNR_sta1, minSNR_sta2 } );

        // calc required baseline scan duration
        double anum = ( maxminSNR / ( SEFD_src * efficiency ) );
        double anu1 = SEFD_sta1 * SEFD_sta2;
        double anu2 = mode->recordingRate( staid1, staid2, band );
        if ( anu2 == 0 ) {
            return std::numeric_limits<size_t>::max();
        }
        double new_duration = anum * anum * anu1 / anu2 + maxCorSynch;
        new_duration = ceil( new_duration );
        auto new_duration_uint = static_cast<unsigned int>( new_duration );

        // check if required baseline scan duration is within min and max scan times of baselines
        unsigned int minScanBl = blPara.minScan;
        if ( new_duration_uint < minScanBl ) {
            new_duration_uint = minScanBl;
        }
        unsigned int maxScanBl = blPara.maxScan;
        if(new_duration_uint > maxDuration) {
            maxDuration = new_duration_uint;
        }
    }
    return (maxDuration + blockLength_ - 1) / blockLength_ - 1;
}
}
