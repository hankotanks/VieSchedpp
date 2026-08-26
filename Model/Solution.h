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

/**
 * @file Model.h
 * @brief class Model
 *
 * @author Hank Lewis
 * @date 22.01.2026
 */

#ifndef SOLUTION_H
#define SOLUTION_H

// system
#include <functional>
#include <unordered_set>

// VieSchedpp
#include "../Source/SourceList.h"
#include "../Station/Network.h"
#include "../Scan/Scan.h"

namespace VieVS {

/**
 * @class Solution
 * @brief computes variable feasibility and stores accumulated solution
 * @author Hank Lewis
 * @date 24.05.2025
 */
class Solution {
    friend class ModelBase;

public:
    /**
     * @brief constructor
     * @author Hank Lewis
     * 
     * @param network network
     * @param sourceList sourceList
     * @param sourceMask sourceMask
     * @param modes modes
     * @param blockLength blockLength
     * @param windowLength windowLength
     */
    Solution(VieVS::Network& network, VieVS::SourceList& sourceList, 
        const std::set<unsigned long>& sourceMask, 
        const std::shared_ptr<const ObservingMode>& modes,
        unsigned int blockLength, unsigned int windowLength);

public:
    /**
    * @class Key
    * @brief key for accessing ILP variables and solutions
    * @author Hank Lewis
    * @date 03.03.2026
    */
    struct Key {
        // represents the type of variable
        enum Type { sta_active, bln_active, sta_coverage };
        Type type;

        // the identifier
        union {
            struct { size_t q, s, t; } sta_active;
            struct { size_t q, b, t; } bln_active;
            struct { size_t s, c; } sta_coverage;
        } key;

        // a name for pretty printing
        std::string name;

        /**
         * @brief StaActive builder
         * @author Hank Lewis
         *
         * @param sol the solution map
         * @param q source
         * @param s station
         * @param t time segment
         */
        static Key StaActive(const Solution* sol, 
            std::shared_ptr<const VieVS::AbstractSource> const q, 
            const Station& s, size_t t);

        /**
         * @brief BlnActive builder
         * @author Hank Lewis
         *
         * @param sol the solution map
         * @param q source
         * @param b baseline
         * @param t time segment
         */
        static Key BlnActive(const Solution* sol, 
            std::shared_ptr<const VieVS::AbstractSource> const q, 
            const Baseline& b, size_t t);

        /**
         * @brief StaCoverage builder
         * @author Hank Lewis
         *
         * @param sol the solution map
         * @param s station
         * @param c cell of sky plot
         */
        static Key StaCoverage(const Solution* sol, const Station& s, size_t c);
        
        /**
         * @brief operator<
         * @author Hank Lewis
         */
        bool operator<(const Key& other) const;

        /**
         * @brief operator==
         * @author Hank Lewis
         */
        bool operator==(const Key& other) const noexcept;

        /**
         * @struct hash implementation
         * @author Hank Lewis
         */
        struct Hash {
            size_t operator()(const Key& k) const noexcept {
                size_t h = std::hash<int>{}(static_cast<int>(k.type));
                auto combine = [&](size_t v) {
                    h ^= std::hash<size_t>{}(v)
                    + 0x9e3779b97f4a7c15ULL
                    + (h << 6)
                    + (h >> 2);
                };
                switch(k.type) {
                case Key::sta_active:
                    combine(k.key.sta_active.q);
                    combine(k.key.sta_active.s);
                    combine(k.key.sta_active.t);
                    break;
                case Key::bln_active:
                    combine(k.key.bln_active.q);
                    combine(k.key.bln_active.b);
                    combine(k.key.bln_active.t);
                    break;
                case Key::sta_coverage:
                    combine(k.key.sta_coverage.s);
                    combine(k.key.sta_coverage.c);
                    break;
                }
                return h;
            }
        };
    };

public:
    /**
     * @brief returns a mutable value from the solution map (if one exists)
     * @author Hank Lewis
     * 
     * @param key key
     */
    boost::optional<bool&> getSol(const Key& key) noexcept;

    /**
     * @brief returns a value from the solution map (if one exists)
     * @author Hank Lewis
     * 
     * @param key key
     */
    boost::optional<const bool&> getSol(const Key& key) const noexcept;

    /**
     * @brief returns the number of adjacent time segments for a BlnActive variable
     * @author Hank Lewis
     * 
     * @param key key
     */
    size_t getMinObs(const Key& key) const noexcept;

    /**
     * @brief returns the number of empty segments required between two StaActive observations
     * @author Hank Lewis
     * 
     * @param from from
     * @param to to
     */
    size_t getSlew(const Key& from, const Key& to) const noexcept;

    /**
     * @brief returns pointing vector at the start of a time segment
     * @author Hank Lewis
     * 
     * @param key key
     */
    boost::optional<PointingVector> getPointingVector(const Key& key) const noexcept;

public:
    /**
     * @brief returns all segments in an interval
     * @author Hank Lewis
     * 
     * @param t0 t0
     * @param tf tf
     */
    std::vector<size_t> getBlocks(size_t t0, size_t tf) const noexcept;

    /**
     * @brief returns all segments in an interval where s can observe q
     * @author Hank Lewis
     * 
     * @param t0 t0
     * @param tf tf
     * @param q q
     * @param s s
     */
    std::vector<size_t> getBlocks(size_t t0, size_t tf, const std::shared_ptr<const AbstractSource>& q, const Station& s) const noexcept;
    
    /**
     * @brief returns all segments in an interval where b can observe q
     * @author Hank Lewis
     * 
     * @param t0 t0
     * @param tf tf
     * @param q q
     * @param b b
     */
    std::vector<size_t> getBlocks(size_t t0, size_t tf, const std::shared_ptr<const AbstractSource>& q, const Baseline& b) const noexcept;
    
    /**
     * @brief returns the number of blocks needed to cover an interval of n seconds
     * @author Hank Lewis
     * 
     * @param seconds seconds
     */
    size_t getBlocks(unsigned int seconds) const noexcept;

    /**
     * @brief returns all stations
     * @author Hank Lewis
     */
    std::vector<std::reference_wrapper<const Station>> getStations() const noexcept;

    /**
     * @brief retrieves the two stations comprising a baseline
     * @author Hank Lewis
     * 
     * @param b b
     */
    std::pair<std::reference_wrapper<const Station>, std::reference_wrapper<const Station>> getStations(const Baseline& b) const noexcept;
    
    /**
     * @brief returns all stations that can observe q at t
     * @author Hank Lewis
     * 
     * @param t t
     * @param q q
     */
    std::vector<std::reference_wrapper<const Station>> getStations(size_t t, const std::shared_ptr<const AbstractSource>& q) const noexcept;
    
    /**
     * @brief returns all baselines
     * @author Hank Lewis
     */
    std::vector<std::reference_wrapper<const Baseline>> getBaselines() const noexcept;
    
    /**
     * @brief returns all baselines at t where q can be observed
     * @author Hank Lewis
     * 
     * @param t t
     * @param q q
     */
    std::vector<std::reference_wrapper<const Baseline>> getBaselines(size_t t, const std::shared_ptr<const AbstractSource>& q) const noexcept;
    
    /**
     * @brief returns all sources
     * @author Hank Lewis
     */
    std::vector<std::shared_ptr<const AbstractSource>> getSources() const noexcept;
    
    /**
     * @brief returns all sources that can be observed by s at t
     * @author Hank Lewis
     * 
     * @param t t
     * @param s s
     */
    std::vector<std::shared_ptr<const AbstractSource>> getSources(size_t t, const Station& s) const noexcept;

    /**
     * @brief returns all sources that can observed by b at t
     * @author Hank Lewis
     * 
     * @param t t
     * @param b b
     */
    std::vector<std::shared_ptr<const AbstractSource>> getSources(size_t t, const Baseline& b) const noexcept;

private:
    /**
     * @brief returns the number of segments needed to make an observation
     * @author Hank Lewis
     * 
     * @param q q
     * @param b b
     * @param t t
     */
    size_t calculateMinObs(std::shared_ptr<const VieVS::AbstractSource> const q, const Baseline& b, size_t t);

    /**
     * @brief returns the number of segments needed to make an observation under a specific mode
     * @author Hank Lewis
     * 
     * @param q q
     * @param b b
     * @param t t
     * @param mode mode
     */
    size_t calculateMinObs(std::shared_ptr<const VieVS::AbstractSource> const q, const Baseline& b, size_t t, 
        std::shared_ptr<const VieVS::Mode> const mode);

private:
    // stations and baselines
    VieVS::Network& network_;
    // complete list of sources
    VieVS::SourceList& sourceList_;
    // source mask
    std::set<unsigned long> sourceMask_;
    // observing modes
    std::shared_ptr<const ObservingMode> modes_;
    // number of blocks in the schedule
    size_t blockCount_;
    // the size of each time segment (in seconds)
    unsigned int blockLength_;
    // number of optimization windows
    size_t windowCount_;
    // the length of the sliding optimization window (in seconds)
    unsigned int windowLength_;
    // number of segments in each window
    size_t windowBlockCount_;
    // bidirectional mapping of indices and ids
    std::map<unsigned long, size_t> sta2idx_;
    std::map<unsigned long, size_t> bln2idx_;
    std::map<unsigned long, size_t> src2idx_;
    std::vector<unsigned long> idx2sta_;
    std::vector<unsigned long> idx2bln_;
    std::vector<unsigned long> idx2src_;
    // vectors containing complete sets of station, source, and baseline references
    std::vector<std::reference_wrapper<const Station>> sta_;
    std::vector<std::reference_wrapper<const Baseline>> bln_;
    std::vector<std::shared_ptr<const AbstractSource>> src_;
    // the solution map (StaActive and BlnActive)
    std::map<Key, bool> sol_;
    // precomputed pointing vectors (StaActive)
    std::unordered_map<Key, PointingVector, Key::Hash> pvs_;
    // station visibility (StaActive)
    std::unordered_set<Key, Key::Hash> vis_;
    // minObs by baseline (BlnActive)
    std::unordered_map<Key, size_t, Key::Hash> snr_;
};
}

#endif // SOLUTION_H