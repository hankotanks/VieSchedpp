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

#ifndef MODEL_H
#define MODEL_H
#include <cassert>
#include <unordered_map>
#ifdef WITH_GUROBI

#include <gurobi_c++.h>

#include "Source/SourceList.h"
#include "Station/Network.h"

namespace VieVS {
/**
 * @class ModelKey
 * @brief key for accessing ILP variables
 * @author Hank Lewis
 * @date 03.03.2026
 */
struct ModelKey {
    enum ModelKeyType {
        sta_active,
        bln_active,
        coverage
    };
    ModelKeyType type;
    union {
        struct {
            unsigned long q;
            unsigned long s;
            size_t t;
        } sta_active;
        struct {
            unsigned long q;
            unsigned long s1;
            unsigned long s2;
            size_t t;
        } bln_active;
        struct {
            unsigned long s;
            size_t c;
        } coverage;
    } key;

    bool operator<(const ModelKey& other) const {
        if(type != other.type) return type < other.type;
        switch(type) {
            case ModelKey::ModelKeyType::sta_active:
                if(key.sta_active.q != other.key.sta_active.q) return key.sta_active.q < other.key.sta_active.q;
                if(key.sta_active.s != other.key.sta_active.s) return key.sta_active.s < other.key.sta_active.s;
                return key.sta_active.t < other.key.sta_active.t;
            case ModelKey::ModelKeyType::bln_active:
                if(key.bln_active.q != other.key.bln_active.q) return key.bln_active.q < other.key.bln_active.q;
                if(key.bln_active.s1 != other.key.bln_active.s1) return key.bln_active.s1 < other.key.bln_active.s1;
                if(key.bln_active.s2 != other.key.bln_active.s2) return key.bln_active.s2 < other.key.bln_active.s2;
                return key.bln_active.t < other.key.bln_active.t;
            case ModelKey::ModelKeyType::coverage:
                if(key.coverage.s != other.key.coverage.s) return key.coverage.s < other.key.coverage.s;
                return key.coverage.c < other.key.coverage.c;
            default: throw std::logic_error("unreachable");
        }
    }

    static ModelKey StaActive(std::shared_ptr<const VieVS::AbstractSource> const q, 
        const Station& s, size_t t) {
        ModelKey key;
        key.type = ModelKey::ModelKeyType::sta_active;
        key.key.sta_active.q = q->getId();
        key.key.sta_active.s = s.getId();
        key.key.sta_active.t = t;
        return key;
    }

    static ModelKey BlnActive(std::shared_ptr<const VieVS::AbstractSource> const q, 
        const Station& s1, const Station& s2, size_t t) {
        ModelKey key;
        key.type = ModelKey::ModelKeyType::bln_active;
        key.key.bln_active.q = q->getId();
        key.key.bln_active.s1 = s1.getId();
        key.key.bln_active.s2 = s2.getId();
        assert(key.key.bln_active.s1 != key.key.bln_active.s2);
        if(key.key.bln_active.s1 > key.key.bln_active.s2) 
            std::swap(key.key.bln_active.s1, key.key.bln_active.s2);
        key.key.bln_active.t = t;
        return key;
    }

    static ModelKey Coverage(const Station& s, size_t c) {
        ModelKey key;
        key.type = ModelKey::ModelKeyType::coverage;
        key.key.coverage.s = s.getId();
        key.key.coverage.c = c;
        return key;
    }
};

/**
 * @class Model
 * @brief wraps ILP implementation
 * @author Hank Lewis
 * @date 24.05.2025
 */
class Model {
public:
    /**
     * @brief constructor
     * @author Hank Lewis
     *
     * @param network network
     * @param sourceList sourceList
     * @param blockLength blockLength
     */
    Model(VieVS::Network& network, VieVS::SourceList& sourceList, size_t blockLength);

    /**
     * @brief destructor
     * @author Hank Lewis
    */
    ~Model() {
        delete model_;
        delete env_;
    }

private:
    bool checkStationVisibility(size_t t, 
        std::shared_ptr<const VieVS::AbstractSource> q, 
        const Station& s) const noexcept;
    size_t calculateCell(size_t t, 
        const std::shared_ptr<const AbstractSource> q,
        Station& s) const noexcept;
    size_t calculateSlewTime(Station& s, 
        const std::shared_ptr<const AbstractSource> q1, 
        const std::shared_ptr<const AbstractSource> q2,
        size_t t1, size_t t2) const noexcept;
private:
    boost::optional<GRBVar> getVar(const ModelKey& key) const noexcept;
private:
    // references from VieVS::Scheduler
    VieVS::Network& network_;
    VieVS::SourceList& sourceList_;
    // the size of each time segment (in seconds)
    unsigned int blockLength_;
    // number of segments in schedule
    unsigned int blockCount_;
private:
    std::map<ModelKey, int> var_;
private:
    // gurobi environment
    GRBEnv* env_ = nullptr;
    GRBModel* model_ = nullptr;
};
}
#endif // WITH_GUROBI
#endif // MODEL_H