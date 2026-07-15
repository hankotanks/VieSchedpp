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

#ifdef WITH_GUROBI
#include <gurobi_c++.h>
#endif // WITH_GUROBI

#include "ModelBase.h"

#include "../Source/SourceList.h"
#include "../Station/Network.h"

namespace VieVS {
/**
 * @class Model
 * @brief wraps ILP implementation
 * @author Hank Lewis
 * @date 24.05.2025
 */
class Model : public ModelBase {
    friend struct ModelCoverage;
    friend struct ModelCoverage13;
public:
    /**
     * @brief constructor
     * @author Hank Lewis
     *
     * @param network network
     * @param sourceList sourceList
     * @param sourceMask sourceMask
     * @param blockLength blockLength
     * @param windowLength windowLength
     */
    Model(VieVS::Network& network, VieVS::SourceList& sourceList, 
        const std::set<unsigned long>& sourceMask, 
        const std::shared_ptr<const ObservingMode>& modes,
        unsigned int blockLength, unsigned int windowLength) : ModelBase(network, sourceList, sourceMask, modes, blockLength, windowLength) { /* STUB */ }

    /**
     * @brief constructor with specified ModelCoverage implementation
     * @author Hank Lewis
     *
     * @param network network
     * @param sourceList sourceList
     * @param sourceMask sourceMask
     * @param blockLength blockLength
     * @param windowLength windowLength
     */
    template<typename T>
    Model(VieVS::Network& network, VieVS::SourceList& sourceList, 
        const std::set<unsigned long>& sourceMask, 
        const std::shared_ptr<const ObservingMode>& modes,
        unsigned int blockLength, unsigned int windowLength) : ModelBase(network, sourceList, sourceMask, modes, blockLength, windowLength) {
        static_assert(std::is_base_of<ModelCoverage, T>::value, "unreachable");
        coverage_ = std::make_unique<T>();
    }

    /**
     * @brief destructor
     * @author Hank Lewis
    */
    ~Model() = default;
protected:
    virtual void prepare(size_t tp, size_t t0, size_t tf, size_t tn) override;
private:
#ifdef WITH_GUROBI
    void constrExclusive(size_t t0, size_t tf);
    void constrBaseline(size_t t0, size_t tf);
    void constrPairwise(size_t t0, size_t tf);
    void constrDuration(size_t tp, size_t t0, size_t tf, size_t tn);
    void constrSNR(size_t tp, size_t t0, size_t tf, size_t tn);
    void constrSlew(size_t tp, size_t t0, size_t tf, size_t tn);
    void constrCoverage(size_t tp, size_t t0, size_t tf, size_t tn);
private:
    GRBLinExpr objSkyCov();
    GRBLinExpr objBaselines(size_t t0, size_t tf);
private:
    std::set<ModelKey> fixed_;
#endif // WITH_GUROBI
};
}
#endif // MODEL_H