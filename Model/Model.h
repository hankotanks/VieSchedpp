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

// gurobi
#ifdef WITH_GUROBI
#include <gurobi_c++.h>
#endif // WITH_GUROBI

// VieSchedpp
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
    /**
     * @brief override of ModelBase::prepare, called once the model is initialized
     * @author Hank Lewis
    */
    virtual void prepare(const Window& window) override;

private:
#ifdef WITH_GUROBI
    /**
     * @brief ensure no station makes two concurrent observations
     * @author Hank Lewis
    */
    void constrExclusive(const Window& window);

    /**
     * @brief ensure StaActive is tied to BlnActive
     * @author Hank Lewis
    */
    void constrBaseline(const Window& window);

    /**
     * @brief ensure minNumberOfSites is respected
     * @author Hank Lewis
    */
    void constrPairwise(const Window& window);

    /**
     * @brief ensure max scan duration is respected
     * @author Hank Lewis
    */
    void constrDuration(const Window& window);

    /**
     * @brief ensure sufficient time for observations
     * @author Hank Lewis
    */
    void constrSNR(const Window& window);

    /**
     * @brief ensure enough time to slew between observations
     * @author Hank Lewis
    */
    void constrSlew(const Window& window);

    /**
     * @brief ensure StaCoverage tied to StaActive
     * @author Hank Lewis
    */
    void constrCoverage(const Window& window);

private:
    /**
     * @brief sky coverage objected
     * @author Hank Lewis
    */
    GRBLinExpr objSkyCov(const Window& window);

    /**
     * @brief baseline objective, weighted by baseline length
     * @author Hank Lewis
    */
    GRBLinExpr objBaselines(const Window& window);
#endif // WITH_GUROBI
};
}
#endif // MODEL_H