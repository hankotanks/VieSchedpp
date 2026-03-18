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
 * @file SchedulerILP.h
 * @brief class SchedulerILP
 *
 * @author Hank Lewis
 * @date 24.05.2025
 */

#ifndef SCHEDULER_ILP_H
#define SCHEDULER_ILP_H
#include <boost/date_time.hpp>
#include <boost/optional.hpp>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <limits>

#include "Initializer.h"
#include "Model.h"
#include "Scheduler.h"

namespace VieVS {
/**
 * @class SchedulerILP
 * @brief this is an overload of the original Scheduler class that implements
 * globally optimal scheduling using integer linear programming (ILP)
 *
 * @author Hank Lewis
 * @date 24.05.2025
 */
class SchedulerILP : public Scheduler {
    friend class Output;
private:
    void initialize(const std::set<unsigned long>& sourceMask) noexcept;
public:
    /**
    * @brief constructor
    * @author Hank Lewis
    *
    * @param init initializer
    * @param path path to VieSchedpp.xml file
    * @param fname file name
    */
    SchedulerILP( Initializer &init, std::string path, std::string fname, std::ofstream& statisticsOf ) : 
        Scheduler(init, path, fname), statisticsOf_(statisticsOf) { /* STUB */ }

    ~SchedulerILP() {
        if(model_) delete model_;
    }

    /**
     * @brief main function that starts the scheduling
     * @author Hank Lewis
     */
    void start() noexcept override;

    /**
     * @brief updates the selected next scans to the schedule
     * @author Hank Lewis
     *
     * @param scan best possible next scans
     * @param of outstream file object
     */
    void update( Scan &scan, std::ofstream &of ) noexcept override;

    /**
     * @brief statistics output
     * @author Hank Lewis
     *
     * @param of output stream
     */
     void statistics( std::ofstream &of ) override;

private:
     Model* model_ = nullptr;
     std::ofstream& statisticsOf_;
};
}
#endif // SCHEDULER_ILP_H