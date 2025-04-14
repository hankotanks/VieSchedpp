#include "NetworkSourcePaths.h"
#include <memory>
#include <utility>
#include <vector>
#include <iostream>

#include "../Scan/PointingVector.h"
#include "../Station/Network.h"
#include "../Source/SourceList.h"

namespace VieVS {
    StationSourcePath::StationSourcePath(Station& sta, std::shared_ptr<const AbstractSource> src) {
        this->staid_ = sta.getId();
        this->srcid_ = src->getId();
        for(unsigned int t = 0; t < TimeSystem::duration + 1800; t += StationSourcePath::step) {
            PointingVector curr(sta.getId(), src->getId());
            curr.setTime(t);
            sta.calcAzEl_rigorous(src, curr);
            this->pvs_.push_back(curr);
        }
    }

    unsigned long StationSourcePath::getSourceId() const {
        return this->srcid_;
    }

    unsigned long StationSourcePath::getStationId() {
        return this->staid_;
    }

    const std::vector<PointingVector>& StationSourcePath::getVectors() const {
        return this->pvs_;
    }

    NetworkSourcePaths::NetworkSourcePaths(Network& network, SourceList& sources) {
        this->source_paths_ = std::map<unsigned long, std::vector<StationSourcePath>>();
        this->source_names_ = std::map<std::string, unsigned long>();
        for(auto& sta : network.refStations()) {
            std::vector<StationSourcePath> sta_paths = std::vector<StationSourcePath>();
            for(auto& src : sources.refSources()) sta_paths.push_back(StationSourcePath(sta, src));
            this->source_paths_.insert(std::pair<unsigned long, std::vector<StationSourcePath>>(sta.getId(), sta_paths));
            this->source_names_.insert(std::pair<std::string, unsigned long>(sta.getName(), sta.getId()));
        }
    }

    std::vector<StationSourcePath>& NetworkSourcePaths::getStationSourcePaths(std::string stationName) {
        unsigned long idx = this->source_names_[stationName];
        return this->source_paths_[idx];
    }
};