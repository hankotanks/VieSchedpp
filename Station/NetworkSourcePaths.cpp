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

    NetworkSourcePaths::NetworkSourcePaths(Network& network, SourceList& sources) {
        this->source_paths_ = std::map<unsigned long, std::vector<StationSourcePath>>();
        for(auto& sta : network.refStations()) {
            std::vector<StationSourcePath> sta_paths = std::vector<StationSourcePath>();
            for(auto& src : sources.refSources()) sta_paths.push_back(StationSourcePath(sta, src));
            this->source_paths_.insert(std::pair<unsigned long, std::vector<StationSourcePath>>(sta.getId(), sta_paths));
        }
    }
};