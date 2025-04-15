#include "NetworkSourcePaths.h"
#include <memory>
#include <utility>
#include <vector>
#include <iostream>

#include "../Scan/PointingVector.h"
#include "../Station/Network.h"
#include "../Source/SourceList.h"

namespace VieVS {
    StationSourcePath::StationSourcePath() {
        this->srcid_ = std::numeric_limits<unsigned long>::max();
        this->staid_ = std::numeric_limits<unsigned long>::max();
        this->pvs_ = std::vector<VieVS::PointingVector>();
    }

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
        this->src_paths_ = std::map<unsigned long, std::map<unsigned long, StationSourcePath>>();
        this->sta_names_ = std::map<std::string, unsigned long>();
        for(auto& src : sources.refSources()) {
            this->src_names_.insert(std::pair<std::string, unsigned long>(src->getName(), src->getId()));
        }
        for(auto& sta : network.refStations()) {
            std::map<unsigned long, StationSourcePath> sta_paths = std::map<unsigned long, StationSourcePath>();
            for(auto& src : sources.refSources()) {
                sta_paths.insert(std::pair<unsigned long, StationSourcePath>(src->getId(), StationSourcePath(sta, src)));
            }
            this->src_paths_.insert(std::pair<unsigned long, std::map<unsigned long, StationSourcePath>>(sta.getId(), sta_paths));
            this->sta_names_.insert(std::pair<std::string, unsigned long>(sta.getName(), sta.getId()));
        }
    }

    std::map<unsigned long, StationSourcePath>& NetworkSourcePaths::getAllPaths(std::string stationName) {
        return getAllPaths(this->sta_names_[stationName]);
    }

    std::map<unsigned long, StationSourcePath>& NetworkSourcePaths::getAllPaths(unsigned long stationId) {
        return this->src_paths_[stationId];
    }

    StationSourcePath& NetworkSourcePaths::getPath(unsigned long stationId, unsigned long sourceId) {
        return getAllPaths(stationId)[sourceId];
    }

    StationSourcePath& NetworkSourcePaths::getPath(std::string stationName, unsigned long sourceId) {
        return getAllPaths(stationName)[sourceId];
    }
};