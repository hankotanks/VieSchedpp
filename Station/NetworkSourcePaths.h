#ifndef NETWORK_SOURCE_PATHS_H
#define NETWORK_SOURCE_PATHS_H

#include <memory>
#include <utility>
#include <vector>

#include "../Scan/PointingVector.h"
#include "../Station/Network.h"
#include "../Source/SourceList.h"

namespace VieVS {
    class StationSourcePath {
    public:
        StationSourcePath();
        StationSourcePath(Station& sta, std::shared_ptr<const AbstractSource> src);
        unsigned long getSourceId() const;
        unsigned long getStationId();
        const std::vector<PointingVector>& getVectors() const;
    private:
        static const unsigned int step = 1800;
        unsigned long srcid_;
        unsigned long staid_;
        std::vector<PointingVector> pvs_;
    };

    class NetworkSourcePaths {
    public:
        NetworkSourcePaths(VieVS::Network& network, VieVS::SourceList& sources);
        std::map<unsigned long, StationSourcePath>& getAllPaths(std::string stationName);
        std::map<unsigned long, StationSourcePath>& getAllPaths(unsigned long stationId);
        StationSourcePath& getPath(unsigned long stationId, unsigned long sourceId);
        StationSourcePath& getPath(std::string stationName, unsigned long sourceId);
    private:
        std::map<unsigned long, std::map<unsigned long, StationSourcePath>> src_paths_;
        std::map<std::string, unsigned long> sta_names_;
        std::map<std::string, unsigned long> src_names_;
    };
};

#endif /* NETWORK_SOURCE_PATHS_H */