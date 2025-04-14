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
        StationSourcePath(Station& sta, std::shared_ptr<const AbstractSource> src);
    private:
        static const unsigned int step = 600;
        unsigned long srcid_;
        unsigned long staid_;
        std::vector<PointingVector> pvs_;
    };

    class NetworkSourcePaths {
    public:
        NetworkSourcePaths(VieVS::Network& network, VieVS::SourceList& sources);
    private:
        std::map<unsigned long, std::vector<StationSourcePath>> source_paths_;
    };
};

#endif /* NETWORK_SOURCE_PATHS_H */