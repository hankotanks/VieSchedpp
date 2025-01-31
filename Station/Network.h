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
 * @file Network.h
 * @brief class Network
 *
 * @author Matthias Schartner
 * @date 28.06.2018
 */

#ifndef NETWORK_H
#define NETWORK_H


#include <utility>
#include <vector>

#include "../Misc/VieVS_Object.h"
#include "../ObservingMode/ObservingMode.h"
#include "Baseline.h"
#include "SkyCoverage.h"
#include "Station.h"


namespace VieVS {

/**
 * @class Network
 * @brief representation of a VLBI network
 *
 * @author Matthias Schartner
 * @date 28.06.2018
 */

class Network : public VieVS_Object {
   public:
    /**
     * @brief constructor
     * @author Matthias Schartner
     */
    Network();


    /**
     * @brief add a new station to the network
     * @author Matthias Schartner
     *
     * @param station new station which should be added
     */
    void addStation( Station station );


    /**
     * @brief get station per id
     * @author Matthias Schartner
     *
     * @param id station id
     * @return station with this id
     */
    const Station &getStation( unsigned long id ) const noexcept;


    /**
     * @brief get station per name
     * @author Matthias Schartner
     *
     * @param name station name
     * @return station with this name
     */
    const Station &getStation( const std::string &name ) const noexcept;


    /**
     * @brief get all stations in this network
     * @author Matthias Schartner
     *
     * @return all station in this network
     */
    const std::vector<Station> &getStations() const noexcept;


    /**
     * @brief ger baseline per id
     * @author Matthias Schartner
     *
     * @param id baseline id
     * @return baseline with this id
     */
    const Baseline &getBaseline( unsigned long id ) const noexcept;


    /**
     * @brief get baselined between two stations
     * @author Matthias Schartner
     *
     * @param staid1 first station id
     * @param staid2 second station id
     * @return baseline between those stations
     */
    const Baseline &getBaseline( unsigned long staid1, unsigned long staid2 ) const noexcept;


    /**
     * @brief get baseline between two stations
     * @author Matthias Schartner
     *
     * @param staids station ids
     * @return baseline between those stations
     */
    const Baseline &getBaseline( const std::pair<unsigned long, unsigned long> &staids ) const noexcept;


    /**
     * @brief get baseline per name
     * @author Matthias Schartner
     *
     * @param name name of the baseline
     * @return baseline with this name
     */
    const Baseline &getBaseline( const std::string &name ) const noexcept;


    /**
     * @brief get all baselines
     * @author Matthias Schartner
     *
     * @return all baselines
     */
    const std::vector<Baseline> &getBaselines() const noexcept;


    /**
     * @brief get sky coverage per id
     * @author Matthias Schartner
     *
     * @param id sky coverage id
     * @return sky coverage with this id
     */
    const SkyCoverage &getSkyCoverage( unsigned long id ) const noexcept;


    /**
     * @brief get all sky coverages
     * @author Matthias Schartner
     *
     * @return all sky coverage id
     */
    const std::vector<SkyCoverage> &getSkyCoverages() const noexcept;


    /**
     * @brief reference to station per id
     * @author Matthias Schartner
     *
     * @param id station id
     * @return station with this id
     */
    Station &refStation( unsigned long id );


    /**
     * @brief reference to station per name
     * @author Matthias Schartner
     *
     * @param name station name
     * @return station with this name
     */
    Station &refStation( const std::string &name );


    /**
     * @brief reference to all stations
     * @author Matthias Schartner
     *
     * @return all stations
     */
    std::vector<Station> &refStations();


    /**
     * @brief reference to baseline per id
     * @author Matthias Schartner
     *
     * @param id baseline id
     * @return baseline with this id
     */
    Baseline &refBaseline( unsigned long id );


    /**
     * @brief reference to baseline between two stations
     * @author Matthias Schartner
     *
     * @param staid1 first station id
     * @param staid2 second station id
     * @return baseline between those stations
     */
    Baseline &refBaseline( unsigned long staid1, unsigned long staid2 );


    /**
     * @brief reference to baseline between two stations
     * @author Matthias Schartner
     *
     * @param staids station ids
     * @return baseline between those stations
     */
    Baseline &refBaseline( const std::pair<unsigned long, unsigned long> &staids );


    /**
     * @brief reference to baseline per name
     * @author Matthias Schartner
     *
     * @param name baseline name
     * @return baseline with this name
     */
    Baseline &refBaseline( const std::string &name );


    /**
     * @brief reference to all baselines
     * @author Matthias Schartner
     *
     * @return all baselines
     */
    std::vector<Baseline> &refBaselines();


    /**
     * @brief reference to sky coverage per id
     * @author Matthias Schartner
     *
     * @param id sky coverage id
     * @return sky coverage with this id
     */
    SkyCoverage &refSkyCoverage( unsigned long id );


    /**
     * @brief reference to all sky coverages
     * @author Matthias Schartner
     *
     * @return all sky coverages
     */
    std::vector<SkyCoverage> &refSkyCoverages();


    /**
     * @brief get baseline id between stations
     * @author Matthias Schartner
     *
     * @param staids station ids
     * @return baseline id between these stations
     */
    unsigned long getBlid( const std::pair<unsigned long, unsigned long> &staids ) const noexcept;


    /**
     * @brief get baseline id between stations
     * @author Matthias Schartner
     *
     * @param staid1 first station id
     * @param staid2 second station id
     * @return baseline id between these stations
     */
    unsigned long getBlid( unsigned long staid1, unsigned long staid2 ) const noexcept;


    /**
     * @brief set maximum distance between corresponding telescopes
     * @author Matthias Schartner
     *
     * corresponding telescopes share the same sky coverage information
     *
     * @param maxDistBetweenCorrespondingTelescopes maximum distance in meters
     */
    void setMaxDistBetweenCorrespondingTelescopes( double maxDistBetweenCorrespondingTelescopes ) {
        maxDistBetweenCorrespondingTelescopes_ = maxDistBetweenCorrespondingTelescopes;
    }


    /**
     * @brief get number of stations
     * @author Matthias Schartner
     *
     * @return number of stations
     */
    unsigned long getNSta() const noexcept { return nsta_; }


    /**
     * @brief get number of baselines
     * @author Matthias Schartner
     *
     * @return number of baselines
     */
    unsigned long getNBls() const noexcept { return nbls_; }


    /**
     * @brief get number of theoretically possible independent closure phases and amplitudes
     * @author Matthias Schartner
     *
     * @return number of theoretically possible independent closure phases and amplitudes
     */
    unsigned long getNClosures_max() const noexcept { return nclosures_max_; }


    /**
     * @brief add new observation to station
     * @author Matthias Schartner
     *
     * @param nObs number of observations
     * @param pointingVector observed pointing vector
     * @param influence flag if this observations should count for upcoming scan selections
     */
    void update( unsigned long nObs, const PointingVector &pointingVector, bool influence = true );


    /**
     * @brief add new observation to baseline
     * @author Matthias Schartner
     *
     * @param blid baseline id
     * @param influence flag if observation should count for upcoming scan selections
     */
    void update( unsigned long blid, bool influence = true );


    /**
     * @brief get baseline vector between two stations
     * @author Matthias Schartner
     *
     * @param staid1 first station id
     * @param staid2 second station id
     * @return baseline vector between these stations
     */
    const std::vector<double> &getDxyz( unsigned long staid1, unsigned long staid2 ) const;


    /**
     * @brief calc total score per sky coverages
     * @author Matthias Schartner
     *
     * @param pvs list of pointing vectors
     * @return total score
     */
    double calcScore_skyCoverage( const std::vector<PointingVector> &pvs ) const;


    /**
     * @brief calc total score per sky coverages
     * @author Matthias Schartner
     *
     * save result in input parameter for further use
     *
     * @param pvs list of pointing vectors
     * @param staids2skyCoverageScore storage for results
     * @return total score
     */
    double calcScore_skyCoverage( const std::vector<PointingVector> &pvs,
                                  std::unordered_map<unsigned long, double> &staids2skyCoverageScore ) const;


    /**
     * @brief calc total score per sky coverages
     * @author Matthias Schartner
     *
     * score is taken from input parameter
     *
     * @param pvs list of pointing vectors
     * @param staids2skyCoverageScore precalculated score
     * @return total score
     */
    double calcScore_skyCoverage_subnetting(
        const std::vector<PointingVector> &pvs,
        const std::unordered_map<unsigned long, double> &staids2skyCoverageScore ) const;


    /**
     * @brief get lookup table for sky coverage id based on station id
     * @author Matthias Schartner
     *
     * @return lookup table
     */
    const std::map<unsigned long, unsigned long> &getStaid2skyCoverageId() const { return staids2skyCoverageId_; };


    /**
     * @brief creates station summary table
     * @author Matthias Schartner
     *
     * @param of outfile stream
     */
    void stationSummary( std::ofstream &of ) const;

    /**
     * @brief add sky coverage object
     * @author Matthias Schartner
     *
     * @param name sky-coverage name (for connection with stations only)
     * @param distance influence distance
     * @param time influence time
     * @param dist_type influence distance function
     * @param time_type influence time function
     */
    void addSkyCoverage( const std::string &name, double distance, double time, SkyCoverage::Interpolation dist_type,
                         SkyCoverage::Interpolation time_type );

    /**
     * @brief connect sky coverage id and station name
     * @author Matthias Schartner
     *
     * @param sta2id key = station name, value = sky-coverage id
     */
    void connectSkyCoverageWithStation( const std::map<std::string, std::string> &sta2id );

    /**
     * @brief connect sky coverage id and station name
     * @author Matthias Schartner
     *
     * @param twinDistance maximum distance distance between telescopes sharing same sky-coverage
     * @param dist influence distance
     * @param time influence time
     * @param distType distance function
     * @param timeType distance time
     */
    void addSkyCoverages( double twinDistance, double dist, double time, SkyCoverage::Interpolation distType,
                          SkyCoverage::Interpolation timeType );

    static int stationIdsToNSites( std::vector<unsigned long> staids ) {
        std::set<int> sites;
        for ( auto staid : staids ) {
            sites.insert( station2site[staid] );
        }
        return sites.size();
    };

    static int nSites() {
        std::set<int> sites;
        for ( const auto &any : station2site ) {
            sites.insert( any.second );
        }
        return sites.size();
    };

    static const std::map<unsigned long, int> &stations2sites() { return station2site; }

    static void addSites( std::map<unsigned long, int> map ) { station2site = std::move( map ); }


   private:
    unsigned long nsta_;                     ///< number of stations
    unsigned long nbls_;                     ///< number of baselines
    unsigned long nclosures_max_;            ///< maximum number of possible closure delays
    std::vector<Station> stations_;          ///< all stations
    std::vector<Baseline> baselines_;        ///< all baselines
    std::vector<SkyCoverage> skyCoverages_;  ///< all sky coverages

    std::map<std::pair<unsigned long, unsigned long>, unsigned long> staids2blid_;  ///< lookup table for baseline id

    static unsigned long nextId;                       ///< next id for this object type
    static std::map<unsigned long, int> station2site;  ///< station to site converter

    std::map<std::pair<unsigned long, unsigned long>, std::vector<double>>
        staids2dxyz_;  ///< lookup table for baseline vectors

    double maxDistBetweenCorrespondingTelescopes_;  ///< maximum distance between corresponding telescopes in meteres
    std::map<unsigned long, unsigned long> staids2skyCoverageId_;  ///< lookup table for sky coverage ids
};
}  // namespace VieVS

#endif  // NETWORK_H
