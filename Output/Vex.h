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
 * @file Vex.h
 * @brief class Vex
 *
 * @author Matthias Schartner
 * @date 07.12.2017
 */

#ifndef VEX_H
#define VEX_H


#include "../ObservingMode/ObservingMode.h"
#include "../Scan/Scan.h"
#include "../Station/Station.h"


namespace VieVS {

/**
 * @class Vex
 * @brief this is the VLBI vex-output class
 *
 * @author Matthias Schartner
 * @date 07.12.2017
 */

class Vex : public VieVS_Object {
   public:
    /**
     * @brief constructor
     * @author Matthias Schartner
     *
     * @param file file name
     */
    explicit Vex( const std::string &file );


    /**
     * @brief writ vex file
     * @author Matthias Schartner
     *
     * @param network station network
     * @param sources list of all sources
     * @param scans list of all scans
     * @param obsModes observing mode
     * @param xml paramters.xml file
     */
    void writeVex( const Network &network, const std::vector<Source> &sources, const std::vector<Scan> &scans,
                   const std::shared_ptr<const ObservingMode> &obsModes, const boost::property_tree::ptree &xml );


   private:
    static unsigned long nextId;  ///< next id for this object type

    std::ofstream of;                   ///< output stream object *filename*.vex
    std::string eol = ";\n";            ///< end of line string
    std::map<int, int> channelNr2Bbc_;  ///< channel number to bbc number

    /**
     * @brief write vex $GLOBAL block
     * @author Matthias Schartner
     *
     * @param expName experiment name
     */
    void global_block( const std::string &expName );


    /**
     * @brief write vex $EXPER block
     * @author Matthias Schartner
     *
     * @param expName experiment name
     * @param expDescription experiment description
     * @param piName pi name
     * @param piEmail pi email
     * @param contactName contact name
     * @param contactEmail contact email
     * @param schedulerName scheduler name
     * @param schedulerEmail scheduler email
     * @param notes additional notes
     * @param targetCorrelator target correlator
     * @param gui_version gui scheduler version
     */
    void exper_block( const std::string &expName, const std::string &expDescription, const std::string &piName,
                      const std::string &piEmail, const std::string &contactName, const std::string &contactEmail,
                      const std::string &schedulerName, const std::string &schedulerEmail, const std::string &notes,
                      const std::string &targetCorrelator, const std::string &gui_version );


    /**
     * @brief write vex $STATION block
     * @author Matthias Schartner
     *
     * @param stations list of all stations
     */
    void station_block( const std::vector<Station> &stations );


    /**
     * @brief write vex $STATION block
     * @author Matthias Schartner
     *
     * @param stations list of all stations
     */
    void sites_block( const std::vector<Station> &stations );


    /**
     * @brief write vex $ANTENNA block
     * @author Matthias Schartner
     *
     * @param stations list of all stations
     */
    void antenna_block( const std::vector<Station> &stations );


    /**
     * @brief write vex $DAS block
     * @author Matthias Schartner
     *
     * @param stations list of all stations
     */
    void das_block( const std::vector<Station> &stations );


    /**
     * @brief write vex $SOURCE block
     * @author Matthias Schartner
     *
     * @param sources list of all sources
     */
    void source_block( const std::vector<Source> &sources );


    /**
     * @brief write vex $MODE block
     * @author Matthias Schartner
     *
     * @param obsModes observing mode
     */
    void mode_block( const std::shared_ptr<const ObservingMode> &obsModes );


    /**
     * @brief write vex $FREQ block
     * @author Matthias Schartner
     *
     * @param obsModes observing mode
     */
    void freq_block( const std::shared_ptr<const ObservingMode> &obsModes );


    /**
     * @brief write vex $BBC block
     * @author Matthias Schartner
     *
     * @param obsModes observing mode
     */
    void bbc_block( const std::shared_ptr<const ObservingMode> &obsModes );


    /**
     * @brief write vex $IF block
     * @author Matthias Schartner
     *
     * @param obsModes observing mode
     */
    void if_block( const std::shared_ptr<const ObservingMode> &obsModes );


    /**
     * @brief write vex $TRACKS block
     * @author Matthias Schartner
     *
     * @param obsModes observing mode
     */
    void tracks_block( const std::shared_ptr<const ObservingMode> &obsModes );


    /**
     * @brief write vex head block (deprecated)
     * @author Matthias Schartner
     */
    void head_pos_block();


    /**
     * @brief write $PASS_ORDER block
     * @author Matthias Schartner
     */
    void pass_order_block();


    /**
     * @brief write $ROLL block
     * @author Matthias Schartner
     */
    void roll_block();


    /**
     * @brief write $PHASE_CAL_DETECT block
     * @author Matthias Schartner
     */
    void phase_cal_detect_block();


    /**
     * @brief write $SCHED block
     * @author Matthias Schartner
     *
     * @param scans list of all scans
     * @param network station network
     * @param sources list of all sources
     * @param obsModes observing mode
     */
    void sched_block( const std::vector<Scan> &scans, const Network &network, const std::vector<Source> &sources,
                      const std::shared_ptr<const ObservingMode> &obsModes );
};
}  // namespace VieVS

#endif  // VEX_H
