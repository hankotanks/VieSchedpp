/**
 * @file Position.h
 * @brief class Position
 *
 *
 * @author Matthias Schartner
 * @date 23.06.2017
 */

#ifndef POSITION_H
#define POSITION_H
#include <cmath>
#include <iostream>
#include <boost/format.hpp>
#include "Constants.h"
#include "VieVS_Object.h"
#include "sofa.h"

namespace VieVS{
    /**
     * @class Position
     * @brief representation of VLBI station position
     *
     * @author Matthias Schartner
     * @date 23.06.2017
     */
    class Position: public VieVS_Object {
    public:

        /**
         * @brief constructor
         *
         * @param x_m x coordinate in meters
         * @param y_m y coordinate in meters
         * @param z_m z coordinate in meters
         */
        Position(double x_m, double y_m, double z_m);


        /**
         * @brief getter for x coordinate
         *
         * @return x coordinate in meters
         */
        double getX() const noexcept {
            return x_;
        }

        /**
         * @brief getter for y coordinate
         *
         * @return y coordinate in meters
         */
        double getY() const noexcept {
            return y_;
        }

        /**
         * @brief getter for z coordinate
         *
         * @return z coordinate in meters
         */
        double getZ() const noexcept {
            return z_;
        }

        /**
         * @brief getter for latitude
         *
         * @return latitude in radians
         */
        double getLat() const noexcept {
            return lat_;
        }

        /**
         * @brief getter for longitude
         *
         * @return longitude in radians
         */
        double getLon() const noexcept {
            return lon_;
        }

        /**
         * @brief calculates distance between two stations
         *
         * @param other second station
         * @return distance between stations
         */
        double getDistance(const Position &other) const noexcept;


        void geodetic2Local(double g2l[3][3]){
            g2l = g2l_;
        }

        const std::vector<std::vector<double>> getGeodetic2Local() const{
            return g2l_2;
        }

    private:
        static unsigned long nextId;

        double x_; ///< x coordinate in meters
        double y_; ///< y coordinate in meters
        double z_; ///< z coordinate in meters
        double lat_; ///< latitude in radians
        double lon_; ///< longitude in radians
        double h_; ///< height in meters

        double g2l_[3][3]; ///< geocentric to local transformation matrix

        std::vector<std::vector<double>> g2l_2; ///< geocentric to local transformation matrix

    };
}
#endif /* POSITION_H */

