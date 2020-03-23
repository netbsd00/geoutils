#ifndef GEOCONV_H
#define GEOCONV_H

#include <QObject>

class GeoConv
{
public:
  static double convertLatitudeDeg(UINT8 LatDegrees, UINT8 LatIntegerMinutes, UINT16 LatFractionalMinutes);
  static double convertLatitudeDeg(LAT_LONG_T *latLon_t);
  static double convertLongitudeDeg(UINT8 LongDegrees, UINT8 LongIntegerMinutes, UINT16 LongFractionalMinutes);
  static double convertLongitudeDeg(LAT_LONG_T *latLon_t);
  static bool calculateRelativeLatLon(const double currentLat, const double currentLon,
                                      const double range, const double azimuthDeg,
                                      double *destinationLat, double *destinationLon);
  static bool calculateRelativeRangeAzimuth(const double currentLat, const double currentLon,
                                            const double destinationLat, const double destinationLon,
                                            double *relativeRange, double *relativeAzimuthDeg);
  static double calculateDistance(const double lat1, const double lon1, const double lat2, const double lon2);
  static double calculateBearingP1toP2(const double lat1, const double lon1, const double lat2, const double lon2);

  /* 근사값 */
  static double meters_to_latitude(const double meter) {return static_cast<double>(meter * 0.000009);}
  static double meters_to_longitude(const double meter) {return static_cast<double>(meter * 0.000011);}
  static double latitude_to_meters(const double latDelta) {return latDelta / 0.000009;}
  static double longitude_to_meters(const double lonDelta) {return lonDelta / 0.000011;}
  static double to_radian(double degree);
  static double to_degree(double radian);
};

#endif // GEOCONV_H
