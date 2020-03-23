#include "GeoConv.h"
#include <QtMath>
#include <QDebug>


double BlighterUtil::to_radian(double degree) {return degree * (M_PI/180.0);}
double BlighterUtil::to_degree(double radian) {return radian * (180.0/M_PI);}

double BlighterUtil::convertLatitudeDeg(UINT8 LatDegrees, UINT8 LatIntegerMinutes, UINT16 LatFractionalMinutes)
{
  //return (double)latLon_t->LatDegrees + (((double)latLon_t->LatIntegerMinutes) / 60) + (((double)latLon_t->LatFractionalMinutes) / 3600);
  double deg = static_cast<double>(LatDegrees);
  double min = stod(std::to_string(LatIntegerMinutes) + "." + std::to_string(LatFractionalMinutes));

  return deg + (min / 60.0);
}
double BlighterUtil::convertLatitudeDeg(LAT_LONG_T *latLon_t)
{
  //return (double)latLon_t->LatDegrees + (((double)latLon_t->LatIntegerMinutes) / 60) + (((double)latLon_t->LatFractionalMinutes) / 3600);
  double deg = static_cast<double>(latLon_t->LatDegrees);
  double min = stod(std::to_string(latLon_t->LatIntegerMinutes) + "." + std::to_string(latLon_t->LatFractionalMinutes));

  return deg + (min / 60.0);
}
double BlighterUtil::convertLongitudeDeg(UINT8 LongDegrees, UINT8 LongIntegerMinutes, UINT16 LongFractionalMinutes)
{
  //return (double)latLon_t->LongDegrees + (((double)latLon_t->LongIntegerMinutes) / 60) + (((double)latLon_t->LongFractionalMinutes) / 3600);
  double deg = static_cast<double>(LongDegrees);
  double min = stod(std::to_string(LongIntegerMinutes) + "." + std::to_string(LongFractionalMinutes));

  return deg + (min / 60.0);
}
double BlighterUtil::convertLongitudeDeg(LAT_LONG_T *latLon_t)
{
  //return (double)latLon_t->LongDegrees + (((double)latLon_t->LongIntegerMinutes) / 60) + (((double)latLon_t->LongFractionalMinutes) / 3600);
  double deg = static_cast<double>(latLon_t->LongDegrees);
  double min = stod(std::to_string(latLon_t->LongIntegerMinutes) + "." + std::to_string(latLon_t->LongFractionalMinutes));

  return deg + (min / 60.0);
}

bool BlighterUtil::calculateRelativeLatLon(const double currentLat, const double currentLon,
                                           const double range, const double azimuthDeg,
                                           double *destinationLat, double *destinationLon)
{
    bool ret = true;

    const double rad = to_radian(azimuthDeg);

    *destinationLat = currentLat + meters_to_latitude(range * qCos(rad));
    *destinationLon = currentLon + meters_to_longitude(range * qSin(rad));

    return ret;
}

bool BlighterUtil::calculateRelativeRangeAzimuth(const double currentLat, const double currentLon,
                                                 const double destinationLat, const double destinationLon,
                                                 double *relativeRange, double *relativeAzimuthDeg)
{
    bool ret = true;

    *relativeRange = calculateDistance(currentLat, currentLon, destinationLat, destinationLon);
    *relativeAzimuthDeg = calculateBearingP1toP2(currentLat, currentLon, destinationLat, destinationLon);

    return ret;
}

double BlighterUtil::calculateDistance(const double lat1, const double lon1, const double lat2, const double lon2)
{
    if ((lat1 == lat2)&&(lon1 == lon2))
    {
        return 0;
    }

    double e10 = lat1 * M_PI / 180.0;
    double e11 = lon1 * M_PI / 180.0;
    double e12 = lat2 * M_PI / 180.0;
    double e13 = lon2 * M_PI / 180.0;

    /* 타원체 GRS80 */
    double c16 = 6356752.314140910;
    double c15 = 6378137.000000000;
    double c17 = 0.0033528107;


    double f15 = c17 + c17 * c17;
    double f16 = f15 / 2.0;
    double f17 = c17 * c17 / 2.0;
    double f18 = c17 * c17 / 8.0;
    double f19 = c17 * c17 / 16.0;

    double c18 = e13 - e11;
    double c20 = (1.0 - c17) * qTan(e10);
    double c21 = qAtan(c20);
    double c22 = qSin(c21);
    double c23 = qCos(c21);
    double c24 = (1.0 - c17) * qTan(e12);
    double c25 = qAtan(c24);
    double c26 = qSin(c25);
    double c27 = qCos(c25);

    double c29 = c18;
    double c31 = (c27 * qSin(c29) * c27 * qSin(c29)) + (c23 * c26 - c22 * c27 * qCos(c29)) * (c23 * c26 - c22 * c27 * qCos(c29));
    double c33 = (c22 * c26) + (c23 * c27 * qCos(c29));
    double c35 = qSqrt(c31) / c33;
    double c36 = qAtan(c35);
    double c38 = 0.0;
    if (c31==0.0)
    {
        c38 = 0.0;
    }else{
        c38 = c23 * c27 * qSin(c29) / qSqrt(c31);
    }

    double c40 = 0.0;
    if ((qCos(qAsin(c38)) * qCos(qAsin(c38))) == 0.0)
    {
        c40 = 0.0;
    }else{
        c40 = c33 - 2.0 * c22 * c26 / (qCos(qAsin(c38)) * qCos(qAsin(c38)));
    }

    double c41 = qCos(qAsin(c38)) * qCos(qAsin(c38)) * (c15 * c15 - c16 * c16) / (c16 * c16);
    double c43 = 1 + c41 / 16384.0 * (4096.0 + c41 * (-768.0 + c41 * (320.0 - 175.0 * c41)));
    double c45 = c41 / 1024.0 * (256.0 + c41 * (-128.0 + c41 * (74.0 - 47.0 * c41)));
    double c47 = c45 * qSqrt(c31) * (c40 + c45 / 4.0 * (c33 * (-1.0 + 2.0 * c40 * c40) - c45 / 6.0 * c40 * (-3.0 + 4.0 * c31) * (-3.0 + 4.0 * c40 * c40)));
    double c50 = c17 / 16.0 * qCos(qAsin(c38)) * qCos(qAsin(c38)) * (4.0 + c17 * (4.0 - 3.0 * qCos(qAsin(c38)) * qCos(qAsin(c38))));
    double c52 = c18 + (1.0 - c50) * c17 * c38 * (qAcos(c33) + c50 * qSin(qAcos(c33)) * (c40 + c50 * c33 * (-1.0+ 2.0 * c40 * c40)));

    double c54 = c16 * c43 * (qAtan(c35) - c47);

    // return distance in meter
    return c54;
}

double BlighterUtil::calculateBearingP1toP2(const double lat1, const double lon1, const double lat2, const double lon2)
{
    // 현재 위치 : 위도나 경도는 지구 중심을 기반으로 하는 각도이기 때문에 라디안 각도로 변환한다.
    double Cur_Lat_radian = lat1 * (M_PI / 180.0);
    double Cur_Lon_radian = lon1 * (M_PI / 180.0);


    // 목표 위치 : 위도나 경도는 지구 중심을 기반으로 하는 각도이기 때문에 라디안 각도로 변환한다.
    double Dest_Lat_radian = lat2 * (M_PI / 180.0);
    double Dest_Lon_radian = lon2 * (M_PI / 180.0);

    // radian distance
    double radian_distance = 0.0;
    radian_distance = qAcos(qSin(Cur_Lat_radian) * qSin(Dest_Lat_radian) + qCos(Cur_Lat_radian) * qCos(Dest_Lat_radian) * qCos(Cur_Lon_radian - Dest_Lon_radian));

    // 목적지 이동 방향을 구한다.(현재 좌표에서 다음 좌표로 이동하기 위해서는 방향을 설정해야 한다. 라디안값이다.
    double radian_bearing = qAcos((qSin(Dest_Lat_radian) - qSin(Cur_Lat_radian) * qCos(radian_distance)) / (qCos(Cur_Lat_radian) * qSin(radian_distance)));		// acos의 인수로 주어지는 x는 360분법의 각도가 아닌 radian(호도)값이다.

    double true_bearing = 0.0;
    if (qSin(Dest_Lon_radian - Cur_Lon_radian) < 0.0)
    {
        true_bearing = radian_bearing * (180.0 / M_PI);
        true_bearing = 360.0 - true_bearing;
    }
    else
    {
        true_bearing = radian_bearing * (180.0 / M_PI);
    }

    return true_bearing;
}
