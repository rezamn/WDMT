from pyproj import Geod


def delaz(lon_1, lat_1, lon_2, lat_2, ellipse='WGS84'):
    """
    Calculates azimuth and back azimuth in degrees and distance in km  from input coordinates

    Parameters
    ----------
    :param lon_1:
        longitude of first point
        :type lon_1: float
    :param lat_1:
        latitude of first point
        :type lat_1: float
    :param lon_2:
        longitude of second point
        :type lon_2: float
    :param lat_2:
        latitude of second point
        :type lat_2: float
    :param ellipse:
        ellipsoid of the earth for more information see
        :type ellipse: str

    Returns
    --------
    :return:
        returns (azimuth, back azimuth, dist)
        :rtype: (float, float, float)
    """
    geod = Geod(ellps=ellipse)
    az, baz, dist = geod.inv(lon_1, lat_1, lon_2, lat_2)
    dist_km = dist / 1000
    return dist_km, az, baz


def find_coordnate(lon_1, lat_1, az12, distance, ellipse='WGS84'):
    """
    Calculates the coordinate of a point base on reference point, azimuth and distance

    Parameters
    ----------
    :param lon_1:
        longitude of first point
        :type lon_1: float
    :param lat_1:
        latitude of first point
        :type lat_1: float
    :param az12:
        azimuth from reference point
        :type az12: float
    :param distance:
        distance in km from reference point
        :type distance: float
    :param ellipse:
        ellipsoid of the earth for more information see
        :type ellipse: str

    Returns
    --------
    :return:
        returns (longitude, latitude)
        :rtype: (float, float)
    """
    geod = Geod(ellps=ellipse)
    lon_2, lat_2, baz = geod.fwd(lon_1, lat_1, az12, distance*1000)

    return lon_2, lat_2, baz


if __name__ == "__main__":
    lon_a = 0
    lat_a = 0
    lon_b = 0
    lat_b = 0
    print(delaz(0, 0, 1, 0))
    print(delaz(0, 0, 0, 1))
