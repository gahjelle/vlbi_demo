#!/usr/bin/env python3
"""@package coords
Library module for handling of different coordinate systems

Example:
import coords
deg = coords.dms_to_deg(59, 54, 52.32)

Description:

This module will provide functions for handling coordinate systems and
transformations between them.

The function names use the following conventions:

deg     - degrees
dms     - degrees, minutes, seconds
hms     - hours, minutes, seconds
rad     - radians

References:
[1] SOFA Tools for Earth Attitude.
    http://www.iausofa.org/sofa_pn.pdf

[2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
    IERS Technical Note No. 36, BKG (2010).
    http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

Author:
    Geir Arne Hjelle <geir.arne.hjelle@kartverket.no>

$Revision: 10525 $
$Date: 2016-05-19 15:19:38 +0200 (Thu, 19 May 2016) $
$LastChangedBy: dahmic $
"""

# Standard library imports
from datetime import datetime

# External library imports
from astropy.time import Time
import numpy as np

# WHERE imports
import sofa


def dms_to_deg(degrees, minutes, seconds, sign=None):
    """Convert degrees, minutes and seconds to degrees.

    The input parametes can be given as numbers or strings. Degrees
    and minutes will be cast to int before calculation, while seconds
    will be interpreted as floats.

    The sign can be specified in several ways. For instance will True,
    '+', '', and 1 all represent a positive sign, while False, '-',
    and -1 will represent a negative sign. See _parse_sign for more
    information.

    If sign is not specified, the sign of degrees will be used. In
    this case, be careful that for instance the sign of +0 or -0 is
    correctly passed on. That is, degrees must be specified as either
    a string or a float, not an int. In general, it is safer to pass
    the sign explicitly.

    Args:
        degrees:   Degrees as int or string representing int.
        minutes:   Minutes as int or string representing int.
        seconds:   Seconds as float or string representing float.
        sign:      Sign of degrees.

    Returns:
        Float representing degrees.

    Examples:
        >>> dms_to_deg('59', '54', '52.32')
        59.91453333333333
        >>> dms_to_deg(12, 34, 56.789, sign='-')
        -12.582441388888888
        >>> dms_to_deg(-0.0, 19, 59.974870)
        -0.3333263527777778
    """
    if sign is None:
        sign = degrees
    return _parse_sign(sign) * (abs(int(degrees)) + int(minutes) / 60 +
                                float(seconds) / 3600)


def dms_to_rad(degrees, minutes, seconds, sign=None):
    """Convert degrees, minutes and seconds to radians.

    The input parametes can be given as numbers or strings. Degrees
    and minutes will be cast to int before calculation, while seconds
    will be interpreted as floats.

    The sign can be specified in several ways. For instance will True,
    '+', '', and 1 all represent a positive sign, while False, '-',
    and -1 will represent a negative sign. See _parse_sign for more
    information.

    If sign is not specified, the sign of degrees will be used. In
    this case, be careful that for instance the sign of +0 or -0 is
    correctly passed on. That is, degrees must be specified as either
    a string or a float, not an int. In general, it is safer to pass
    the sign explicitly.

    Args:
        degrees:   Degrees as int or string representing int.
        minutes:   Minutes as int or string representing int.
        seconds:   Seconds as float or string representing float.

    Returns:
        Float representing radians.

    Examples:
        >>> dms_to_rad(59, 54, '52.32')
        1.04570587646256
        >>> dms_to_rad('-12', '34', 56.789)
        -0.21960503017531938
        >>> dms_to_rad(0, 19, 59.974870, '-')
        -0.005817642339636369
    """
    return np.radians(dms_to_deg(degrees, minutes, seconds, sign))


def hms_to_deg(hours, minutes, seconds):
    """Convert hours, minutes and seconds to degrees.

    The input parametes can be given as numbers or strings. Hours and
    minutes will be cast to int before calculation, while seconds will
    be interpreted as floats.

    Args:
        hours:     Hours as int or string representing int.
        minutes:   Minutes as int or string representing int.
        seconds:   Seconds as float or string representing float.

    Returns:
        Float representing degrees.

    Examples:
        >>> hms_to_deg('17', '7', '17.753427')
        256.82397261250003
        >>> hms_to_deg(12, 0, 0)
        180.0
        >>> hms_to_deg(-12, 34, 56.789)
        Traceback (most recent call last):
        ValueError: hours must be non-negative
    """
    if float(hours) < 0:
        raise ValueError('hours must be non-negative')
    return 15 * (int(hours) + int(minutes) / 60 + float(seconds) / 3600)


def hms_to_rad(hours, minutes, seconds):
    """Convert hours, minutes and seconds to radians.

    The input parametes can be given as numbers or strings. Hours and
    minutes will be cast to int before calculation, while seconds will
    be interpreted as floats.

    Args:
        hours:     Hours as int or string representing int.
        minutes:   Minutes as int or string representing int.
        seconds:   Seconds as float or string representing float.

    Returns:
        Float representing radians.

    Examples:
        >>> hms_to_rad(17, 7, 17.753427)
        4.482423920139868
        >>> hms_to_rad('12', '0', '0.00')
        3.141592653589793
    """
    return np.radians(hms_to_deg(hours, minutes, seconds))


def arcsec_to_rad(arcsec):
    """Convert arc seconds to radians.

    Args:
        arcsec:   Arcseconds as float.

    Returns:
        Float representing radians.

    Examples:
        >>> arcsec_to_rad(1)
        4.84813681109536e-06
        >>> arcsec_to_rad(12345.6)
        0.05985315781505887
    """
    return arcsec * np.pi / (180 * 60 * 60)


def gcrs_to_itrs(obsdt, eop=None):
    """Calculate the GCRS to ITRS transformation matrix.

    The calculation is based on the "IAU 2006/2000A, CIO based, using
    X, Y series" example in section 5.6 of SOFA Tools for Earth
    Attitude [1]. More information is found in Chapter 5, and in
    particular section 5.5 of the IERS Conventions [2].

    Warning: This simplified calculation does not handle time scales correctly.

    Args:
        obsdt_utc:   Datetime object representing the observation date.

    Returns:
        (3, 3)-ndarray with the GCRS to ITRS transformation matrix.

    Examples:
        >>> from datetime import datetime
        >>> gcrs_to_itrs(datetime(2007, 4, 5, 12, 0, 0), \
                         eop=dict(x=0.0349282, y=0.4833163, \
                                  ut1_utc=-0.072073685))
        array([[  9.73104318e-01,   2.30363826e-01,  -7.03163482e-04],
               [ -2.30363800e-01,   9.73104571e-01,   1.18545368e-04],
               [  7.11560163e-04,   4.66264024e-05,   9.99999746e-01]])
    """
    # Modified Julian dates
    obsjd, obsfrac = _dt_to_jd(obsdt)

    # CIP and CIO
    x, y = sofa.iau_xy06(obsjd, obsfrac)
    s = sofa.iau_s06(obsjd, obsfrac, x, y)

    # Add Celestial Intermediate Pole corrections
    x += arcsec_to_rad(eop['dx'])
    y += arcsec_to_rad(eop['dy'])

    # Celestial pole motion
    transformation = sofa.iau_c2ixys(x, y, s)

    # Earth rotation angle
    era = sofa.iau_era00(obsjd, obsfrac)
    sofa.iau_rz(era, transformation)

    # Polar motion    # apriori/eop/eopc04_iau
    xp = arcsec_to_rad(eop['x'])
    yp = arcsec_to_rad(eop['y'])
    polar_motion = sofa.iau_pom00(xp, yp, sofa.iau_sp00(obsjd, obsfrac))

    # Return final matrices
    gcrs2itrs = np.dot(polar_motion, transformation)

    return gcrs2itrs


def itrs_to_gcrs(obsdt_utc, eop=None):
    """Calculate the ITRS to GCRS transformation matrix.

    The calculation is based on the "IAU 2006/2000A, CIO based, using
    X, Y series" example in section 5.6 of SOFA Tools for Earth
    Attitude (http://www.iausofa.org/sofa_pn.pdf).

    Args:
        obsdt_utc: Datetime object representing the observation date in UTC.

    Returns:
        (3, 3)-ndarray with the ITRS to GCRS transformation matrix.

    Examples:
        >>> from datetime import datetime
        >>> itrs_to_gcrs(datetime(2007, 4, 5, 12, 0, 0), \
                         eop=dict(x=0.0349282, y=0.4833163, \
                                  ut1_utc=-0.072073685))
        array([[  9.73104318e-01,  -2.30363800e-01,   7.11560163e-04],
               [  2.30363826e-01,   9.73104571e-01,   4.66264024e-05],
               [ -7.03163482e-04,   1.18545368e-04,   9.99999746e-01]])
    """
    return gcrs_to_itrs(obsdt_utc, eop=eop).T


def _parse_sign(sign):
    """Interpret sign as either +1 or -1.

    The sign can be specified in several ways. For instance will True,
    '+', '', and 1 all represent a positive sign, while False, '-',
    and -1 will represent a negative sign.

    Args:
        sign:   A bool, string or number to interpret as a sign

    Returns:
        An int, either 1 to represent a positive sign or -1 for negative.

    Examples:
        >>> _parse_sign(True)
        1
        >>> _parse_sign('')
        1
        >>> _parse_sign('+')
        1
        >>> _parse_sign(42)
        1
        >>> _parse_sign(0)
        1
        >>> _parse_sign(False)
        -1
        >>> _parse_sign('-')
        -1
        >>> _parse_sign(-1)
        -1
        >>> _parse_sign(-0.0)
        -1
    """
    # Deal with bools on their own
    if isinstance(sign, bool):
        return 1 if sign else -1

    number = float(str(sign) + '1')
    return (number > 0) - (number < 0)


# Julian day
def _dt_to_jd(given_date, scale='utc', scale_out=None):
    """Calculate the Julian day for a given date.

    Use the astropy module to calculate the Julian day for a given
    date. The given date is assumed to be given in UTC, but this can
    be changed by specifying scale. If the output should be in a
    different timescale than input, then scale_out can be specified.
    Timescales can be any of 'tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1',
    and 'utc'.

    Args:
        given_date:   Datetime object representing the given date.
        scale:        Timescale of given date.
        scale_out:    Timescale of returned Julian Day.

    Returns:
        Tuple of floats representing the Julian day for the given date.

    Example:
        >>> _dt_to_jd(datetime(2014, 5, 6))
        (2456783.5, 0.0)
        >>> dt_to_jd(datetime(2014, 5, 6), scale='utc', scale_out='tt')
        (2456783.5, 0.0007775924168527126)
    """
    if scale_out is None:
        jd = Time(given_date, scale=scale).jd
    else:
        jd = Time(Time(given_date, scale=scale), scale=scale_out).jd

    # Offset by 0.5 since Julian day starts at noon
    jd_int, jd_frac = divmod(jd - 0.5, 1)
    return jd_int + 0.5, jd_frac


# Run doctests if module is run as script
if __name__ == '__main__':
    import doctest
    import sys
    sys.exit(doctest.testmod()[0])  # Non-zero exit code if tests fail
