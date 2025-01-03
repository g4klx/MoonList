/*============================================================================*/
/*!\file
 *
 * \brief main.c - Simple demo program for Sun position using rectangular
 *                   coordinates
 *
 * \author  David Hoadley <vcrumble@westnet.com.au>
 *
 * \details
 *      See [the main page](@ref index) (or the end of this source file) for a
 *      description of this file.
 *
 * \copyright
 * \parblock
 * Copyright (c) 2020, David Hoadley <vcrumble@westnet.com.au>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * \endparblock
 *
 *==============================================================================
 */

 /* ANSI includes etc. */
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>

#include <map>

/* Local and project includes */
#include "general.h"
#include "sky.h"
#include "skyio.h"
#include "sun.h"
#include "moon.h"
#include "skyfast.h"
#include "vectors3d.h"

#define SITE_PRESSURE_hPa      1013.0   /*!< Atmospheric pressure (hPa = mbar)*/
#define SITE_TEMPERATURE_degC  3.0      /*!< Temperature (degrees Celsius)    */

char callsign[25];
double latitude = 0.0;
double longitude = 0.0;
double height = 0.0;
int increment = 10;
double timeZone = 0.0;
char outputFile[255];
int excludeHourStart = 0;
int excludeHourEnd = 0;

std::map<double, double> map;

bool readIniFile();
bool processLocator(const char* locator);
bool readElevationFile(const char* locator);
bool interpolate(double azimuth, double elevation);

void sunTopocentricFast(double             j2kUtc_d,
    const Sky_DeltaTs* dTs,
    const Sky_SiteProp* site,
    Sky_SiteHorizon* topo);


const char* days[] = {"Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"};

const char* months[] = {"January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"};

int main(int argc, char** argv)
{
    if (!readIniFile())
        return 1;

    FILE* fp = ::fopen(outputFile, "wt");
    if (fp == NULL) {
        ::fprintf(stderr, "Cannot open the Output File\n");
        return 1;
    }

    Sky_DeltaTs     deltaTs;        /* Differences between timescales */
    Sky_SiteProp    site;           /* Properties of the site */
    Sky_SiteHorizon topo;           /* Moon position, as seen from the site */
    double          j2kUtc_d;       /* Time, expressed as a count of days */

    /* Set up the various delta T values. */
    sky_initTime(37, 0.0, &deltaTs);

    /* Set up the fixed site data */
    sky_setSiteLocation(latitude, longitude, height, &site);
    sky_setSiteTempPress(SITE_TEMPERATURE_degC, SITE_PRESSURE_hPa, &site);
    sky_setSiteTimeZone(timeZone, &site);

    time_t now = ::time(NULL);

    if (argc > 1) {
        struct tm* tm = ::gmtime(&now);

        int month = atoi(argv[1]);
        tm->tm_mon = month - 1;

        if (argc > 2) {
            int year = atoi(argv[2]);
            tm->tm_year = year - 1900;
        }

        tm->tm_mday = 1;
        tm->tm_hour = 0;
        tm->tm_min = 0;

        now = ::mktime(tm);
    }

    struct tm* tm = ::gmtime(&now);

    int month = tm->tm_mon;
    int year = tm->tm_year;

    tm->tm_min /= increment;
    tm->tm_min *= increment;

    tm->tm_sec = 0;

    now = ::mktime(tm);

    bool isUp = false;

    ::fprintf(fp, "<!DOCTYPE html>\n<html>\n<head>\n<style>\n");
    ::fprintf(fp, "table, th, td {\nborder: 1px solid black;\nborder-collapse: collapse;\n}\n");
    ::fprintf(fp, "tr:nth-child(even) {\nbackground-color: #D6EEEE;\n}\n");
    ::fprintf(fp, "td {\ntext-align: center;\n}\n");
    ::fprintf(fp, "</style>\n<title>\nMoon Predictions for %s %d for %s\n</title>\n</head>\n<body>\n", months[month], year + 1900, callsign);
    ::fprintf(fp, "<h1>\nMoon Predictions for %s %d for %s\n</h1>\n", months[month], year + 1900, callsign);

    for (;;) {
        /* Get the current time in days since 2000-01-01 Greenwich noon */
        j2kUtc_d = sky_unixTimeToJ2kd(now);

        /* Get moon position */
        moon_nrelTopocentric(j2kUtc_d, &deltaTs, &site, &topo);

        if (topo.elevation_rad > degToRad(0.0)) {
            int y, m, d, h, mins;
            double s;
            sky_j2kdToCalTime(j2kUtc_d /* + site.timeZone_d */, &y, &m, &d, &h, &mins, &s);

            if ((h < excludeHourStart) || (h > excludeHourEnd)) {

                if (topo.azimuth_rad < 0.0)
                    topo.azimuth_rad += 2.0 * PI;

                if (interpolate(radToDeg(topo.azimuth_rad), radToDeg(topo.elevation_rad))) {
                    if (!isUp) {
                        ::fprintf(fp, "<h4>Starting on %s %02d/%02d</h4>\n", days[tm->tm_wday], d, m);
                        ::fprintf(fp, "<table style=\"width:35%%\"><tr>\n<th width=\"5%%\">Date</th>\n<th width=\"6%%\">UTC</th>\n<th width=\"9%%\">Moon Az.</th>\n<th width=\"9%%\">Moon El.</th>\n<th width=\"9%%\">Sun Az.</th>\n<th width=\"9%%\">Sun El.</th>\n</tr>\n");
                    }

                    isUp = true;

                    ::fprintf(fp, "<tr>\n<td>%02d/%02d</td>\n<td>%02d%02d</td>\n<td>%.0f</td>\n<td>%.0f</td>\n", d, m, h, mins, radToDeg(topo.azimuth_rad), radToDeg(topo.elevation_rad));

                    /* Get sun position */
                    sun_nrelTopocentric(j2kUtc_d, &deltaTs, &site, &topo);

                    if (topo.elevation_rad >= 0.0) {
                        /* Print out time and place */
                        if (topo.azimuth_rad < 0.0)
                            topo.azimuth_rad += 2.0 * PI;
                        ::fprintf(fp, "<td>%.0f</td>\n<td>%.0f</td>\n", radToDeg(topo.azimuth_rad), radToDeg(topo.elevation_rad));
                    } else {
                        ::fprintf(fp, "<td>-</td>\n<td>-</td>\n");
                    }

                    ::fprintf(fp, "</tr>\n");
                } else {
                    if (isUp) {
                        ::fprintf(fp, "</table>\n");
                        isUp = false;
                    }
                }
            } else {
                if (isUp) {
                    ::fprintf(fp, "</table>\n");
                    isUp = false;
                }
            }
        } else {
            if (isUp) {
                ::fprintf(fp, "</table>\n");
                isUp = false;
            }
        }

        now += increment * 60;

        tm = ::gmtime(&now);

        if ((tm->tm_mon != month) || (tm->tm_year != year))
            break;
    }

    ::fprintf(fp, "</body>\n</html>\n");

    ::fclose(fp);

    return 0;
}

void sunTopocentricFast(double             j2kUtc_d,
    const Sky_DeltaTs* dTs,
    const Sky_SiteProp* site,
    Sky_SiteHorizon* topo)
    /*! Uses previously calculated values of the Sun's apparent coordinates and
        interpolation to quickly calculate the Sun's topocentric position.
     \param[in]  j2kUtc_d  Days since Greenwich noon, 1-Jan-2000, UTC timescale
     \param[in]  deltas    Delta T values, as set by the sky_initTime() (or
                           sky_initTimeSimple() or sky_initTimeDetailed()) routines
     \param[in]  site      Properties of the observing site, particularly its
                           geometric location on the surface of the Earth, as set by
                           the sky_setSiteLocation() function (or sky_setSiteLoc2())
     \param[out] topo      Topocentric position, in both rectangular (unit vector)
                           form, and as Azimuth and Elevation (altitude) angles

     \par When to call this function
        This function is designed to be called at high frequency (more than once
        per second, even up to, say, 20 Hz).
     \par
        Use this function when you have made a previous call to skyfast_init().
        Depending on how long a period you want to track the sun using this
        routine, you may also need to run the function skyfast_backgroundUpdate() as
        a low frequency, low priority routine, to update the interpolation ranges.
        If you have not done this, call sun_nrelTopocentric() instead.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky_Times           atime;
    Sky_TrueEquatorial  approx;
    V3D_Vector          terInterV;// unit vector in Terrestrial Intermed Ref Sys

    REQUIRE_NOT_NULL(dTs);
    REQUIRE_NOT_NULL(site);
    REQUIRE_NOT_NULL(topo);

    sky_updateTimes(j2kUtc_d, dTs, &atime);

    /* Get the approximate apparent coordinates by interpolation */
    skyfast_getApprox(atime.j2kTT_cy, &approx);

    /* Convert apparent coordinates to topocentric coordinates at the site */
    sky0_appToTirs(&approx.appCirsV,
        atime.j2kUT1_d,
        approx.eqEq_rad,
        &terInterV);
    sky_siteTirsToTopo(&terInterV, approx.distance_au, site, topo);
}

/*!
 *  \mainpage
 *
 * \author  David Hoadley <vcrumble@westnet.com.au>
 *
 * \version $(DOXYGEN_APP_VER)
 *
 * \details
 *      This program calls two demo routines to calculate the Sun's position.
 *      The first (called demo1()) calculates it using the NREL SPA algorithm
 *      (see reference), but implemented in rectangular coordinates.
 *      The second (called demo2()) runs a loop to repeatedly calculate the
 *      position, but using interpolation between two previously calculated
 *      positions.
 *
 *      The interpolation will introduce a very small error in position. For the
 *      value of 720 minutes supplied to routine skyfast_init(), the maximum
 *      error that is introduced is less than 0.075 arcseconds.
 *      Given that the NREL-SPA algorithm itself is only
 *      accurate to approximately 1 arcsecond, this is pretty good.
 *
 *      This demo shows some of the use of the routines in sky.h and sun.h,
 *      and for interpolation, skyfast.h.
 *
 * \par Reference:
 *          Reda, Ibrahim and Andreas, Afshin.
 *          _Solar Position Algorithm for Solar Radiation Applications._
 *          National Renewable Energy Laboratory, publication no.
 *          NREL/TP-560-34302, June 2003, revised January 2008
 *
 * \copyright
 * \parblock
 * Copyright (c) 2020, David Hoadley <vcrumble@westnet.com.au>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * \endparblock
 */
 /*! \example demo1_sun.c Single calculation of the Sun's position
  *  \example demo2_sun_tracking.c
  *                      Repeated calculation (tracking) of the Sun using
  *                      interpolation
  *  \example demo3_moon.c
  *                      Single calculation of the Moon's position
  *  \example demo4_multiple_sites.c
  *                      Repeated calculation of the Sun's position for three
  *                      European sites simultaneously. This demo also shows the
  *                      use of the sky_initTimeDetailed() function.
  */

bool readIniFile()
{
    FILE* fp = ::fopen("MoonList.ini", "rt");
    if (fp == NULL) {
        ::fprintf(stderr, "Cannot open MoonList.ini\n");
        return false;
    }

    char buffer[255];
    while (::fgets(buffer, 255, fp) != NULL) {
        char* keyword = ::strtok(buffer, "=");
        char* value = ::strtok(NULL, "\r\n");

        if (keyword == NULL || value == NULL) {
            ::fprintf(stderr, "Malformed ini file\n");
            return false;
        }

        if (::strcmp(keyword, "Callsign") == 0) {
            ::strcpy(callsign, value);
        } else if (::strcmp(keyword, "Locator") == 0) {
            bool ret = ::processLocator(value);
            if (!ret)
                return false;
        } else if (::strcmp(keyword, "Height") == 0) {
            height = ::atof(value);
        } else if (::strcmp(keyword, "TimeZone") == 0) {
            timeZone = ::atof(value);
        } else if (::strcmp(keyword, "Increment") == 0) {
            increment = ::atoi(value);
        } else if (::strcmp(keyword, "ElevationFile") == 0) {
            bool ret = readElevationFile(value);
            if (!ret)
                return false;
        } else if (::strcmp(keyword, "OutputFile") == 0) {
            ::strcpy(outputFile, value);
        }  else if (::strcmp(keyword, "ExcludeHourStart") == 0) {
            excludeHourStart = ::atoi(value);
        } else if (::strcmp(keyword, "ExcludeHourEnd") == 0) {
            excludeHourEnd = ::atoi(value);
        }
    }

    return true;
}

bool processLocator(const char* locator)
{
    assert(locator != NULL);

    if (::strlen(locator) != 6)
        return false;

    if (locator[0] < 'A' || locator[0] > 'R' || locator[1] < 'A' || locator[1] > 'R' ||
        locator[2] < '0' || locator[2] > '9' || locator[3] < '0' || locator[3] > '9' ||
        locator[4] < 'A' || locator[4] > 'X' || locator[5] < 'A' || locator[5] > 'X')
        return false;

    longitude  = double(locator[0] - 'A') * 20.0;
    latitude   = double(locator[1] - 'A') * 10.0;
    longitude += double(locator[2] - '0') * 2.0;
    latitude  += double(locator[3] - '0') * 1.0;
    longitude += double(locator[4] - 'A') * 0.0833;
    latitude  += double(locator[5] - 'A') * 0.0417;

    latitude  -= 90.0;
    longitude -= 180.0;

    return true;
}

bool readElevationFile(const char* fileName)
{
    assert(fileName != NULL);

    FILE* fp = ::fopen(fileName, "rt");
    if (fp == NULL) {
        ::fprintf(stderr, "Cannot open %s\n", fileName);
        return false;
    }

    char buffer[255];
    while (::fgets(buffer, 255, fp) != NULL) {
        if (buffer[0] == '#')
            continue;

        char* pAzimuth   = ::strtok(buffer, " \t");
        char* pElevation = ::strtok(NULL, "\r\n");

        if (pAzimuth == NULL || pElevation == NULL) {
            ::fprintf(stderr, "Malformed elevation file\n");
            return false;
        }

        double azimuth   = ::atof(pAzimuth);
        double elevation = ::atof(pElevation);

        map.insert(std::pair<double, double>(azimuth, elevation));
    }

    ::fclose(fp);

    return true;
}

bool interpolate(double azimuth, double elevation)
{
    if (map.empty())
        return true;

    std::map<double, double>::const_iterator it = map.lower_bound(azimuth);

    if (it == map.end()) {
        double el = map.rbegin()->second;
        return el < elevation;
    }

    if (it == map.begin()) {
        double el = it->second;
        return el < elevation;
    }

    double x2 = it->first;
    double y2 = it->second;

    --it;
    double x1 = it->first;
    double y1 = it->second;

    double p = (azimuth - x1) / (x2 - x1);

    double el = ((1.0 - p) * y1) + (p * y2);

    return el < elevation;
}
