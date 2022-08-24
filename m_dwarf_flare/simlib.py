import re
import pandas as pd
import astropy_healpix
from pytz import country_timezones
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS


class SIMLIB_OBS:

    def __init__(self, data):

        # Finding the LIBID in the data
        self.LIBID = re.search('LIBID:(.+?)\n', data)
        if self.LIBID:
            self.LIBID = int(self.LIBID.group(1))
        
        # Finding the RA in the data
        self.RA = re.search('RA:(.+?) DECL', data)
        if self.RA:
            self.RA = float(self.RA.group(1))

        # Finding the DECL in the data
        self.DECL = re.search('DECL:(.+?) NOBS', data)
        if self.DECL:
            self.DECL = float(self.DECL.group(1))

        # Finding the NOBS in the data
        self.NOBS = re.search('NOBS:(.+?) MWEBV', data)
        if self.NOBS:
            self.NOBS = int(self.NOBS.group(1))
        
        # Finding the MWEBV in the data
        self.MWEBV = re.search('MWEBV:(.+?) PIXSIZE', data)
        if self.MWEBV:
            self.MWEBV = float(self.MWEBV.group(1))

        # Finding the PIXSIZE in the data
        self.PIXSIZE = re.search('PIXSIZE:(.+?)\n', data)
        if self.PIXSIZE:
            self.PIXSIZE = float(self.PIXSIZE.group(1))

        # Store the row information
        self.lines = data.split('\n')
        self.columns = self.lines[3]
        self.obs_rows = self.lines[4:-2]

    def isInCIPixel(self, nside, healpix_indices):
        
        self.nside = nside

        # Find the healpix pixel of this observation
        coordinates = SkyCoord(ra=self.RA, dec=self.DECL, frame=ICRS, unit='deg')
        map = astropy_healpix.HEALPix(self.nside, frame=ICRS, order="nested")
        self.HEALPIX_PIXEL = map.skycoord_to_healpix(coordinates, return_offsets=False)

        # Check if the healpix pixel is in the high ci region
        return self.HEALPIX_PIXEL in healpix_indices

    
    def countInTimeRange(self, timeStart, timeEnd):
        
        count = 0

        # Check to see if any observation is in the time range
        for i in range(len(self.obs_rows)):

            vals = self.obs_rows[i].split(' ')
            MJD = float(vals[1])


            # If any observation lies within the range, count
            if  MJD >= timeStart and MJD <= timeEnd:
                count += 1

        # Return the count
        return count


    def writeSubSampledObs(self, timeStart, timeEnd, new_ncount, libid, file):
        
        file.write("# --------------------------------------------\n")
        # Writing the LIBID
        file.write('LIBID: {}\n'.format(self.LIBID))

        # Writing the first row
        file.write("RA: {RA} DECL: {DECL} NOBS: {NOBS} MWEBV: {MWEBV} PIXSIZE: {PIXSIZE}\n".\
                    format(RA=self.RA, DECL=self.DECL, NOBS=new_ncount, MWEBV=self.MWEBV, PIXSIZE=self.PIXSIZE) )

        file.write(self.lines[2]+'\n')
        file.write(self.lines[3]+'\n')

        for i in range(len(self.obs_rows)):

            vals = self.obs_rows[i].split(' ')
            MJD = float(vals[1])
            file.write(self.obs_rows[i]+'\n')

        file.write('END_LIBID: {}\n'.format(libid))


class SIMLIB_OBJ:

    def __init__(self, header):

        self.header = header

    def getObservation(self, data):
        
       obs = SIMLIB_OBS(data)
       return obs

    def writeCustomHeader(self, file, NLIBID):

        custom_header = ('DOCUMENTATION:\n'
                        'PURPOSE: simulate LSST based on mock opsim version baseline_v2.0_10yrs. Subsampled for KN. Modified by Ved Shah\n'
                        'INTENT:   Nominal\n'
                        'USAGE_KEY: SIMLIB_FILE\n'
                        'USAGE_CODE: snlc_sim.exe\n'
                        'VALIDATION_SCIENCE:\n'
                        'FIELD: WFD\n'
                        'NOTES:\n'
                        '    PARAMS MINMJD: 60218.0018\n'
                        '    PARAMS MAXMJD: 63870.1061\n'
                        '    PARAMS TOTAL_AREA: 2622.793\n'
                        '    PARAMS SOLID_ANGLE: 0.799\n'
                        'VERSIONS:\n'
                        '- DATE : 21-11-22\n'
                        'AUTHORS : catarina, OpSimSummary version 3.0.0.dev1\n'
                        'DOCUMENTATION_END:\n\n\n\n'
                        'SURVEY: LSST    FILTERS: ugrizY  TELESCOPE: LSST\n'
                        'USER: catarina     HOST: b\'Catarinas-MacBook-Pro.local\'\n'
                        'NLIBID: {NLIBID}\n'
                        'NPE_PIXEL_SATURATE:   100000\n'
                        'PHOTFLAG_SATURATE:    1024\n\n\n'
                        'BEGIN LIBGEN\n').format(NLIBID = NLIBID)

        file.write(custom_header)