import numpy as np
import astropy
from astropy import units as u
import astropy.coordinates as coord
import astropy.io.ascii 
import math

FLARE_DATA_PATH = 'flare_data/apjaa8ea2t3_mrt.txt'

TOTAL_KEPLER_M_DWARF_COUNT = 4664
TOTAL_KEPLER_DURATION_IN_DAYS = 1460
TOTAL_FLARE_INSTANCES_COUNT = 103187

LOCAL_RADIUS = 20
NUMBER_OF_LOCAL_M_DWARFS = 1082

KEPLER_MEAN_EFFECTIVE_TEMP_FOR_M_DWARFS = 3743.4117647058824
KEPLER_STD_EFFECTIVE_TEMP_FOR_M_DWARFS = 161.37182827551771



def get_flare_number_of_expected_flares(radius, duration):

    flare_density = (NUMBER_OF_LOCAL_M_DWARFS / (4/3 * math.pi * LOCAL_RADIUS**2)) * (TOTAL_FLARE_INSTANCES_COUNT / (TOTAL_KEPLER_M_DWARF_COUNT * TOTAL_KEPLER_DURATION_IN_DAYS))
    volume = 4/3 * math.pi * (radius ** 3)
    flare_count = math.floor(flare_density * volume * duration)

    print('Flare Density:', flare_density, 'flares/pc^3/day')
    print('Expected Number of flares in {radius} pc sphere during the {duration} day period:'.format(radius = radius, duration = duration), flare_count)

    return flare_count

def get_uniformly_distributed_spherical_coordinates(radius, count):

    R = 1000 # Big Box edge in pc
    ra = []
    dec = []
    dist = []

    while len(ra) != count:
        x = (np.random.uniform(-1 * R,R) * u.pc)
        y = (np.random.uniform(-1 * R,R) * u.pc)
        z = (np.random.uniform(-1 * R,R) * u.pc)
        if (x ** 2 + y ** 2 + z ** 2) ** (0.5) < radius * u.pc:
            r, lat, lon = coord.cartesian_to_spherical(x,y,z)
            ra.append(lon.degree)
            dec.append(lat.degree)
            dist.append(r.value)


def get_normally_distributed_star_temp(count):
    return np.random.normal(KEPLER_MEAN_EFFECTIVE_TEMP_FOR_M_DWARFS, KEPLER_STD_EFFECTIVE_TEMP_FOR_M_DWARFS, count)

    
get_flare_number_of_expected_flares(20, 3650)