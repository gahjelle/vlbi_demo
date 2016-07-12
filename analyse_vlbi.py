from datetime import datetime
import numpy as np

import coords

# Constants
c_light = 299792458


def analyse_ngs_file(filename):

    # Earth orientation parameters for Aug 4 2015
    eop = dict(x=0.202223, y=0.414629, dx=0.000181, dy=-0.000178)

    stations = dict()
    sources = dict()
    residual = list()

    with open(filename) as fid:
        # Skip header lines
        next(fid)
        next(fid)

        # Stations
        for line in fid:
            fields = line.split()
            name = fields[0]
            if name == '$END':
                break
            stations[name] = np.array((float(fields[1]), float(fields[2]), float(fields[3])))

        # Sources
        for line in fid:
            fields = line.split()
            name = fields[0]
            if name == '$END':
                break
            right_ascension = coords.hms_to_rad(fields[1], fields[2], fields[3])
            declination = coords.dms_to_rad(line[29:32].replace(' ', ''), line[33:35], line[36:48])
            print(line[29:48], declination)
            sources[name] = np.array((np.cos(declination) * np.cos(right_ascension),
                                      np.cos(declination) * np.sin(right_ascension),
                                      np.sin(declination)))

        # Auxilliary info thrown away
        for line in fid:
            if line.startswith('$END'):
                break

        # Observations
        for line in fid:
            fields = line.strip().split()

            # Calculate delay
            if fields[-1].endswith('01'):
                sta_1 = fields[0]
                sta_2 = fields[1]
                src = fields[2]
                obsdt = datetime(int(fields[3]), int(fields[4]), int(fields[5]),
                                 int(fields[6]), int(fields[7]), int(float(fields[8])))

                K = sources[src]
                b = np.dot(coords.itrs_to_gcrs(obsdt, eop), (stations[sta_2] - stations[sta_1]))
                theoretical_delay = -np.dot(K, b)


            if fields[-1].endswith('02'):
                observed_delay = float(fields[0]) * c_light * 1e-9
                residual.append(observed_delay - theoretical_delay)

        return residual


residual = analyse_ngs_file('15AUG04XA_N004')

import matplotlib.pyplot as plt
plt.plot(residual, 'x')
plt.show()
