import numpy as np
import h5py
import sys
import os
import scipy
from scipy import integrate
import shutil
import itertools
import argparse


## Camera geodesic evaluations

def GetCameraData(p):
    """ Get data about geodesic positions and fates
        in the plane of the camera """
    
    file = p + '/RefinementMethodData.h5'
    f = h5py.File(file, 'r')
    data = f['LensingCore.dat']
    
    tags = data[:,0]
    tags = tags.astype(int)
    x_pos = data[:,1]
    y_pos = data[:,2]
    surface = data[:,4]
    
    return tags, x_pos, y_pos, surface

def GrabSurfaceIndices(surface, fate):
    
    i = np.where(surface == fate)[0]
    return i

def SelectCameraGeodesics(p, fate):
    
    tags, x_pos, y_pos, surface = GetCameraData(p)
    
    i = np.where(surface == fate)[0]
    
    return tags[i], x_pos[i], y_pos[i], surface[i]
    
def GetInfinityGeodesics(p):
    """ 
    inf = np.where(surface == 1.)[0]
    infty = np.where(surface == 7.)[0]
    aha = np.where(surface == 2.)[0]
    ahb = np.where(surface == 3.)[0]
    ahc = np.where(surface == 4.)[0]
    """
    file = p + '/RefinementMethodData.h5'
    f = h5py.File(file, 'r')
    data = f['LensingCore.dat']
    surface = data[:,4]
    return GrabSurfaceIndices(surface, 7)

def GetCameraPosition(p):
    p = p.split('Trace')[1]
    x = p.split('_')[1]
    y = p.split('_')[2]
    z = p.split('_')[3]
    return '[' + x + ',' + y + ',' + z + ']'

def GetTime(p):
    p = p.split('Trace')[1]
    t = p.split('_')[4]
    t = t.split('/')[0]
    return t

## Geodesic data processing

def ReadGeodesicData(p, t_start, t_end):
    """ Read in an array of times and positions for all geodesics at once, 
        and return the trajectories indexed by geodesic """
        
    def AppendGeodesicsTime(Lev):

        print("\n Reading Geodesic data for " + Lev)

        ## Use the first segement to grab the number of geodesics 
        file = p + '/' + Lev + '/Run_AA/Run/Node0.h5'
        f = h5py.File(file, 'r')
        ## grab the .dat files
        keys = [k for k in f.keys() if 'dat' in k]
        ## Array of times from the .dat files
        times = [float(k.split('.dat')[0]) for k in keys]
        ## sort keys according to times
        times, keys = zip(*sorted(zip(times, keys)))
        # grab the number of geodesics
        N_geodesics = len(f[keys[-1]][:,0])
        print("Total geodesics: ", N_geodesics, "Time steps: ", len(times))
        ## Minimum index
        m = int(f[keys[-1]][:,0][0])
        print("Geodesic index offset of this refinement iteration: ", m)

        ## Now make the arrays 
    
        X = [ [] for _ in range(N_geodesics)]
        Y = [ [] for _ in range(N_geodesics)]
        Z = [ [] for _ in range(N_geodesics)]
        L = [ [] for _ in range(N_geodesics)]
        T = [ [] for _ in range(N_geodesics)]
    
        Segments = [el for el in os.listdir(p + '/' + lev) if "Run" in el][::-1]
        print(lev + " Segments:", Segments)

        for segment in Segments:
                
            print("\n Working on segment " + segment)

            file = p + '/' + Lev + '/' + segment + '/Run/Node0.h5'
            f = h5py.File(file, 'r')
            ## grab the .dat files
            keys = [k for k in f.keys() if 'dat' in k]
            ## Array of times from the .dat files
            times = [float(k.split('.dat')[0]) for k in keys]
            ## sort keys according to times
            times, keys = zip(*sorted(zip(times, keys)))

            print("Going to go through the keys and the times")
            for k, t in zip(keys, times): 
                ## Check that time is within the bounds
                if ((t > t_start) and (t < t_end)):
                    ## Check that time is printed at the correct interval
                    if ((t % 0.1 == 0.0) or (abs(t % 0.1 - 0.09999999999999) < 1e-9)):
                        print("%.1f  " % t, end = '', flush=True)
                        data = f[k]
                        ## indices and positions for all geodesics at this time
                        indices = data[:,0]
                        if (len(data[0]) == 7):
                            """ If we've run the numerical evolution so that we haven't 
                                evolved lapse p0 -- just set lapsep0 to 1 in the files. 
                                This field need to be in the resulting .dat file for the 
                                paraview trajectory visualization"""
                            x = data[:,4]
                            y = data[:,5]
                            z = data[:,6]
                            l = np.ones(len(x))
                        elif (len(data[0]) == 8):
                            """ If we have evolved lapse p0"""
                            l = data[:,1]
                            x = data[:,5]
                            y = data[:,6]
                            z = data[:,7]
                        else:
                            print("Unrecognized number of variables")
                        ## fill in the array for each index
                        for i, j in zip(indices.astype(int), range(len(indices))):
                            X[i-m] = np.append(X[i-m], x[j])
                            Y[i-m] = np.append(Y[i-m], y[j])
                            Z[i-m] = np.append(Z[i-m], z[j])
                            L[i-m] = np.append(L[i-m], l[j])
                            T[i-m] = np.append(T[i-m], t)

        print("\n")
        ## Once we have the arrays constructed, append them to the file
        print('Read the geodesic data, now writing the files')
        
        hf = h5py.File(p + '/Trajectories.h5', 'a')
        print(p)
        
        for a in range(len(T)):
            if (len(T[a]) > 1): 
                my_data = np.c_[T[a],X[a],Y[a],Z[a],L[a]]
                ## Remember to add in the minimum index since the starting 
                ## geodesic index just gets incremented during reach refinement
                ## level (by the number of geodesics that came from the levels before)
                ## Hence a + m 
                hf.create_dataset('Geodesic' + str(a+m) + '.dat', data=my_data)
                
        hf.close()
        
        ## Previous code for writing .dat files
        #for a in range(len(T)):
        #    if (len(T[a]) > 1):
        #        ## Remember to add in the minimum index since the starting 
        #        ## geodesic index just gets incremented during reach refinement
        #        ## level (by the number of geodesics that came from the levels before)
        #        ff = open(p + '/Trajectories/' + str(a + m) + '.dat','ab')
        #        np.savetxt(ff, np.c_[T[a][2:],X[a][2:],Y[a][2:],Z[a][2:],L[a][2:]])
        #        ff.close()
        
        print('Finished writing the files')
            
    ## Go through the refinement levels and the segments
    RefinementLevs = [el for el in os.listdir(p) if "Lev" in el]
    print("RefinementLevs:", RefinementLevs)
    for lev in RefinementLevs:

        AppendGeodesicsTime(lev)

def MakeGeodesicDatFiles(p, t_start, t_end):
    """ Print the result of ReadGeodesicData to files """
    ReadGeodesicData(p, t_start, t_end)
    
def GetGeodesicTrajectory(p, n):
    """ Read in the post-processed trajectory for the nth geodesic """
    f = p + '/Trajectories.h5'
    hf = h5py.File(f, 'r')
    data = hf['Geodesic' + str(n) + '.dat']
    t = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]
    hf.close()
    ## Make sure the times are sorted
    t, x, y, z = zip(*sorted(zip(t, x, y, z)))
    t = np.array(t)
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    # Former method of reading in 
    #t, x, y, z, lapse = np.loadtxt(f, comments="#",usecols=([0,1,2,3,4]),unpack=True)
    return t, x, y, z

def GetGeodesicIndices(p, infinity=True):
    """ Return the indices of all of the geodesics we have printed to file 
        If infinity == True, only returnt the geodesics that make it out to infinity """
    Files = os.listdir(p + '/Trajectories')
    f = p + '/Trajectories.h5'
    hf = h5py.File(f, 'r')
    Indices = [int(k.split('Geodesic')[1].split('.dat')[0]) for k in hf.keys()]
    Indices = sorted(Indices)
    if infinity:
        indices_infinity = GetInfinityGeodesics(p)
        Indices = list(set(indices_infinity) & set(Indices))
    return Indices

## Zero crossings computations

def ComputeZeroCrossings(p, n):
    """ Compute the number of zero crossings in the y-z plane - 
        useful for head-on runs where the collision happens
        along the x axis """
    t, x, y, z = GetGeodesicTrajectory(p, n)
    theta = np.arctan2(z, y)
    zero_crossings = len(np.where(np.diff(np.sign(theta)))[0])
    return zero_crossings

def MakeZeroCrossingsFile(p):
    """ Make a file with the format [geodesic index, number of zero crossings] so 
        that we only have to compute the number of zero crossings once """
    ns = GetGeodesicIndices(p)
    crossings = [ComputeZeroCrossings(p, n) for n in ns]
    np.savetxt(p + 'ZeroCrossings.dat', np.c_[ns, crossings], fmt = '%d %d')
    
def GetGeodesicsZeroCrossings(p):
    """ Return the zero crossings for each index """
    f = p + 'ZeroCrossings.dat'
    ns, zero_crossings = np.loadtxt(f, comments="#",usecols=([0,1]),unpack=True,dtype=int)
    return ns, zero_crossings

def GetGeodesicsZeroCrossingsIndices(p, N, infinity=True):
    """ Return the indices of the geodesics that make N zero-crossings 
        If infinity == True, return only the ones that make it to infinity"""
    ns, zero_crossings = GetGeodesicsZeroCrossings(p)
    indices = ns[np.where(zero_crossings == N)[0]]
    if infinity:
        indices_infinity = GetInfinityGeodesics(p)
        indices = list(set(indices_infinity) & set(indices))
    return indices

def GetGeodesicsZeroCrossingsIndicesGreater(p, N, infinity=True):
    """ Return the indices of the geodesics that make N or more zero-crossings 
        If infinity == True, return only the ones that make it to infinity"""
    ns, zero_crossings = GetGeodesicsZeroCrossings(p)
    indices = ns[np.where((zero_crossings >= N))[0]]
    if infinity:
        indices_infinity = GetInfinityGeodesics(p)
        indices = list(set(indices_infinity) & set(indices))
    return indices

def GetGeodesicsMaxZeroCrossings(p, infinity=True):
    """ Return the maximum number N of zero crossings for a given run, 
        along with the indices of the geodesics that make N zero-crossings 
        If infinity == Ture, only return the ones that make it to infinity"""
    ns, zero_crossings = GetGeodesicsZeroCrossings(p)
    max_N = max(zero_crossings)
    indices = GetGeodesicsZeroCrossingsIndices(p, max_N, infinity=infinity)
    return max_N, indices


## X Turns computations 

def ComputeXTurns(p, n):
    """ Compute tbe number of turns a geodesic makes along the x direction
        Useful for head-on runs where the collision happens along the x axis"""
    print(n, end=' ')
    t, x, y, z = GetGeodesicTrajectory(p, n)
    x_diff = x[1:] - x[:-1]
    turns = len(list(itertools.groupby(x_diff, lambda x_diff: x_diff > 0)))
    ## Set to zero turns if the trajectory does not change in x - these end up being
    ## artifical turns hovering around zero (and unresolved)
    if np.max(np.abs(x)) < 1e-3:
        turns = 0
    return turns 

def MakeXTurnsFile(p):
    """ Make a file with the format [geodesic index, number of X turns] so 
        that we only have to compute the number of X turns once """
    ns = GetGeodesicIndices(p)
    print('TOTAL NUMBER: ', str(len(ns)))
    turns = [ComputeXTurns(p, n) for n in ns]
    np.savetxt(p + 'XTurns.dat', np.c_[ns, turns], fmt = '%d %d')
    
def GetGeodesicsXTurns(p):
    """ Return the x turns for each index """
    f = p + 'XTurns.dat'
    ns, turns = np.loadtxt(f, comments="#",usecols=([0,1]),unpack=True,dtype=int)
    return ns, turns

def GetGeodesicsXTurnsIndices(p, N, infinity=True):
    """ Return the indices of the geodesics that make N x turns 
        If infinity == True, return only the ones that make it to infinity"""
    ns, turns = GetGeodesicsXTurns(p)
    indices = ns[np.where(turns == N)[0]]
    if infinity:
        indices_infinity = GetInfinityGeodesics(p)
        indices = list(set(indices_infinity) & set(indices))
    return indices

def GetGeodesicsXTurnsIndicesGreater(p, N, infinity=True, N_less = 1e6):
    """ Return the indices of the geodesics that make N or more x turns 
        If infinity == True, return only the ones that make it to infinity"""
    ns, turns = GetGeodesicsXTurns(p)
    indices = ns[np.where(turns >= N)[0]]
    indices_less = ns[np.where(turns <= N_less)[0]]
    indices = list(set(indices) & set(indices_less))
    if infinity:
        indices_infinity = GetInfinityGeodesics(p)
        indices = list(set(indices_infinity) & set(indices))
    return indices
    
def GetGeodesicsMaxXTurns(p, infinity=True):
    """ Return the maximum number N of zero crossings for a given run, 
        along with the indices of the geodesics that make N zero-crossings 
        If infinity == Ture, only return the ones that make it to infinity"""
    f = p + 'ZeroCrossings.dat'
    ns, turns = GetGeodesicsXTurns(p)
    max_N = max(turns)
    indices = GetGeodesicsXTurnsIndices(p, max_N, infinity=infinity)
    return max_N, indices

## Frenet-Serret computations

def ComputeArcLength(t, x, y, z):
    """ For a given trajectory {x(t), y(t), z(t)}, compute and return
        the arclength s(t) as a function of time"""
    dx = np.gradient(x, t)
    dy = np.gradient(y, t)
    dz = np.gradient(z, t)

    ## Compute the magnitude
    mag = np.sqrt(dx**2 + dy**2 + dz**2)
    
    ## Take a running integral 
    s = integrate.cumtrapz(mag, t, initial = 0.0)
    
    ## check that s(t) is monotonically increasing
    check = np.where((s[1:] - s[:-1]) <= 0.0)[0]
    if len(check) > 0:
        print("s(t) is non-monotonically increasing!")
 
    return s

def ComputeFrenetSerretTNB(t, x, y, z):
    s = ComputeArcLength(t, x, y, z)
    
    T_x = np.gradient(x, s)
    T_y = np.gradient(y, s)
    T_z = np.gradient(z, s)

    T_mag = np.sqrt(T_x**2 + T_y**2 + T_z**2)

    T_x = T_x / T_mag
    T_y = T_y / T_mag
    T_z = T_z / T_mag

    dTds_x = np.gradient(T_x, s)
    dTds_y = np.gradient(T_y, s)
    dTds_z = np.gradient(T_z, s)
    
    dTds_mag = np.sqrt(dTds_x**2 + dTds_y**2 + dTds_z**2)
    
    N_x = dTds_x / dTds_mag
    N_y = dTds_y / dTds_mag
    N_z = dTds_z / dTds_mag
    
    B_x = T_y * N_z - T_z * N_y
    B_y = - T_x * N_z + T_z * N_x
    B_z = T_x * N_y - T_y * N_x
    
    kappa = dTds_mag ## curvature

    return T_x, T_y, T_z, N_x, N_y, N_z, B_x, B_y, B_z, kappa

def MakeFrenetSerretDatFiles(p):
    """ For all geodesics .dat files in a given directory, dump the Frenet-Serret Frame 
        data (so we only have to compute it once)"""
    ns = GetGeodesicIndices(p)
    for n in ns: 
        t, x, y, z = GetGeodesicTrajectory(p, n)
        t, x, y, z = zip(*sorted(zip(t, x, y, z)))
        T_x, T_y, T_z, N_x, N_y, N_z, B_x, B_y, B_z, kappa = ComputeFrenetSerretTNB(t, x, y, z)
        np.savetxt(p + '/FrenetSerret/' + str(n) + '.dat', np.c_[t, kappa, T_x, T_y, T_z, \
                                                                        N_x, N_y, N_z, \
                                                                        B_x, B_y, B_z])
        
def GetFrenetSerretVectors(p, n):
    """ Read in the post-processed trajectory for the nth geodesic """
    f = p + '/FrenetSerret/' + str(n) + '.dat'
    t, kappa, T_x, T_y, T_z, N_x, N_y, N_z, B_x, B_y, B_z \
        = np.loadtxt(f, comments="#",usecols=(range(11)),unpack=True)
    return t, kappa, T_x, T_y, T_z, N_x, N_y, N_z, B_x, B_y, B_z
  
def GetFrenetSerretCurvature(p, n):
    """ Read in the post-processed curvature for the nth geodesic """
    f = p + '/FrenetSerret/' + str(n) + '.dat'
    t, kappa = np.loadtxt(f, comments="#",usecols=([0,1]),unpack=True)
    return t, kappa

def MakeMaxCurvatureFile(p):
    """ Make file constraining the geodesic ID, the maximum value of the curvature, 
        and the x, y, z, coordinates of this point"""
    ns = GetGeodesicIndices(p)
    N = len(ns)
    MaxKappa = np.zeros(N)
    T = np.zeros(N)
    X = np.zeros(N)
    Y = np.zeros(N)
    Z = np.zeros(N)
    
    for i in range(N):
        n = ns[i]
        t, x, y, z = GetGeodesicTrajectory(p, n)
        time, kappa = GetFrenetSerretCurvature(p, n)
        index = np.argmax(kappa)
        MaxKappa[i] = kappa[index]
        T[i] = t[index]
        X[i] = x[index]
        Y[i] = y[index]
        Z[i] = z[index]
    
    np.savetxt(p + '/MaxCurvatures.dat', np.c_[ns, MaxKappa, T, X, Y, Z], fmt = '%d %f %f %f %f %f')

def main():

    p = argparse.ArgumentParser(description="Post-process geodesics evolution data")
    p.add_argument("--dir", required=True, \
           help="Root directory of geodesic evolution data to post-process")
    p.add_argument("--trajectories", help='Post-process the trajectories', \
        dest='trajectories', action='store_true')
    p.add_argument("--zerocrossings", help='Post-process the zero-crossings', \
        dest='zerocrossings', action='store_true')
    p.add_argument("--xturns", help='Post-process the x-turns', \
        dest='xturns', action='store_true')
    p.add_argument("--frenetserret", help='Perform the Frenet-Serret analysis', \
        dest='frenetserret', action='store_true')
    
    p.set_defaults(trajectories=False)
    p.set_defaults(zerocrossings=False)
    p.set_defaults(xturns=False)
    p.set_defaults(frenetserret=False)
   
    args = p.parse_args()

    print("Processing geodesics for " + args.dir)
    
    ## Make geodesic dat files
    if (args.trajectories):
        print("Making the geodesic trajectory file")
        try:
            os.remove(args.dir + '/Trajectories.h5')
            print("Removed previous Trajectories.h5 file")
        except:
            print("Trajectories.h5 file did not previous exist")
        MakeGeodesicDatFiles(args.dir, 0, 1000)

    ## Compute zero crossings
    if (args.zerocrossings):
        print("Computing zero crossings")
        MakeZeroCrossingsFile(args.dir)
    
    ## Compute x turns
    if (args.xturns):
        print("Computing x turns")
        MakeXTurnsFile(args.dir)

    ## Compute Frenet-Serret
    if (args.frenetserret):
        print("Computing Frenet-Serret")
        try:
            os.mkdir(p + '/FrenetSerret')
        except:
            print("FrenetSerret directory already exists")
        MakeFrenetSerretDatFiles(p)
        print("Computing max curvatures")
        MakeMaxCurvatureFile(p)

if __name__ == "__main__":
    main()




