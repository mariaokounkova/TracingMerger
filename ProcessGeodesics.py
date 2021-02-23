import numpy as np
import h5py
import sys
import os
import scipy
from scipy import integrate
import shutil

def ReadGeodesicData(p, t_start, t_end):
    """ Read in an array of times and positions for all geodesics at once, 
        and return the trajectories indexed by geodesic """
        
    def AppendGeodesicsTime(Lev, Run):

        print("Reading Geodesic data for " + Lev + " and " + Run)

        file = p + '/' + Lev + '/' + Run + '/Run/Node0.h5'
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
    
        X = [ [] for _ in range(N_geodesics)]
        Y = [ [] for _ in range(N_geodesics)]
        Z = [ [] for _ in range(N_geodesics)]
        L = [ [] for _ in range(N_geodesics)]
        T = [ [] for _ in range(N_geodesics)]
    
        print("Going to go through the keys and the times")
        for k, t in zip(keys, times):                      
            print(k, t)
            if ((t > t_start) and (t < t_end)):
                print("%.1f  " % t, end = '')
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
        for a in range(len(T)):
            if (len(T[a]) > 1):
                ## Remember to add in the minimum index since the starting 
                ## geodesic index just gets incremented during reach refinement
                ## level (by the number of geodesics that came from the levels before)
                ff = open(p + '/Trajectories/' + str(a + m) + '.dat','ab')
                np.savetxt(ff, np.c_[T[a][2:],X[a][2:],Y[a][2:],Z[a][2:],L[a][2:]])
                ff.close()
        print('Finished writing the files')
            
    ## Go through the refinement levels and the segments
    RefinementLevs = [el for el in os.listdir(p) if "Lev" in el]
    ##RefinementLevs = ["Lev_AA", "Lev_AB", "Lev_AC"]
    print("RefinementLevs:", RefinementLevs)
    for lev in RefinementLevs:

        Segments = [el for el in os.listdir(p + '/' + lev) if "Run" in el] 
        print(lev + " Segments:", Segments)
        for segment in Segments:
            AppendGeodesicsTime(lev, segment)

def MakeGeodesicDatFiles(p, t_start, t_end):
    """ Print the result of ReadGeodesicData to files """
    ReadGeodesicData(p, t_start, t_end)

def GetGeodesicTrajectory(p, n):
  """ Read in the post-processed trajectory for the nth geodesic """
  f = p + '/Trajectories/' + str(n) + '.dat'
  t, x, y, z, lapse = np.loadtxt(f, comments="#",usecols=([0,1,2,3,4]),unpack=True)
  return t, x, y, z, lapse

def GetGeodesicIndices(p):
    """ Return the indices of all of the geodesics we have printed to file """
    Files = os.listdir(p + '/Trajectories')
    Indices = [int(file.split('.dat')[0]) for file in Files]
    Indices = sorted(Indices)
    return Indices

def ComputeZeroCrossings(p, n):
    """ Compute the number of zero crossings in the y-z plane - 
        useful for head-on runs where the collision happens
        along the x axis """
    t, x, y, z, lapse = GetGeodesicTrajectory(p, n)
    theta = np.arctan2(z, y)
    zero_crossings = len(np.where(np.diff(np.sign(theta)))[0])
    return zero_crossings

def MakeZeroCrossingsFile(p):
    """ Make a file with the format [geodesic index, number of zero crossings] so 
        that we only have to compute the number of zero crossings once """
    ns = GetGeodesicIndices(p)
    crossings = [ComputeZeroCrossings(p, n) for n in ns]
    np.savetxt(p + '/ZeroCrossings.dat', np.c_[ns, crossings], fmt = '%d %d')
    
def GetGeodesicsZeroCrossingsIndices(p, N):
    """ Return the indices of the geodesics that make N zero-crossings """
    f = p + '/ZeroCrossings.dat'
    ns, zero_crossings = np.loadtxt(f, comments="#",usecols=([0,1]),unpack=True,dtype=int)
    indices = ns[np.where(zero_crossings == N)[0]]
    return indices

def GetGeodesicsZeroCrossings(p):
    """ Return the zero crossings for each index """
    f = p + 'ZeroCrossings.dat'
    ns, zero_crossings = np.loadtxt(f, comments="#",usecols=([0,1]),unpack=True,dtype=int)
    return ns, zero_crossings

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
        t, x, y, z, lapse = GetGeodesicTrajectory(p, n)
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
        t, x, y, z, lapse = GetGeodesicTrajectory(p, n)
        time, kappa = GetFrenetSerretCurvature(p, n)
        index = np.argmax(kappa)
        MaxKappa[i] = kappa[index]
        T[i] = t[index]
        X[i] = x[index]
        Y[i] = y[index]
        Z[i] = z[index]
    
    np.savetxt(p + '/MaxCurvatures.dat', np.c_[ns, MaxKappa, T, X, Y, Z], fmt = '%d %f %f %f %f %f')

def GetFrenetSerretMaxCurvatures(p, loc = False):
    """ Read in the post-processed max curvature for all geodesics
        if loc == True, then also return the location"""
    f = p + '/MaxCurvatures.dat'
    ns, kappa, t, x, y, z = np.loadtxt(f, comments="#",usecols=([0,1,2,3,4,5]),unpack=True)
    if loc == False:
        return ns.astype(int), kappa
    else:
        return ns.astype(int), kappa, x, y, z
    
def GetFrenetSerretMaxCurvature(p, n, loc = False):
    """ Read in the post-processed max curvature for the nth geodesic 
        if loc == True, then also return the location"""
    f = p + '/MaxCurvatures.dat'
    ns, kappa, t, x, y, z = np.loadtxt(f, comments="#",usecols=([0,1,2,3,4,5]),unpack=True)
    index = np.where(ns == n)[0][0]
    print(index)
    if loc == False:
        return kappa[index]
    else:
        return kappa[index], x[index], y[index], z[index]

def ZeroCrossingsFilesIndices(p):
    """ Return the ids of the geodesics that have an interesting 
        number of zero-crossings """
    Indices = []
    for n in range(4, 12):
        print("Copying geodesics with " + str(n) + " zero crossings")
        idx = GetGeodesicsZeroCrossingsIndices(p, n)
        print(str(len(idx)) + " such geodesics")
        Indices = np.concatenate((Indices, idx), axis = 0)
    Indices = Indices.astype(int)
    print(Indices)
    return Indices

def CopyIndicesFiles(p):
    """ Copy over the files that have a set number of zero crossings into
    a new directory (this is useful if we want to delete the files that
    have only one zero crossing in order to save space) """
    Indices = ZeroCrossingsFilesIndices(p)
    for n in Indices:
        shutil.copy(p + '/Trajectories/' + str(n) + '.dat', p + '/Trajectories_Temp/' + str(n) + '.dat')
        #shutil.copy(p + '/FrenetSerret/' + str(n) + '.dat', p + '/FrenetSerret_Temp/' + str(n) + '.dat')

def main():

    p = sys.argv[1]

    ## Make required directories
    try:
        os.mkdir(p + '/Trajectories')
    except:
        print("Trajectory directory already exists")

    #try:
    #    os.mkdir(p + '/FrenetSerret')
    #except:
    #    print("FrenetSerret directory already exists")

    print("Processing geodesics for " + p)
    ## Make geodesic dat files
    print("Making the geodesic dat files")
    #MakeGeodesicDatFiles(p, -1, 1000)
    MakeGeodesicDatFiles(p, 0, 1000)

    ## Compute zero crossings
    print("Computing zero crossings")
    MakeZeroCrossingsFile(p)

    ## Compute Frenet-Serret
    #print("Computing Frenet-Serret")
    #MakeFrenetSerretDatFiles(p)

    ## Compute max curvature file
    #print("Computing max curvatures")
    #MakeMaxCurvatureFile(p)

    ## Now copy over files based on the zero crossings
    print("Copying zero crossings files")

    try:
        os.mkdir(p + '/Trajectories_Temp')
    except:
        print("Trajectories_Temp directory already exists")
    #try:
    #    os.mkdir(p + '/FrenetSerret_Temp')
    #except:
    #    print("FrenetSerret_Temp directory already exists")
    
    CopyIndicesFiles(p)

    ## Now swap the dirs with all of the trajectories with 
    ## the directory with just the interesting trajectories
    shutil.move(p + '/Trajectories', p + '/Trajectories_All')
    shutil.move(p + '/Trajectories_Temp', p + '/Trajectories')

    #shutil.move(p + '/FrenetSerret', p + '/FrenetSerret_All')
    #shutil.move(p + '/FrenetSerret_Temp', p + '/FrenetSerret')

    ## Now remove the _All directories to save space
    #shutil.rmtree(p + '/Trajectories_All')
    #shutil.rmtree(p + '/FrenetSerret_All')


    ## Now clean up the files based on the max curvature

if __name__ == "__main__":
    main()




