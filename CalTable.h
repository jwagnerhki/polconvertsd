/* CALTABLE - calibration table interface to PolConvert

             Copyright (C) 2013  Ivan Marti-Vidal
             Nordic Node of EU ALMA Regional Center (Onsala, Sweden)
             Max-Planck-Institut fuer Radioastronomie (Bonn, Germany)
  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
  
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
  
You should have received a copy of the GNU General Public License   
along with this program.  If not, see <http://www.gnu.org/licenses/>
  
*/



#include <sys/types.h>
#include <iostream>
#include <complex> 

/* Class to read a calibration table and interpolate
   its gains in frequency and time. */
class CalTable {
  public:
     CalTable(int kind, double **R1,double **I1,double **R2,double **I2, double *freqs, double **times, int Na, long *Nt, long Nc, bool **flag, bool islinear, FILE *logF);
     ~CalTable();
     int getNant() const;
     long getNchan() const;
     bool isBandpass() const;
     void getTimeRange(double *JD) const;
     void getFreqRange(double *Fr) const;
     long getNEntries(int ant) const;
     void getFrequencies(double *freqs) const;
     void getTimes(int ant, double *times) const;
     void getGains(int ant, long timeidx, double *gain[4]) const;
     bool getChanged() const;
     void setChanged(bool ch);

/* Prepares the instance for the frequency interpolation. User must provide
   the array of frequencies to interpolate (*freqs) and the number 
   of channels in that array (mschan). */
     void setMapping(long mschan, double *freqs);

/* Prepares the instance for the time interpolation. Returns False if the interpolation coefficients have not changed (so it would be a waste of resources to recompute everything). */
     bool setInterpolationTime(double itime);

/* Interpolates the gains of an antenna (iant) at a given time, itime, 
   and applies them to the array gain[2] (elements of this array
   are the different polarizations (X,Y)). User should have run "setMapping" before,
   in case that the table channel frequencies do not match with the array of frequencies
   where the interpolation is desired. Then,
   If mode==0, the gains are just written to gain. 
   If mode==1, the gains are ADDED to the already-existing values in gain.
   If mode==2, the gains are MULTIPLIED to the already-existing values in gain. 
 -------------------
*/
     void applyInterpolation(int iant, int mode, std::complex<float>* gain[2]);

// Same as above, but one gain (for one channel) is returned:
     bool getInterpolation(int iant, int ichan, std::complex<float> gain[2]);



  private:

     const int Nants;
     const long Nchan;

     FILE *logFile;
     char message[512];
     void interpolateFlaggedData();  // Fills flagged gains with interpolated values.
     void interpolateFlaggedData_Time(const int ant);
     void interpolateFlaggedData_Channel(const int ant);
     static const int Nmax = 256; // Maximum number of antennas.
     std::string name;
     long *Ntimes;

     bool SignFreq, success;
     double ***GainAmp[2];
     double ***GainPhase[2];
     bool **Flags, *firstTime;
     double *BuffPhase[2];
     double *BuffAmp[2];
     double **Time;
     double currTime;
     double *Freqs;
     double JDRange[2];
     double *K0;
     long *I0;
     long *I1;
     double *preKt;
     long *pret0, *pret1;
     long MSChan;

     bool isDelay, gainChanged, isDterm, isTsys;
     const bool isLinear;

     double deltaNu0, deltaNu;
     std::complex<float>** bufferGain[2];
};




