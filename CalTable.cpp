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

// Refactoring 6/2018 JanW

#include <sys/types.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <complex>

#include "CalTable.h"

#define PI 3.141592653589793  // todo: math.h M_PI, M_2PI ?
#define TWOPI 6.283185307179586

CalTable::~CalTable()
{
    // TODO: delete all the allocated stuff!
}

CalTable::CalTable(
    int kind,
    double **R1, double **P1, double **R2, double **P2,
    double *freqs, double **times,
    int Na, long *Nt, long Nc, bool **flag, bool islinear, FILE *logF)
    : Nants(Na), Nchan(Nc), isLinear(islinear)
{
    // Initiate variables and declare auxiliary variables
    logFile = logF;
    MSChan = Nchan;
    isDelay = (kind==1); // ToDO: enum! CalTableKind { Delay, DTerm, Tsys, ... }
    isDterm = (kind==2);
    isTsys = (kind==3);
    gainChanged = true;
    currTime = -1.0;

    // Array allocations
    firstTime = new bool[Nants];
    for (int ia=0; ia<Nants; ia++) {
        firstTime[ia] = true;
    }

    // Copy input array args
    Ntimes = new long[Nants];
    std::memcpy(Ntimes,Nt,sizeof(long)*Nants);

    Freqs = new double[Nchan];
    std::memcpy(Freqs,freqs,sizeof(double)*Nchan);
    SignFreq = (Freqs[1]>Freqs[0]);

    Flags = new bool*[Nants];
    for (int ia=0; ia<Nants; ia++) {
        Flags[ia] = new bool[Ntimes[ia]*Nchan];
        std::memcpy(Flags[ia],flag[ia],sizeof(bool)*Ntimes[ia]*Nchan);
    }

    Time = new double*[Nants];
    for (int ia=0; ia<Nants; ia++) {
        Time[ia] = new double[Ntimes[ia]];
        std::memcpy(Time[ia],times[ia],sizeof(double)*Ntimes[ia]);
    }

    for (int polnr=0; polnr<=1; polnr++) {
        GainAmp[polnr] = new double**[Nants];
        GainPhase[polnr] = new double**[Nants];
        for (long i=0; i<Nants; i++) {
            GainAmp[polnr][i] = new double*[Nchan];
            GainPhase[polnr][i] = new double*[Nchan];
            for (long j=0; j<Nchan; j++) {
                GainAmp[polnr][i][j] = new double[Ntimes[i]];
                GainPhase[polnr][i][j] = new double[Ntimes[i]];
                for (long k=0; k<Ntimes[i]; k++){
                    long auxI = j*Ntimes[i]+k;
                    GainAmp[polnr][i][j][k] = R1[i][auxI];
                    GainPhase[polnr][i][j][k] = P1[i][auxI];
                }
            }
        }
    }

    // Detect time range (min T, max T)
    JDRange[0] = 1.e20;
    JDRange[1] = 0.0;
    for (long i=0; i<Nants; i++) {
        if (Time[i][0]<JDRange[0]) { JDRange[0] = Time[i][0]; }
        if (Time[i][Ntimes[i]-1]>JDRange[1]) { JDRange[1] = Time[i][Ntimes[i]-1]; }
    }

    // Interpolate across any data gaps
    interpolateFlaggedData();

    // Set default channel mapping to 1->1:
    preKt = new double[Nants];
    pret0 = new long[Nants];
    pret1 = new long[Nants];
    for (long i=0; i<Nants; i++) {
       pret0[i] = 0;
       pret1[i] = 0;
       preKt[i] = -1.0;
    }

    for (int polnr=0; polnr<=1; polnr++) {
        bufferGain[polnr] = new std::complex<float>*[Nants];
        for (long i=0; i<Nants; i++) {
            bufferGain[polnr][i] = new std::complex<float>[Nchan];
        }
    }

    K0 = new double[Nchan];
    I0 = new long[Nchan];
    I1 = new long[Nchan];
    for (long i = 0; i<Nchan; i++) {
      K0[i] = 1.0;
      I0[i] = i;
      I1[i] = 0;
    }

}


bool CalTable::isBandpass() const
{
    return Nchan > 1;
}

int CalTable::getNant() const
{
    return Nants;
}

long CalTable::getNchan() const
{
    return Nchan;
}

long CalTable::getNEntries(int ant) const
{
    assert((ant >= 0) && (ant < Nant));
    return Ntimes[ant];
}

void CalTable::getTimeRange(double *JD) const
{
    assert(JD != NULL);
    JD[0] = JDRange[0];
    JD[1] = JDRange[1];
}

bool CalTable::getChanged() const
{
    return gainChanged;
}

void CalTable::setChanged(bool ch)
{
    gainChanged = ch;
}

void CalTable::getFreqRange(double *Fr) const
{
    assert(Fr != NULL);
    if (Freqs[0] < Freqs[Nchan-1]) {
        Fr[0] = Freqs[0];
        Fr[1] = Freqs[Nchan-1];
    } else {
        Fr[0] = Freqs[Nchan-1];
        Fr[1] = Freqs[0];
    }
}

void CalTable::getFrequencies(double *freqs) const
{
    assert(freqs != NULL);
    std::memcpy(freqs,Freqs,sizeof(double)*Nchan);
}

void CalTable::getTimes(int ant, double *times) const
{
    assert((ant >= 0) && (ant < Nant));
    assert(times != NULL);
    std::memcpy(times,Time[ant],sizeof(double)*Ntimes[ant]);
    // printf("TIMES: %.8e -  %.8e\n",Time[ant][0],times[0]);
}

void CalTable::getGains(int ant, long timeidx, double *gain[4]) const
{
    assert((ant >= 0) && (ant < Nant));
    assert((timeidx >= 0) && (timeidx < Ntimes[ant]));
    for (long i=0; i< Nchan; i++) {
        gain[0][i] = GainAmp[0][ant][i][timeidx];
        gain[1][i] = GainAmp[1][ant][i][timeidx];
        gain[2][i] = GainPhase[0][ant][i][timeidx];
        gain[3][i] = GainPhase[1][ant][i][timeidx];
    }
}

/** Interpolate failed frequency channels of each antenna and for each time */
void CalTable::interpolateFlaggedData()
{
#ifdef CT_DEBUG
    FILE *gainFile = fopen("GAINS.ASSESS.orig","ab");
    fwrite(&Nchan,sizeof(int),1,gainFile);
    for (int chan=0; chan<Nchan; chan++) {
        fwrite(&GainAmp[0][0][chan][0],sizeof(double),1,gainFile);
        fwrite(&GainAmp[1][0][chan][0],sizeof(double),1,gainFile);
        fwrite(&GainPhase[0][0][chan][0],sizeof(double),1,gainFile);
        fwrite(&GainPhase[1][0][chan][0],sizeof(double),1,gainFile);
        fwrite(&Flags[0][chan*Ntimes[0]+0],sizeof(bool),1,gainFile);
    }
    fclose(gainFile);
#endif

    if (Nchan==1) {
        for (int ant=0; ant<Nants; ant++) {
            interpolateFlaggedData_Time(ant);
        }
    } else {
        for (int ant=0; ant<Nants; ant++) {
            interpolateFlaggedData_Channel(ant);
        }
    }

#ifdef CT_DEBUG
    FILE *gainFile2 = fopen("GAINS.ASSESS.interp","ab");
    fwrite(&Nchan,sizeof(int),1,gainFile2);
    for (chan=0; chan<Nchan; chan++) {
        fwrite(&GainAmp[0][0][chan][0],sizeof(double),1,gainFile2);
        fwrite(&GainAmp[1][0][chan][0],sizeof(double),1,gainFile2);
        fwrite(&GainPhase[0][0][chan][0],sizeof(double),1,gainFile2);
        fwrite(&GainPhase[1][0][chan][0],sizeof(double),1,gainFile2);
        fwrite(&Flags[0][chan*Ntimes[0]+0],sizeof(bool),1,gainFile2);
     }
     fclose(gainFile2);
#endif
}

void CalTable::interpolateFlaggedData_Time(const int ant)
{
    bool allflagged = true;
    for (long tidx=0; tidx<Ntimes[ant]; tidx++) {
        if (!Flags[ant][tidx]) {
            allflagged = false;
            break;
        }
    }

    // If all data were flagged, set gains to zero and remove flagging
    if (allflagged) {
        snprintf(message,sizeof(message)-1,"\nWARNING: ALMA ANTENNA #%i HAS ALL ITS TIMES FLAGGED!\n",ant);
        fprintf(logFile,"%s",message);
        snprintf(message,sizeof(message)-1,"SETTING ITS GAIN TO ZERO. CHECK RESULTS CAREFULLY!\n\n");
        fprintf(logFile,"%s",message);
        fflush(logFile);

        for (long tidx=0; tidx<Ntimes[ant]; tidx++) {
            if (isDterm) {
                GainAmp[0][ant][0][tidx] = 0.0;
                GainAmp[1][ant][0][tidx] = 0.0;
            } else {
                GainAmp[0][ant][0][tidx] = 0.0;
                GainAmp[1][ant][0][tidx] = 0.0;
            }
            GainPhase[0][ant][0][tidx] = 0.0;
            GainPhase[1][ant][0][tidx] = 0.0;
            Flags[ant][tidx] = false;
        }
        return;
    }

    // Determine first unflagged time
    long firstIdx = 0;
    while ((firstIdx < Ntimes[ant]) && Flags[ant][firstIdx]) {
        ++firstIdx;
    }
    assert(firstIdx < Ntimes[ant] /* should never happen since at least one non-flagged original or one inserted-zero entry must exist*/);

    // Fill the first flagged times if any. No interpolation.
    for (long index=0; index<firstIdx; index++) {
        GainAmp[0][ant][0][index] = GainAmp[0][ant][0][firstIdx];
        GainAmp[1][ant][0][index] = GainAmp[1][ant][0][firstIdx];
        GainPhase[0][ant][0][index] = GainPhase[0][ant][0][firstIdx];
        GainPhase[1][ant][0][index] = GainPhase[1][ant][0][firstIdx];
        Flags[ant][index] = false;
    }

    // Determine last unflagged time
    long lastIdx = Ntimes[ant]-1;
    while ((lastIdx >= 0) && Flags[ant][lastIdx]) {
        --lastIdx;
    }
    assert(firstIdx <= lastIdx);

    // Fill the last flagged channels if any. No interpolation.
    for (long index=lastIdx+1; index<Ntimes[ant]; index++) {
        GainAmp[0][ant][0][index] = GainAmp[0][ant][0][lastIdx];
        GainAmp[1][ant][0][index] = GainAmp[1][ant][0][lastIdx];
        GainPhase[0][ant][0][index] = GainPhase[0][ant][0][lastIdx];
        GainPhase[1][ant][0][index] = GainPhase[1][ant][0][lastIdx];
        Flags[ant][index] = false;
    }

    // Now look for flagged intermediate time ranges. Linearly interpolate.
    bool firstflag = false;
    for (long tidx=1; tidx<Ntimes[ant]; tidx++) {
        if (!firstflag && Flags[ant][tidx]) {
            firstflag = true;
            firstIdx = tidx-1;
        }
        if (firstflag && !Flags[ant][tidx]) {
            firstflag = false;
            lastIdx = tidx;
            for (long intpIdx=firstIdx+1; intpIdx<lastIdx; intpIdx++) {
                double frchan = ((double) (intpIdx - firstIdx)) / ((double) (lastIdx - firstIdx));
                GainAmp[0][ant][0][intpIdx] = GainAmp[0][ant][0][lastIdx]*frchan + GainAmp[0][ant][0][firstIdx]*(1.0 - frchan);
                GainAmp[1][ant][0][intpIdx] = GainAmp[1][ant][0][lastIdx]*frchan + GainAmp[1][ant][0][firstIdx]*(1.0 - frchan);
                GainPhase[0][ant][0][intpIdx] = GainPhase[0][ant][0][lastIdx]*frchan + GainPhase[0][ant][0][firstIdx]*(1.0 - frchan);
                GainPhase[1][ant][0][intpIdx] = GainPhase[1][ant][0][lastIdx]*frchan + GainPhase[1][ant][0][firstIdx]*(1.0 - frchan);
                Flags[ant][intpIdx] = false;
            }
        }
    }
}

void CalTable::interpolateFlaggedData_Channel(const int ant)
{
    for (long tidx=0; tidx<Ntimes[ant]; tidx++) {

        long offset = tidx;

         bool allflagged = true;
         for (long chan=0; chan<Nchan; chan++) {
            if (!Flags[ant][chan*Ntimes[ant]+offset]) {
                allflagged = false;
                break;
            }
        }

        // If all data were flagged, set gains to zero and remove flagging
        if (allflagged) {
            snprintf(message,sizeof(message)-1,"\nWARNING: ALMA ANTENNA #%i HAS ALL CHANNELS FLAGGED AT TIME #%li\n",ant,tidx);
            fprintf(logFile,"%s",message);
            snprintf(message,sizeof(message)-1,"SETTING ITS GAIN TO DUMMY. CHECK RESULTS CAREFULLY!\n\n");
            fprintf(logFile,"%s",message);
            fflush(logFile);
            for (long chan=0; chan<Nchan; chan++) {
                if(isDterm) {
                    GainAmp[0][ant][chan][tidx] = 0.0;
                    GainAmp[1][ant][chan][tidx] = 0.0;
                } else {
                    GainAmp[0][ant][chan][tidx] = 1.0;
                    GainAmp[1][ant][chan][tidx] = 1.0;
                }
                GainPhase[0][ant][chan][tidx] = 0.0;
                GainPhase[1][ant][chan][tidx] = 0.0;
                Flags[ant][chan*Ntimes[ant]+offset] = false;
            }
            continue;
        }

        // POLARIZATION-WISE

        // Determine first unflagged channel (or first resetted-flag channel)
        long firstIdx = 0;
        while ((firstIdx < Nchan) && Flags[ant][firstIdx*Ntimes[ant]+offset]) {
            ++firstIdx;
        }
        assert(firstIdx < Nchan /* should never happen since at least one non-flagged original or one inserted-zero entry must exist*/);

        // Fill the first flagged channels if any. No interpolation.
        for (long chan=0; chan<firstIdx; chan++){
            GainAmp[0][ant][chan][tidx] = GainAmp[0][ant][firstIdx][tidx];
            GainAmp[1][ant][chan][tidx] = GainAmp[1][ant][firstIdx][tidx];
            GainPhase[0][ant][chan][tidx] = GainPhase[0][ant][firstIdx][tidx];
            GainPhase[1][ant][chan][tidx] = GainPhase[1][ant][firstIdx][tidx];
            Flags[ant][chan*Ntimes[ant]+offset] = false;
        }

        // Determine last unflagged channel
        long lastIdx = Ntimes[ant]-1;
        while ((lastIdx >= 0) && Flags[ant][lastIdx*Ntimes[ant]+offset]) {
            --lastIdx;
        }
        assert(firstIdx <= lastIdx);

        // Fill the last flagged channels
        for (long chan=lastIdx+1; chan<Nchan; chan++) {
            GainAmp[0][ant][chan][tidx] = GainAmp[0][ant][lastIdx][tidx];
            GainAmp[1][ant][chan][tidx] = GainAmp[1][ant][lastIdx][tidx];
            GainPhase[0][ant][chan][tidx] = GainPhase[0][ant][lastIdx][tidx];
            GainPhase[1][ant][chan][tidx] = GainPhase[1][ant][lastIdx][tidx];
            Flags[ant][chan*Ntimes[ant]+offset] = false;
        }

        // Look for flagged regions within the band. Linearly interpolate.
        bool firstflag = false;
        for (long chan=1; chan<Nchan; chan++) {
            if (!firstflag && Flags[ant][chan*Ntimes[ant]+offset]) {
                firstflag = true;
                firstIdx = chan-1;
            }
            if (firstflag && !Flags[ant][chan*Ntimes[ant]+offset]) {
                firstflag = false;
                lastIdx = chan;
                for (long intpIdx=firstIdx+1; intpIdx<lastIdx; intpIdx++) {
                    double frchan = ((double) (intpIdx - firstIdx)) / ((double) (lastIdx - firstIdx));
                    // std::cout << frchan<<" "<<GainAmp[0][ant][chan][lastIdx]<<" "<<GainAmp[0][ant][firstIdx][lastIdx]<<"\n";
                    GainAmp[0][ant][intpIdx][lastIdx] = GainAmp[0][ant][chan][lastIdx]*frchan + GainAmp[0][ant][firstIdx][lastIdx]*(1.0 - frchan);
                    GainAmp[1][ant][intpIdx][lastIdx] = GainAmp[1][ant][chan][lastIdx]*frchan + GainAmp[1][ant][firstIdx][lastIdx]*(1.0 - frchan);
                    GainPhase[0][ant][intpIdx][lastIdx] = GainPhase[0][ant][chan][lastIdx]*frchan + GainPhase[0][ant][firstIdx][lastIdx]*(1.0 - frchan);
                    GainPhase[1][ant][intpIdx][lastIdx] = GainPhase[1][ant][chan][lastIdx]*frchan + GainPhase[1][ant][firstIdx][lastIdx]*(1.0 - frchan);
                    Flags[ant][intpIdx*Ntimes[ant]+offset] = false;
                }
            }
        }

   }//for(times)
}

/** Prepares the instance for the frequency interpolation. User must provide
 the array of frequencies to interpolate to (freqs[0..Nchan-1]) and the number
 of channels in that array (mschan).
 TODO: misleadlingly named function
*/
void CalTable::setMapping(long mschan, double *freqs)
{
    assert(freqs != NULL);

    MSChan = mschan;
    gainChanged = true;
    currTime = -1.0;
    deltaNu = (freqs[1]-freqs[0])*1e-9;
    deltaNu0 = (freqs[0]-Freqs[0])*1e-9;

    delete[] K0;
    delete[] I0;
    delete[] I1;
    K0 = new double[mschan];
    I0 = new long[mschan];
    I1 = new long[mschan];

    for (long ant=0; ant<Nants; ant++) {
        preKt[ant] = -1.0 ;
        delete[] bufferGain[0][ant];
        delete[] bufferGain[1][ant];
        bufferGain[0][ant] = new std::complex<float>[mschan];
        bufferGain[1][ant] = new std::complex<float>[mschan];
    }

    if (SignFreq) {
        for (long i=0; i<mschan; i++) {
            if (freqs[i] <= Freqs[0]) {
                I0[i] = 0;
                I1[i] = 0;
                K0[i] = 1.0;
            } else if (freqs[i] >= Freqs[Nchan-1]) {
                I0[i] = Nchan-1;
                I1[i] = 0;
                K0[i] = 1.0;
           } else {
               for (long auxI=Nchan-1; auxI>=0; auxI--) {
                   if (freqs[i] >= Freqs[auxI]) {
                       I0[i] = auxI;
                       I1[i] = auxI+1;
                       double mmod = Freqs[auxI+1] - Freqs[auxI];
                       K0[i] = 1.0 - (freqs[i]-Freqs[auxI])/mmod;
                       break;
                   }
               }
           }
        }
    } else {
        for (long i=0; i<mschan; i++) {
            if (freqs[i] >= Freqs[0]) {
                I0[i] = 0;
                I1[i] = 0;
                K0[i] = 1.0;
            } else if (freqs[i] <= Freqs[Nchan-1]) {
                I0[i] = Nchan-1;
                I1[i] = 0;
                K0[i] = 1.0;
            } else {
                for (long auxI=Nchan-1; auxI>=0; auxI--) {
                    if (freqs[i] <= Freqs[auxI]) {
                        I0[i] = auxI;
                        I1[i] = auxI+1;
                        double mmod = Freqs[auxI+1] - Freqs[auxI];
                        K0[i] = 1.0 - (freqs[i]-Freqs[auxI])/mmod;
                        break;
                    }
                }
            }
        }
    }//if(SignFreq)

}


bool CalTable::setInterpolationTime(double itime)
{
  if (itime == currTime) {gainChanged = false; return gainChanged;};

  long i; //, auxI;
  long ti0 = 0;
  long ti1 = 0; 
  double Kt = 0.0;
  long Nts; 
  double auxD, auxD2;
  int iant;

  gainChanged = false;

// Find the right time interpolation:  

 for (iant=0; iant<Nants; iant++) {

  Nts = Ntimes[iant];

  if (Nts == 1) {
    pret0[iant] = 0;
    pret1[iant] = 0;
    if (preKt[iant]<0.0){gainChanged=true;};
    preKt[iant] = 1.0;
  } else {
  
  if (itime<=Time[iant][0]) {
    ti0 = 0; ti1 = 0; Kt = 1.0;} 
  else if (itime>=Time[iant][Nts-1]) {
    ti0 = Nts-1; ti1 = 0; Kt = 1.0;}
  else {
    for (i=Nts-1; i>=0; i--) {
      if (itime>Time[iant][i]) {
        ti1 = i+1;
        ti0 = i;
        auxD = Time[iant][i];
        auxD2 = Time[iant][i+1];
        if (isLinear){
          Kt = (1.0 - (itime - auxD)/(auxD2-auxD));
        } else {
          Kt = 1.0;
        };
        break;    
      };
    };
  };



   pret0[iant] = ti0;
   pret1[iant] = ti1;
   preKt[iant] = Kt;
   gainChanged = true;
   };

 };


   currTime = itime;
   return gainChanged;



};




/* Interpolates the gains of an antenna (iant) at a given time, itime, 
   and applies them to the arrays re[2] and im[2] (elements of these arrays
   are the different polarizations (X,Y)). User should have run "setMapping" before,
   in case that the table channel frequencies do not match with the array of frequencies
   where the interpolation is desired. 
   Then, if mode==0, the gains are just written in re and im. 
   If mode==1, the gains are ADDED to the already-existing values in re and im.
   If mode==2, the gains are MULTIPLIED to the already-existing values in re and im. 
 -------------------
*/

void CalTable::applyInterpolation(int iant, int mode, std::complex<float> *gain[2]) {

  long i;
  long ti0, ti1;
  double auxD;
  double Kt, Kt2;

  double auxF0, auxF1, auxF2, auxF3, auxT0, auxT1, auxT2, auxT3;

  ti0 = pret0[iant];
  ti1 = pret1[iant];
  Kt = preKt[iant];
  Kt2 = 1.0-Kt;

 // printf("CALLED!\n");

// Interpolate in frequency (first) and time (second):


  if (gainChanged || firstTime[iant]) {

   firstTime[iant]=false;

 // if (true) {

  for (i=0; i<MSChan; i++) {

     auxF0 = GainAmp[0][iant][I0[i]][ti0]*K0[i];
     auxF1 = GainAmp[1][iant][I0[i]][ti0]*K0[i];
     auxF2 = GainPhase[0][iant][I0[i]][ti0]*K0[i];
     auxF3 = GainPhase[1][iant][I0[i]][ti0]*K0[i];

     if (I1[i] >0) {
       auxF0 += GainAmp[0][iant][I1[i]][ti0]*(1.-K0[i]);
       auxF1 += GainAmp[1][iant][I1[i]][ti0]*(1.-K0[i]);
       auxF2 += GainPhase[0][iant][I1[i]][ti0]*(1.-K0[i]);
       auxF3 += GainPhase[1][iant][I1[i]][ti0]*(1.-K0[i]);
     };

  //  if (i==250){
  //    printf("ORIG: %3e %3e %3e %3e | %3e %3e\n",auxF0,auxF1,auxF2,auxF3, GainAmp[0][iant][I0[i]][ti0], K0[i]);
  //    holadola[0]=auxF0;holadola[1]=auxF1;holadola[2]=auxF2;holadola[3]=auxF3;
   // };

     if (ti1 > 0) {

     auxT0 = GainAmp[0][iant][I0[i]][ti1]*K0[i];
     auxT1 = GainAmp[1][iant][I0[i]][ti1]*K0[i];
     auxT2 = GainPhase[0][iant][I0[i]][ti1]*K0[i];
     auxT3 = GainPhase[1][iant][I0[i]][ti1]*K0[i];

     if (I1[i] >0) {
       auxT0 += GainAmp[0][iant][I1[i]][ti1]*(1.-K0[i]);
       auxT1 += GainAmp[1][iant][I1[i]][ti1]*(1.-K0[i]);
       auxT2 += GainPhase[0][iant][I1[i]][ti1]*(1.-K0[i]);
       auxT3 += GainPhase[1][iant][I1[i]][ti1]*(1.-K0[i]);
     };
        auxF0 = auxF0*Kt+auxT0*Kt2;
        auxF1 = auxF1*Kt+auxT1*Kt2;
        auxF2 = auxF2*Kt+auxT2*Kt2;
        auxF3 = auxF3*Kt+auxT3*Kt2;


     };


       if (isDelay){
         auxD = TWOPI*(((double)i + 0.5)*deltaNu+deltaNu0);
         bufferGain[0][iant][i] = (std::complex<float>) std::polar(1.0,auxF0*auxD);
         bufferGain[1][iant][i] = (std::complex<float>) std::polar(1.0,auxF1*auxD);
       } else if (isTsys) {
         bufferGain[0][iant][i].real(1./sqrt(auxF0));
         bufferGain[0][iant][i].imag(0.0);
         bufferGain[1][iant][i].real(1./sqrt(auxF1));
         bufferGain[1][iant][i].imag(0.0);
       } else if (isDterm) {
         bufferGain[0][iant][i].real(auxF0);
         bufferGain[0][iant][i].imag(auxF2);
         bufferGain[1][iant][i].real(auxF1);
         bufferGain[1][iant][i].imag(auxF3);
       } else {
         bufferGain[0][iant][i] = (std::complex<float>) std::polar(auxF0,auxF2);
         bufferGain[1][iant][i] = (std::complex<float>) std::polar(auxF1,auxF3);
       };

        
    //     if(i==0){printf("\nPH INT: %.2f %.2f",std::abs(bufferGain[0][iant][i]),std::arg(bufferGain[0][iant][i]));};


  };

  };  // Comes from if(gainChanged)



  switch (mode) {
    case 0: // Normal mode.

      for (i=0; i< MSChan; i++) {
          gain[0][i] = bufferGain[0][iant][i];
          gain[1][i] = bufferGain[1][iant][i];
      };

        break;

    case 1:   // Addition mode. 
      for (i=0; i< MSChan; i++) {
         gain[0][i] +=  bufferGain[0][iant][i];
         gain[1][i] +=  bufferGain[1][iant][i];
      };

       break;

    case 2:  // Product mode.
      for (i=0; i< MSChan; i++) {
         gain[0][i] *=  bufferGain[0][iant][i]; 
         gain[1][i] *=  bufferGain[1][iant][i];
     };

       break;

  };

 // if (gainChanged) {
 //   printf("MODE: %i  %i  %i | %3e  %3e | %3e  %3e\n",mode,Nchan, iant, std::abs(gain[0][250]), std::abs(bufferGain[0][iant][250]), std::arg(gain[0][250]), std::arg(bufferGain[0][iant][250]));
 //  printf("MODE: %i  %i  %i | %3e  %3e  %3e %3e\n",mode,I0[250],I1[250], K0[250], GainAmp[0][iant][I0[i]][ti0], GainAmp[1][iant][I0[i]][ti0], Kt );
 //  printf("%i %.5e  %.5e  %.5e  %.5e |  %.5e  %.5e  %.5e %.5e  \n",iant, gain[0][250].real(), gain[0][250].imag(), gain[1][250].real(), gain[1][250].imag(), bufferGain[0][iant][250].real(), bufferGain[0][iant][250].imag(), bufferGain[1][iant][250].real(), bufferGain[1][iant][250].imag());
//  };


};




bool CalTable::getInterpolation(int iant, int ichan, std::complex<float> gain[2]){

  success = true;
  gain[0] = 1.0 ; gain[1] = 1.0 ;
  if(iant>=Nants || ichan > Nchan){success = false;} else {
    gain[0] = bufferGain[0][iant][ichan];
    gain[1] = bufferGain[1][iant][ichan];
    success = true;
  };

//  int i;
//  for (i=0;i<Nants;i++){
//  printf("\nIa: %i  PH1: %.2f  PH2: %.2f",i,std::abs(bufferGain[0][i][ichan]),std::arg(bufferGain[0][i][ichan]));
//  };
//  if(!success){printf("\n %i %i %i %i\n",iant,Nants,ichan,Nchan);};

  return success;

};

