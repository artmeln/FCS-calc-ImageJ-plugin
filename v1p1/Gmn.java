/*
 * Copyright (C) 2016-2017 Artem Melnykov
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.io.FileInputStream;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Gmn {

    // parameters for correlation function calculation
    private int t0=1;     // lowest time base value
    private int nc=16;     // number of cascades
    private int np=16;     // number of points in a cascade

    // data arrays and array sizes for both channels
    private byte[] bufferInA;
    private byte[] bufferInB;
    private int bufferInSizeA;
    private int bufferInSizeB;

    // photon arrays used in correlation
    private int[] photonsIat;   // interarrival times
    private int[] wAint;        // integer weights (channel A)
    private int[] wBint;        // integer weights (channel B)
    private int iatA;           // current inter-arrival time (ch. A)
    private int iatB;           // current inter-arrival time (ch. B)
    private int nEvents;        // current number of events in photonsIat
    
    // calculated correlation function values
    private double[] g11; // regular auto-correlation
    private double[] g12; // higher oder auto-correlation
    private double[] g21; // higher oder auto-correlation
    private double[] g13; // higher oder auto-correlation
    private double[] g31; // higher oder auto-correlation
    private double[] g22; // higher oder auto-correlation
    
    private double[] gAA; // auto and cross-correlations
    private double[] gAB; // auto and cross-correlations
    private double[] gBA; // auto and cross-correlations
    private double[] gBB; // auto and cross-correlations

    // internal parameters used in data processing
    private int evA;            // trackers for event numbers in
    private int evB;            // data blocks read from files
    private int ccA;            // trackers for event numbers in 
    private int ccB;            // data blocks obtained from custom hardware
    private int[] ciatarray;    // array of cumulative inter-arrival times
    private int[] wAintArray;   // corresponding array of weights (channel A) 
    private int[] wBintArray;   // corresponding array of weights (channel B)
    private int[] mta;          // current mta values for each cascade
    private int[] iatnew;       // current inter-arrival time for each cascade
    private int[] wAintTemp;    // current bin weights for each cascade (channel A)
    private int[] wBintTemp;    // current bin weights for each cascade (channel B)
    private byte[] sp;           // stop index for each cascade
    
    private double[] tlast;     // last mta value for each cascade
    private double[] wAtotal;   // total number of photons for each cascade (channel A)
    private double[] wAtotal2;   // total number of photons for each cascade (channel A) ^2
    private double[] wAtotal3;   // total number of photons for each cascade (channel A) ^3
    private double[] wBtotal;   // total number of photons for each cascade (channel B)
    private double[] s01a;      // right normalization constant (ch. A)
    private double[] s01b;      // right normalization constant (ch. B)
    private double[] s11aa;     // correlation function values (not normalized)
    private double[] s11bb;     // correlation function values (not normalized)
    private double[] s11ab;     // correlation function values (not normalized)
    private double[] s11ba;     // correlation function values (not normalized)

    private double[] s02a;      // higher order moments (ch. A)
    private double[] s03a;      // higher order moments (ch. A)
    private double[] s12aa;     // correlation function values (not normalized)
    private double[] s21aa;     // correlation function values (not normalized)
    private double[] s13aa;     // correlation function values (not normalized)
    private double[] s31aa;     // correlation function values (not normalized)
    private double[] s22aa;     // correlation function values (not normalized)

    private double countrateA;  // averagecountrate for all data (ch. A)
    private double countrateB;  // average countrate for all data (ch. B)
    
    // various time values
    private double dt;              // time resolution of the data (in seconds)
    private double bint;            // binning time for intensity traces
    
    // intensity trace values
    private double[] intensityTraceX;        // the entire intensity trace (time value)
    private double[] intensityTraceYA;        // the entire intensity trace (intensity value)
    private double[] intensityTraceYB;        // the entire intensity trace (intensity value)
    private int intensityTracePosition;  // position of the last intensity trace entry
    private int binSize;    // determines number of bins in photonsIat
    private double macroTime;   // macro time value from the start of experiment
    
    public int initializeGmn(String correlationType) {
        int ff=0;
        try {
            if (correlationType=="auto") {
                evA = 32769;
                evB = 32769;
                ccA = 0;
                bufferInA = new byte[32768];
                bufferInSizeA = 32768;                
                iatA=0;
                photonsIat = new int[32768];
                wAint = new int[32768];
                wBint = null;
                nEvents=0;
                ciatarray = new int[256*nc];
                for (int kk=0; kk<256*nc; kk++) {ciatarray[kk]=10000;}
                wAintArray = new int[256*nc];
                mta = new int[nc];
                iatnew = new int[nc];
                wAintTemp = new int[nc];
                sp = new byte[nc];
                for (byte k:sp) {k=-128;}

                g11 = new double[np*(nc+1)];
                for (int kk=0; kk<np*(nc+1); kk++) {g11[kk]=0;}
                tlast = new double[nc];
                wAtotal = new double[nc];
                s01a = new double[np*(nc+1)];
                s11aa = new double[np*(nc+1)];

                intensityTraceX = new double[32768];
                intensityTraceYA = new double[32768];
                intensityTracePosition = 0;
                binSize = 4096; // so that there are 32768/4096 = 8 bins
                macroTime = 0;

                return ff=1;
            } else if (correlationType=="cross") {
                evA = 32769;
                evB = 32769;
                bufferInA = new byte[32768];
                bufferInSizeA = 32768;                
                bufferInB = new byte[32768];
                bufferInSizeB = 32768;                
                iatA=0;
                iatB=0;
                photonsIat = new int[32768];
                wAint = new int[32768];
                wBint = new int[32768];
                nEvents=0;
                ciatarray = new int[256*nc];
                for (int kk=0; kk<256*nc; kk++) {ciatarray[kk]=10000;}
                wAintArray = new int[256*nc];
                wBintArray = new int[256*nc];
                mta = new int[nc];
                iatnew = new int[nc];
                wAintTemp = new int[nc];
                wBintTemp = new int[nc];
                sp = new byte[nc];
                for (byte k:sp) {k=-128;}

                gAA = new double[np*(nc+1)];
                gAB = new double[np*(nc+1)];
                gBA = new double[np*(nc+1)];
                gBB = new double[np*(nc+1)];
                tlast = new double[nc];
                wAtotal = new double[nc];
                wBtotal = new double[nc];
                s01a = new double[np*(nc+1)];
                s01b = new double[np*(nc+1)];
                s11aa = new double[np*(nc+1)];
                s11ab = new double[np*(nc+1)];
                s11ba = new double[np*(nc+1)];
                s11bb = new double[np*(nc+1)];

                intensityTraceX = new double[32768];
                intensityTraceYA = new double[32768];
                intensityTraceYB = new double[32768];
                intensityTracePosition = 0;
                binSize = 4096; // so that there are 32768/4096 = 8 bins
                macroTime = 0;

                return ff=1;
            } else if (correlationType=="autoHOmlt") {
                evA = 32769;
                evB = 32769;
                ccA = 0;
                bufferInA = new byte[32768];
                bufferInSizeA = 32768;                
                iatA=0;
                photonsIat = new int[32768];
                wAint = new int[32768];
                wBint = null;
                nEvents=0;
                ciatarray = new int[256*nc];
                for (int kk=0; kk<256*nc; kk++) {ciatarray[kk]=10000;}
                wAintArray = new int[256*nc];
                mta = new int[nc];
                iatnew = new int[nc];
                wAintTemp = new int[nc];
                sp = new byte[nc];
                for (byte k:sp) {k=-128;}

                g11 = new double[np*(nc+1)];
                for (int kk=0; kk<np*(nc+1); kk++) {g11[kk]=0;}
                g12 = new double[np*(nc+1)];
                for (int kk=0; kk<np*(nc+1); kk++) {g12[kk]=0;}
                g21 = new double[np*(nc+1)];
                for (int kk=0; kk<np*(nc+1); kk++) {g21[kk]=0;}
                g13 = new double[np*(nc+1)];
                for (int kk=0; kk<np*(nc+1); kk++) {g13[kk]=0;}
                g31 = new double[np*(nc+1)];
                for (int kk=0; kk<np*(nc+1); kk++) {g31[kk]=0;}
                g22 = new double[np*(nc+1)];
                for (int kk=0; kk<np*(nc+1); kk++) {g22[kk]=0;}
                tlast = new double[nc];
                wAtotal = new double[nc];
                wAtotal2 = new double[nc];
                wAtotal3 = new double[nc];
                s01a = new double[np*(nc+1)];
                s02a = new double[np*(nc+1)];
                s03a = new double[np*(nc+1)];
                s11aa = new double[np*(nc+1)];
                s12aa = new double[np*(nc+1)];
                s21aa = new double[np*(nc+1)];
                s13aa = new double[np*(nc+1)];
                s31aa = new double[np*(nc+1)];
                s22aa = new double[np*(nc+1)];
                
                intensityTraceX = new double[32768];
                intensityTraceYA = new double[32768];
                intensityTracePosition = 0;
                binSize = 4096; // so that there are 32768/4096 = 8 bins
                macroTime = 0;

                return ff=1;
            } else {
                ff=1;
                return -1;
            }
            
        } finally { // this would happen if there is not enough memory
            if (ff==0) {
                System.out.println("Unable to initialize");
                return -2;
            }
        }
            
    }
    
        public double getMin(double[] v) {
        double vMin=1e100;
        for (double el:v) {
            if (el<vMin) {vMin=el;}
        }
        return vMin;
    }

    public double getMax(double[] v) {
        double vMax=-1e100;
        for (double el:v) {
            if (el>vMax) {vMax=el;}
        }
        return vMax;
    }

    public int byte4ToInt(byte[] byte4) {

        int intout=0;

        int temp=0;
        if (byte4[0]>=0) {temp = byte4[0];} else {temp = 256+byte4[0];}
        intout += temp;
        if (byte4[1]>=0) {temp = byte4[1];} else {temp = 256+byte4[1];}
        intout += 256*temp;
        if (byte4[2]>=0) {temp = byte4[2];} else {temp = 256+byte4[2];}
        intout += 256*256*temp;
        if (byte4[3]>=0) {temp = byte4[3];} else {temp = 256+byte4[3];}
        intout += 256*256*256*temp; 

        return intout;
    }

    public void setParams(int baseTime, int nCascades, int nPointsPerCascade) {
        t0 = baseTime;
        nc = nCascades;
        np = nPointsPerCascade;
    }

    public int[] getParams() {
        int[] a;
        a = new int[3];
        a[0] = t0;
        a[1] = nc;
        a[2] = np;
        return a;
    }
    
    public int getNofEvents() {
        return nEvents;
    }
    
    public int[] getPhotonsIat() {
        return photonsIat;
    }
    
    public int[] getWeightsA() {
        return wAint;
    }

    public int[] getWeightsB() {
        return wBint;
    }
    
    public double[] getWaTotal() {
        return wAtotal;
    }

    public double[] getG11() {
        return g11;
    }
    
    public double[] getG12() {
        return g12;
    }
    public double[] getG21() {
        return g21;
    }
    public double[] getG13() {
        return g13;
    }
    public double[] getG31() {
        return g31;
    }
    public double[] getG22() {
        return g22;
    }

    public double[] getGaa() {
        return gAA;
    }
    
    public double[] getGab() {
        return gAB;
    }
    
    public double[] getGba() {
        return gBA;
    }
    
    public double[] getGbb() {
        return gBB;
    }
    
    public int getIntensityTracePosition() {
        return intensityTracePosition;
    }
    
    public double[] getIntensityTraceX() {
        double[] truncatedTrace;
        truncatedTrace = new double[intensityTracePosition];
        int ii;
        for (ii=0; ii<intensityTracePosition; ii++) {
            truncatedTrace[ii] = intensityTraceX[ii];
        }
        return truncatedTrace;
    }
    
    public double[] getIntensityTraceYA() {
        double[] truncatedTrace;
        truncatedTrace = new double[intensityTracePosition];
        int ii;
        for (ii=0; ii<intensityTracePosition; ii++) {
            truncatedTrace[ii] = intensityTraceYA[ii];
        }
        return truncatedTrace;
    }
    
    public double[] getIntensityTraceYB() {
        double[] truncatedTrace;
        truncatedTrace = new double[intensityTracePosition];
        int ii;
        for (ii=0; ii<intensityTracePosition; ii++) {
            truncatedTrace[ii] = intensityTraceYB[ii];
        }
        return truncatedTrace;
    }
    
    public int updateOnlineCustomHardwareOneCh(byte[] inputData, int nDataPoints) {

        // a chunk of a large data block is converted into inter-arrival time array
        // all events are assigned a weight of 1
        nEvents=0;
        while ( nEvents<32768 && ccA<nDataPoints ) {
            photonsIat[nEvents]=inputData[ccA];
            wAint[nEvents] = 1;
            nEvents++;
            ccA++;
        }
        if (ccA==nDataPoints) {
            ccA=0;
            return nDataPoints;
        } else {
            return ccA;
        }
    }

    public int updateDataFlexOneCh(FileInputStream inputFile) {

        nEvents=0;
        try {
            // Flex8 data format uses unsigned byte data type
            // 0 is 0 in both signed and unsigned formats
            // unsigned and signed coincide up to 127
            // after that the signed number rolls over to -128
            // eventually unsigned 255 corresponds to signed -1
            if (evA>=bufferInSizeA) { // read in the next chunk of data
                bufferInSizeA = inputFile.read(bufferInA);
                evA=0;
            }
            while ( nEvents<32768 && bufferInSizeA!=-1 ) {
                if (bufferInA[evA]==-1) {
                    iatA += 255;
                } else {
                    if (bufferInA[evA]>=0) {
                        iatA += bufferInA[evA]+1;
                    } else {
                        iatA += 256+bufferInA[evA]+1;
                    }
                    photonsIat[nEvents]=iatA;
                    wAint[nEvents] = 1;
                    iatA=0;
                    nEvents++;
                }
                evA++;
                if (evA>=bufferInSizeA) { // read in the next chunk of data
                    bufferInSizeA = inputFile.read(bufferInA);
                    evA=0;
                }
            }
            return bufferInSizeA;
        } finally {
            return bufferInSizeA;
        }
    }
    
    public int updateDataFlexTwoCh(FileInputStream inputFileA, FileInputStream inputFileB) {
        
        nEvents=0;
        try {
            if (evA>=bufferInSizeA) { // read in the next chunk of data (ch. A)
                bufferInSizeA = inputFileA.read(bufferInA);
                evA=0;
            }
            if (evB>=bufferInSizeB) { // read in the next chunk of data (ch. B)
                bufferInSizeB = inputFileB.read(bufferInB);
                evB=0;
            }
                        
            while ( nEvents<32768 && bufferInSizeA!=-1 && bufferInSizeB!=-1 ) {
                // Channel A
                if (iatA==0) {
                    // Channel A
                    while (bufferInA[evA]==-1 && bufferInSizeA!=-1) {  // sum up overrun clock cycles
                        iatA += 255;
                        evA++;
                        if (evA>=bufferInSizeA) { // read in the next chunk of data
                            bufferInSizeA = inputFileA.read(bufferInA);
                            evA=0;
                        }
                    }
                    if ( bufferInSizeA!=-1 ) { // add the photon-interrupted clock cycle
                        if (bufferInA[evA]>=0)  { iatA += bufferInA[evA]+1; }
                        else { iatA += 256+bufferInA[evA]+1; }
                        if (evA>=bufferInSizeA) { // read in the next chunk of data
                            bufferInSizeA = inputFileA.read(bufferInA);
                            evA=0;
                        }
                    }
                }
                if (iatB==0) {
                    // Channel B
                    while (bufferInB[evB]==-1 && bufferInSizeB!=-1) {  // sum up overrun clock cycles
                        iatB += 255;
                        evB++;
                        if (evB>=bufferInSizeB) { // read in the next chunk of data
                            bufferInSizeB = inputFileB.read(bufferInB);
                            evB=0;
                        }
                    }
                    if ( bufferInSizeB!=-1 ) { // add the photon-interrupted clock cycle
                        if (bufferInB[evB]>=0)  { iatB += bufferInB[evB]+1; }
                        else { iatB += 256+bufferInB[evB]+1; }
                        if (evB>=bufferInSizeB) { // read in the next chunk of data
                            bufferInSizeB = inputFileB.read(bufferInB);
                            evB=0;
                        }
                    }
                }
                
                // decide which iat is smaller (A or B) and record this event
                if ( bufferInSizeA!=-1 && bufferInSizeB!=-1 ) {
                    if (iatA>iatB) {
                        photonsIat[nEvents]=iatB;
                        wBint[nEvents] = 1;
                        wAint[nEvents] = 0;
                        iatA=iatA-iatB;
                        iatB=0;
                        evB++;
                    } else if (iatA<iatB) {
                        photonsIat[nEvents]=iatA;
                        wBint[nEvents] = 0;
                        wAint[nEvents] = 1;
                        iatB=iatB-iatA;
                        iatA=0;
                        evA++;
                    } else {
                        photonsIat[nEvents]=iatA;
                        wBint[nEvents] = 1;
                        wAint[nEvents] = 1;
                        iatA=0;
                        iatB=0;
                        evA++;
                        evB++;
                    }
                    nEvents++;
                }
            }
            // finish the calculation on the remaining data
            if (bufferInSizeA==-1 && bufferInSizeB!=-1 ) { // channel B only
                // carry out one channel conversion
                if (evB>=bufferInSizeB) { // read in the next chunk of data
                    bufferInSizeB = inputFileB.read(bufferInB);
                    evB=0;
                }
                while ( nEvents<32768 && bufferInSizeB!=-1 ) {
                    if (bufferInB[evB]==-1) {
                        iatB += 255;
                    } else {
                        if (bufferInB[evB]>=0) {
                            iatB += bufferInB[evB]+1;
                        } else {
                            iatB += 256+bufferInB[evB]+1;
                        }
                        photonsIat[nEvents]=iatB;
                        wAint[nEvents] = 0;
                        wBint[nEvents] = 1;
                        iatB=0;
                        nEvents++;
                    }
                    evB++;
                    if (evB>=bufferInSizeB) { // read in the next chunk of data
                        bufferInSizeB = inputFileB.read(bufferInB);
                        evB=0;
                    }
                }
            }
            if (bufferInSizeA!=-1 && bufferInSizeB==-1 ) { // channel A only
                // carry out one channel conversion
                if (evA>=bufferInSizeA) { // read in the next chunk of data
                    bufferInSizeA = inputFileA.read(bufferInA);
                    evA=0;
                }
                while ( nEvents<32768 && bufferInSizeA!=-1 ) {
                    if (bufferInA[evA]==-1) {
                        iatA += 255;
                    } else {
                        if (bufferInA[evA]>=0) {
                            iatA += bufferInA[evA]+1;
                        } else {
                            iatA += 256+bufferInA[evA]+1;
                        }
                        photonsIat[nEvents]=iatA;
                        wAint[nEvents] = 1;
                        wBint[nEvents] = 0;
                        iatA=0;
                        nEvents++;
                    }
                    evA++;
                    if (evA>=bufferInSizeA) { // read in the next chunk of data
                        bufferInSizeA = inputFileA.read(bufferInA);
                        evA=0;
                    }
                }
            }            
            return bufferInSizeA;
        } finally {
            return bufferInSizeA;
        }
    }
    
    public int readHeaderConfocor2(FileInputStream inputFile) {
        int aa;
        byte[] bb;
        bb = new byte[30];
        
        try {
            aa = inputFile.read(bb);
            return 1;
        } finally {
            return 1;
        }
        
    }

    public int updateDataConfocor2chA(FileInputStream inputFile) {
        
        nEvents=0;
        try {
            if (evA>=bufferInSizeA) { // read in the next chunk of data
                bufferInSizeA = inputFile.read(bufferInA);
                evA=0;
            }
             
            while ( nEvents<32768 && bufferInSizeA!=-1 ) {
                if (bufferInA[evA]!=0) {
                    if (bufferInA[evA]==-1) {
                        iatA += 255;
                        if (bufferInA[evA+1]==0) {           // 00000000
                            iatA += 3;
                        } else {
                            if ( (bufferInA[evA+1] & 1)!=0 ) {    // 00000001
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                nEvents++;
                            }
                            if ( (bufferInA[evA+1] & (1 << 1))!=0 ) {    // 00000010
                                //photonsIat[nEvents]=iat;
                                //iat=0;
                                //wAint[nEvents] = 1;                   
                                //nEvents++;
                            } 
                            iatA += 1;
                            if ( (bufferInA[evA+1] & (1 << 2))!=0 ) {    // 00000100
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                nEvents++;
                            }
                            if ( (bufferInA[evA+1] & (1 << 3))!=0 ) {    // 00001000
                                //photonsIat[nEvents]=iat;
                                //iat=0;
                                //wAint[nEvents] = 1;                   
                                //nEvents++;
                            }
                            iatA += 1;
                            if ( (bufferInA[evA+1] & (1 << 4))!=0 ) {    // 00010000
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                nEvents++;
                            }   
                            if ( (bufferInA[evA+1] & (1 << 5))!=0 ) {    // 00100000
                                //photonsIat[nEvents]=iat;
                                //iat=0;
                                //wAint[nEvents] = 1;                   
                                //nEvents++;
                            }
                            iatA += 1;
                            if ( (bufferInA[evA+1] & (1 << 6))!=0 ) {    // 01000000
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                nEvents++;
                            }
                            if ( (bufferInA[evA+1] & (1 << 7))!=0 ) {    // 10000000
                                //photonsIat[nEvents]=iat;
                                //iat=0;
                                //wAint[nEvents] = 1;                   
                                //nEvents++;
                            }
                        }
                    } else {
                        if (bufferInA[evA]>=0) {
                            iatA += bufferInA[evA];
                        } else {
                            iatA += 256+bufferInA[evA];
                        }
                        if ( (bufferInA[evA+1] & 1)!=0 ) {    // 00000001
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            nEvents++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 1))!=0 ) {    // 00000010
                            //photonsIat[nEvents]=iat;
                            //iat=0;
                            //wAint[nEvents] = 1;                   
                            //nEvents++;
                        }
                        iatA += 1;
                        if ( (bufferInA[evA+1] & (1 << 2))!=0 ) {    // 00000100
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            nEvents++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 3))!=0 ) {    // 00001000
                            //photonsIat[nEvents]=iat;
                            //iat=0;
                            //wAint[nEvents] = 1;                   
                            //nEvents++;
                        }
                        iatA += 1;
                        if ( (bufferInA[evA+1] & (1 << 4))!=0 ) {    // 00010000
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            nEvents++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 5))!=0 ) {    // 00100000
                            //photonsIat[nEvents]=iat;
                            //iat=0;
                            //wAint[nEvents] = 1;                   
                            //nEvents++;
                        }
                        iatA += 1;
                        if ( (bufferInA[evA+1] & (1 << 6))!=0 ) {    // 01000000
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            nEvents++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 7))!=0 ) {    // 10000000
                            //photonsIat[nEvents]=iat;
                            //iat=0;
                            //wAint[nEvents] = 1;                   
                            //nEvents++;
                        }
                    }
                    
                    evA += 2;
                    if (evA>=bufferInSizeA) { // read in the next chunk of data
                        bufferInSizeA = inputFile.read(bufferInA);
                        evA=0;
                    }
                
                } else { bufferInSizeA=-1; } // end of data has been reached
            } 
            return bufferInSizeA;
        } finally {
            return bufferInSizeA;
        }
    }
    
    public int updateDataConfocor2chB(FileInputStream inputFile) {
       
        nEvents=0;
        try {
            if (evA>=bufferInSizeA) { // read in the next chunk of data
                bufferInSizeA = inputFile.read(bufferInA);
                evA=0;
            }
             
            while ( nEvents<32768 && bufferInSizeA!=-1 ) {
                if (bufferInA[evA]!=0) {
                    if (bufferInA[evA]==-1) {
                        iatA += 255;
                        if (bufferInA[evA+1]==0) {           // 00000000
                            iatA += 3;
                        } else {
                            if ( (bufferInA[evA+1] & 1)!=0 ) {    // 00000001
                                //photonsIat[nPhotonsA]=iat;
                                //iat=0;
                                //wAint[nPhotonsA] = 0;                   
                                //nPhotonsA++;
                            }
                            if ( (bufferInA[evA+1] & (1 << 1))!=0 ) {    // 00000010
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                nEvents++;
                            } 
                            iatA += 1;
                            if ( (bufferInA[evA+1] & (1 << 2))!=0 ) {    // 00000100
                                //photonsIat[nPhotonsA]=iat;
                                //iat=0;
                                //wAint[nPhotonsA] = 0;                   
                                //nPhotonsA++;
                            }
                            if ( (bufferInA[evA+1] & (1 << 3))!=0 ) {    // 00001000
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                nEvents++;
                            }
                            iatA += 1;
                            if ( (bufferInA[evA+1] & (1 << 4))!=0 ) {    // 00010000
                                //photonsIat[nPhotonsA]=iat;
                                //iat=0;
                                //wAint[nPhotonsA] = 0;                   
                                //nPhotonsA++;
                            }   
                            if ( (bufferInA[evA+1] & (1 << 5))!=0 ) {    // 00100000
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                nEvents++;
                            }
                            iatA += 1;
                            if ( (bufferInA[evA+1] & (1 << 6))!=0 ) {    // 01000000
                                //photonsIat[nPhotonsA]=iat;
                                //iat=0;
                                //wAint[nPhotonsA] = 0;                   
                                //nPhotonsA++;
                            }
                            if ( (bufferInA[evA+1] & (1 << 7))!=0 ) {    // 10000000
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                nEvents++;
                            }
                        }
                    } else {
                        if (bufferInA[evA]>=0) {
                            iatA += bufferInA[evA];
                        } else {
                            iatA += 256+bufferInA[evA];
                        }
                        if ( (bufferInA[evA+1] & 1)!=0 ) {    // 00000001
                            //photonsIat[nPhotonsA]=iat;
                            //iat=0;
                            //wAint[nPhotonsA] = 0;                   
                            //nPhotonsA++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 1))!=0 ) {    // 00000010
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            nEvents++;
                        }
                        iatA += 1;
                        if ( (bufferInA[evA+1] & (1 << 2))!=0 ) {    // 00000100
                            //photonsIat[nPhotonsA]=iat;
                            //iat=0;
                            //wAint[nPhotonsA] = 0;                   
                            //nPhotonsA++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 3))!=0 ) {    // 00001000
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            nEvents++;
                        }
                        iatA += 1;
                        if ( (bufferInA[evA+1] & (1 << 4))!=0 ) {    // 00010000
                            //photonsIat[nPhotonsA]=iat;
                            //iat=0;
                            //wAint[nPhotonsA] = 0;                   
                            //nPhotonsA++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 5))!=0 ) {    // 00100000
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            nEvents++;
                        }
                        iatA += 1;
                        if ( (bufferInA[evA+1] & (1 << 6))!=0 ) {    // 01000000
                            //photonsIat[nPhotonsA]=iat;
                            //iat=0;
                            //wAint[nPhotonsA] = 0;                   
                            //nPhotonsA++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 7))!=0 ) {    // 10000000
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            nEvents++;
                        }
                    }
                    
                    evA += 2;
                    if (evA>=bufferInSizeA) { // read in the next chunk of data
                        bufferInSizeA = inputFile.read(bufferInA);
                        evA=0;
                    }
                
                } else { bufferInSizeA=-1; } // end of data has been reached
            } 
            return bufferInSizeA;
        } finally {
            return bufferInSizeA;
        }
    }
    
    public int updateDataConfocor2chAB(FileInputStream inputFile) {
        
        nEvents=0;
        try {
            if (evA>=bufferInSizeA) { // read in the next chunk of data
                bufferInSizeA = inputFile.read(bufferInA);
                evA=0;
            }
             
            while ( nEvents<32768 && bufferInSizeA!=-1 ) {
                if (bufferInA[evA]!=0) {
                    if (bufferInA[evA]==-1) {
                        iatA += 255;
                        if (bufferInA[evA+1]==0) {           // 00000000
                            iatA += 3;
                        } else {
                            if ( (bufferInA[evA+1] & 1)!=0 ) {    // 00000001
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                wBint[nEvents] = 0;                   
                                nEvents++;
                            }
                            if ( (bufferInA[evA+1] & (1 << 1))!=0 ) {    // 00000010
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 0;                   
                                wBint[nEvents] = 1;                   
                                nEvents++;
                            } 
                            iatA += 1;
                            if ( (bufferInA[evA+1] & (1 << 2))!=0 ) {    // 00000100
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                wBint[nEvents] = 0;                   
                                nEvents++;
                            }
                            if ( (bufferInA[evA+1] & (1 << 3))!=0 ) {    // 00001000
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 0;                   
                                wBint[nEvents] = 1;                   
                                nEvents++;
                            }
                            iatA += 1;
                            if ( (bufferInA[evA+1] & (1 << 4))!=0 ) {    // 00010000
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                wBint[nEvents] = 0;                   
                                nEvents++;
                            }   
                            if ( (bufferInA[evA+1] & (1 << 5))!=0 ) {    // 00100000
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 0;                   
                                wBint[nEvents] = 1;                   
                                nEvents++;
                            }
                            iatA += 1;
                            if ( (bufferInA[evA+1] & (1 << 6))!=0 ) {    // 01000000
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 1;                   
                                wBint[nEvents] = 0;                   
                                nEvents++;
                            }
                            if ( (bufferInA[evA+1] & (1 << 7))!=0 ) {    // 10000000
                                photonsIat[nEvents]=iatA;
                                iatA=0;
                                wAint[nEvents] = 0;                   
                                wBint[nEvents] = 1;                   
                                nEvents++;
                            }
                        }
                    } else {
                        if (bufferInA[evA]>=0) {
                            iatA += bufferInA[evA];
                        } else {
                            iatA += 256+bufferInA[evA];
                        }
                        if ( (bufferInA[evA+1] & 1)!=0 ) {    // 00000001
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            wBint[nEvents] = 0;                   
                            nEvents++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 1))!=0 ) {    // 00000010
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 0;                   
                            wBint[nEvents] = 1;                   
                            nEvents++;
                        }
                        iatA += 1;
                        if ( (bufferInA[evA+1] & (1 << 2))!=0 ) {    // 00000100
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            wBint[nEvents] = 0;                   
                            nEvents++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 3))!=0 ) {    // 00001000
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 0;                   
                            wBint[nEvents] = 1;                   
                            nEvents++;
                        }
                        iatA += 1;
                        if ( (bufferInA[evA+1] & (1 << 4))!=0 ) {    // 00010000
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            wBint[nEvents] = 0;                   
                            nEvents++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 5))!=0 ) {    // 00100000
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 0;                   
                            wBint[nEvents] = 1;                   
                            nEvents++;
                        }
                        iatA += 1;
                        if ( (bufferInA[evA+1] & (1 << 6))!=0 ) {    // 01000000
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 1;                   
                            wBint[nEvents] = 0;                   
                            nEvents++;
                        }
                        if ( (bufferInA[evA+1] & (1 << 7))!=0 ) {    // 10000000
                            photonsIat[nEvents]=iatA;
                            iatA=0;
                            wAint[nEvents] = 0;                   
                            wBint[nEvents] = 1;                   
                            nEvents++;
                        }
                    }
                    
                    evA += 2;
                    if (evA>=bufferInSizeA) { // read in the next chunk of data
                        bufferInSizeA = inputFile.read(bufferInA);
                        evA=0;
                    }
                
                } else { bufferInSizeA=-1; } // end of data has been reached
            } 
            return bufferInSizeA;
        } finally {
            return bufferInSizeA;
        }
    }

    public int readHeaderConfocor3(FileInputStream inputFile) {
        int aa;
        byte[] bb;
        bb = new byte[128];
        
        try {
            aa = inputFile.read(bb);
            return 1;
        } finally {
            return 1;
        }
        
    }

    public int updateDataConfocor3OneCh(FileInputStream inputFile) {
        
        nEvents=0;
        try {
            if (evA>=bufferInSizeA) { // read in the next chunk of data
                bufferInSizeA = inputFile.read(bufferInA);
                evA=0;
            }
             
            while ( nEvents<32768 && bufferInSizeA!=-1 ) {
                if (bufferInA[evA]>=0) {
                    iatA += bufferInA[evA];
                } else {
                    iatA += 256+bufferInA[evA];
                }
                photonsIat[nEvents]=iatA;
                iatA=0;
                evA++;
                if (bufferInA[evA]>=0) {
                    iatA += bufferInA[evA];
                } else {
                    iatA += 256+bufferInA[evA];
                }
                photonsIat[nEvents]=photonsIat[nEvents]+256*iatA;
                iatA=0;
                evA++;
                if (bufferInA[evA]>=0) {
                    iatA += bufferInA[evA];
                } else {
                    iatA += 256+bufferInA[evA];
                }
                photonsIat[nEvents]=photonsIat[nEvents]+256*256*iatA;
                iatA=0;
                evA++;
                if (bufferInA[evA]>=0) {
                    iatA += bufferInA[evA];
                } else {
                    iatA += 256+bufferInA[evA];
                }
                photonsIat[nEvents]=photonsIat[nEvents]+256*256*256*iatA;
                iatA=0;
                evA++;
                
                wAint[nEvents] = 1;
                nEvents++;

                if (evA>=bufferInSizeA) { // read in the next chunk of data
                    bufferInSizeA = inputFile.read(bufferInA);
                    evA=0;
                }
            }
            return bufferInSizeA;
        } finally {
            return bufferInSizeA;
        }
    }

    public int readHeaderLongPicoHarpPT3(FileInputStream inputFile) {
                
        int aa;
        byte[] bb2, bb4, bb6, bb12, bb16, bb18, bb20, bb60, bb256;
        bb2 = new byte[2];
        bb4 = new byte[4];
        bb6 = new byte[6];
        bb12 = new byte[12];
        bb16 = new byte[16];
        bb18 = new byte[18];
        bb20 = new byte[20];
        bb60 = new byte[60];
        bb256 = new byte[256];
        String Ident;
        String FormatVersion;
        String CreatorName;
        String CreatorVersion;
        String FileTime;
        String CRLF;
        String Comment;
        String ScriptName;   
        String BrdHardwareIdent;
        
        int NumberOfCurves, BitsPerRecord, RoutingChannels, NumberOfBoards;
        int ActiveCurve, MeasurementMode, SubMode, RangeNo, Offset;
        int AcquisitionTime, StopAt, StopOnOvfl, Restart;
        int DisplayLinLog, DisplayTimeAxisFrom, DisplayTimeAxisTo;
        int DisplayCountAxisFrom, DisplayCountAxisTo;
        
        int DisplayCurve1MapTo, DisplayCurve1Show;
        int DisplayCurve2MapTo, DisplayCurve2Show;
        int DisplayCurve3MapTo, DisplayCurve3Show;
        int DisplayCurve4MapTo, DisplayCurve4Show;
        int DisplayCurve5MapTo, DisplayCurve5Show;
        int DisplayCurve6MapTo, DisplayCurve6Show;
        int DisplayCurve7MapTo, DisplayCurve7Show;
        int DisplayCurve8MapTo, DisplayCurve8Show;
        
        float Param1Start, Param1Step, Param1End;
        float Param2Start, Param2Step, Param2End;
        float Param3Start, Param3Step, Param3End;
        
        int RepeatMode, RepeatsPerCurve, RepeatTime, RepeatWaitTime;

        int ExtDevices, Reserved1, Reserved2;
        int InpRate0, InpRate1, StopAfter, StopReason;
        int NumRecords, ImgHdrSize;
        
        try {
            aa = inputFile.read(bb16);
            Ident = new String(bb16, "UTF8");
            aa = inputFile.read(bb6);
            FormatVersion = new String(bb6, "UTF8");
            aa = inputFile.read(bb18);
            CreatorName = new String(bb18, "UTF8");
            aa = inputFile.read(bb12);
            CreatorVersion = new String(bb12, "UTF8");
            aa = inputFile.read(bb18);
            FileTime = new String(bb18, "UTF8");
            aa = inputFile.read(bb2);
            CRLF = new String(bb2, "UTF8");
            aa = inputFile.read(bb256);
            Comment = new String(bb256, "UTF8");

            aa = inputFile.read(bb4); // NumberOfCurves
            aa = inputFile.read(bb4); // BitsPerRecord
            BitsPerRecord = bb4[0];
            aa = inputFile.read(bb4); //  RoutingChannels
            RoutingChannels = bb4[0];
            aa = inputFile.read(bb4); //  NumberOfBoards
            NumberOfBoards = bb4[0];
            aa = inputFile.read(bb4); //  ActiveCurve
            ActiveCurve = bb4[0];
            aa = inputFile.read(bb4); //  MeasurementMode
            MeasurementMode = bb4[0];
            aa = inputFile.read(bb4); //  SubMode
            SubMode = bb4[0];
            aa = inputFile.read(bb4); //  RangeNo
            RangeNo = bb4[0];
            aa = inputFile.read(bb4); //  Offset
            Offset = byte4ToInt(bb4);
            aa = inputFile.read(bb4); //  AcquisitionTime
            aa = inputFile.read(bb4); //  StopAt
            StopAt = byte4ToInt(bb4);
            aa = inputFile.read(bb4); //  StopOnOvfl
            StopOnOvfl = bb4[0];
            aa = inputFile.read(bb4); //  Restart
            Restart = bb4[0];
            aa = inputFile.read(bb4); //  DisplayLinLog
            DisplayLinLog = bb4[0];
            aa = inputFile.read(bb4); //  DisplayTimeAxisFrom
            aa = inputFile.read(bb4); //  DisplayTimeAxisTo
            aa = inputFile.read(bb4); //  DisplayCountAxisFrom
            aa = inputFile.read(bb4); //  DisplayCountAxisTo

            aa = inputFile.read(bb4); //  DisplayCurve1MapTo
            aa = inputFile.read(bb4); //  DisplayCurve1Show
            aa = inputFile.read(bb4); //  DisplayCurve2MapTo
            aa = inputFile.read(bb4); //  DisplayCurve2Show
            aa = inputFile.read(bb4); //  DisplayCurve3MapTo
            aa = inputFile.read(bb4); //  DisplayCurve3Show
            aa = inputFile.read(bb4); //  DisplayCurve4MapTo
            aa = inputFile.read(bb4); //  DisplayCurve4Show
            aa = inputFile.read(bb4); //  DisplayCurve5MapTo
            aa = inputFile.read(bb4); //  DisplayCurve5Show
            aa = inputFile.read(bb4); //  DisplayCurve6MapTo
            aa = inputFile.read(bb4); //  DisplayCurve6Show
            aa = inputFile.read(bb4); //  DisplayCurve7MapTo
            aa = inputFile.read(bb4); //  DisplayCurve7Show
            aa = inputFile.read(bb4); //  DisplayCurve8MapTo
            aa = inputFile.read(bb4); //  DisplayCurve8Show

            aa = inputFile.read(bb4); //  Param1Start
            aa = inputFile.read(bb4); //  Param1Step
            aa = inputFile.read(bb4); //  Param1End
            aa = inputFile.read(bb4); //  Param2Start
            aa = inputFile.read(bb4); //  Param2Step
            aa = inputFile.read(bb4); //  Param2End
            aa = inputFile.read(bb4); //  Param3Start
            aa = inputFile.read(bb4); //  Param3Step
            aa = inputFile.read(bb4); //  Param3End

            aa = inputFile.read(bb4); //  RepeatMode
            aa = inputFile.read(bb4); //  RepeatsPerCurve
            aa = inputFile.read(bb4); //  RepeatTime
            aa = inputFile.read(bb4); //  RepeatWaitTime
            aa = inputFile.read(bb20); // ScriptName
            ScriptName = new String(bb20, "UTF8");
            
            for (int jj=0; jj<NumberOfBoards; jj++) { 
                aa = inputFile.read(bb60);
                //BrdHardwareIdent = new String(bb16, "UTF8");
                //System.out.println(BrdHardwareIdent);
            }

            int NumberOfRouterChannels=4;
            for (int jj=0; jj<NumberOfRouterChannels; jj++) { 
                aa = inputFile.read(bb4); // collective block
                //aa = inputFile.read(bb4); // collective block
                aa = inputFile.read(bb20); // collective block
            }

            aa = inputFile.read(bb4); //  ExtDevices
            aa = inputFile.read(bb4); //  Reserved1
            aa = inputFile.read(bb4); //  Reserved2
            aa = inputFile.read(bb4); //  InpRate0
            aa = inputFile.read(bb4); //  InpRate1
            aa = inputFile.read(bb4); //  StopAfter
            StopAfter = byte4ToInt(bb4);
            System.out.println(StopAfter);
            aa = inputFile.read(bb4); //  StopReason
            StopReason = bb4[0];
            aa = inputFile.read(bb4); //  NumRecords
            NumRecords = byte4ToInt(bb4);
            //System.out.println(NumRecords);
            aa = inputFile.read(bb4); //  ImgHdrSize
            ImgHdrSize = byte4ToInt(bb4);
            //System.out.println(ImgHdrSize);
            for (int jj=0; jj<ImgHdrSize; jj++) { aa = inputFile.read(bb4); }
            
            return 1;
        } finally {
            return 1;
        }
        
    }

    public int readHeaderPicoHarpPT3(FileInputStream inputFile) {
             
        int aa;
        int ii,jj;
        byte[] bb4;
        bb4 = new byte[4];

        int NumberOfBoards;
        int NumberOfRouterChannels=4;
        int ImgHdrSize;
        
        try {
            
            for (ii=0; ii<85; ii++) { aa = inputFile.read(bb4); };

            aa = inputFile.read(bb4); //  NumberOfBoards
            NumberOfBoards = bb4[0];

            for (ii=0; ii<48; ii++) { aa = inputFile.read(bb4); };

            for (jj=0; jj<NumberOfBoards; jj++) {
                for (ii=0; ii<15; ii++) { aa = inputFile.read(bb4); }
            }

            for (jj=0; jj<NumberOfRouterChannels; jj++) { 
                for (ii=0; ii<6; ii++) { aa = inputFile.read(bb4); }
            }

            for (ii=0; ii<8; ii++) { aa = inputFile.read(bb4); }

            aa = inputFile.read(bb4); //  ImgHdrSize
            ImgHdrSize = byte4ToInt(bb4);

            for (jj=0; jj<ImgHdrSize; jj++) { aa = inputFile.read(bb4); }
            
            return 1;
        } finally {
            return 1;
        }
        
    }

    public int updateDataPicoHarpPT3chA(FileInputStream inputFile) {

        byte byteconst15 = 15;
        int chNumber;
        int tempmta;
        
        nEvents=0;
        try {
            if (evA>=bufferInSizeA) { // read in the next chunk of data
                bufferInSizeA = inputFile.read(bufferInA);
                evA=0;
            }

            while ( nEvents<32768 && bufferInSizeA!=-1 ) {
                chNumber = (bufferInA[evA+3]>>4) & byteconst15; // 4 bits with channel number
                if ( chNumber==1) { 
                    if (bufferInA[evA]>=0) {
                        tempmta = bufferInA[evA];
                    } else {
                        tempmta = 256+bufferInA[evA];
                    }
                    if (bufferInA[evA+1]>=0) {
                        tempmta += 256*bufferInA[evA+1];
                    } else {
                        tempmta += 256*(256+bufferInA[evA+1]);
                    }                
                    iatA = tempmta - iatA;
                    photonsIat[nEvents]=iatA;
                    wAint[nEvents] = 1;
                    nEvents++;
                    iatA=tempmta; // iatA holds the last mta value now
                } else if (chNumber==15) {
                    iatA = iatA - 65535;
                }
                evA += 4;

                if (evA>=bufferInSizeA) { // read in the next chunk of data
                    bufferInSizeA = inputFile.read(bufferInA);
                    evA=0;
                }
            }
            return bufferInSizeA;
        } finally {
            return bufferInSizeA;
        }
    }
    
    public int updateDataPicoHarpPT3chB(FileInputStream inputFile) {

        byte byteconst15 = 15;
        int chNumber;
        int tempmta;
        
        nEvents=0;
        try {
            if (evA>=bufferInSizeA) { // read in the next chunk of data
                bufferInSizeA = inputFile.read(bufferInA);
                evA=0;
            }

            while ( nEvents<32768 && bufferInSizeA!=-1 ) {
                chNumber = (bufferInA[evA+3]>>4) & byteconst15; // 4 bits with channel number
                if ( chNumber==2) { 
                    if (bufferInA[evA]>=0) {
                        tempmta = bufferInA[evA];
                    } else {
                        tempmta = 256+bufferInA[evA];
                    }
                    if (bufferInA[evA+1]>=0) {
                        tempmta += 256*bufferInA[evA+1];
                    } else {
                        tempmta += 256*(256+bufferInA[evA+1]);
                    }                
                    iatA = tempmta - iatA;
                    photonsIat[nEvents]=iatA;
                    wAint[nEvents] = 1;
                    nEvents++;
                    iatA=tempmta; // iatA holds the last mta value now
                } else if (chNumber==15) {
                    iatA = iatA - 65535;
                }
                evA += 4;

                if (evA>=bufferInSizeA) { // read in the next chunk of data
                    bufferInSizeA = inputFile.read(bufferInA);
                    evA=0;
                }
            }
            return bufferInSizeA;
        } finally {
            return bufferInSizeA;
        }
    }

    public int updateDataPicoHarpPT3chAB(FileInputStream inputFile) {

        byte byteconst15 = 15;
        int chNumber;
        int tempmta;
        
        nEvents=0;
        try {
            if (evA>=bufferInSizeA) { // read in the next chunk of data
                bufferInSizeA = inputFile.read(bufferInA);
                evA=0;
            }

            while ( nEvents<32768 && bufferInSizeA!=-1 ) {
                chNumber = (bufferInA[evA+3]>>4) & byteconst15; // 4 bits with channel number
                if ( chNumber==1) { 
                    if (bufferInA[evA]>=0) {
                        tempmta = bufferInA[evA];
                    } else {
                        tempmta = 256+bufferInA[evA];
                    }
                    if (bufferInA[evA+1]>=0) {
                        tempmta += 256*bufferInA[evA+1];
                    } else {
                        tempmta += 256*(256+bufferInA[evA+1]);
                    }                
                    iatA = tempmta - iatA;
                    photonsIat[nEvents]=iatA;
                    wAint[nEvents] = 1;
                    wBint[nEvents] = 0;
                    nEvents++;
                    iatA=tempmta; // iatA holds the last mta value now
                } else if ( chNumber==2) { 
                    if (bufferInA[evA]>=0) {
                        tempmta = bufferInA[evA];
                    } else {
                        tempmta = 256+bufferInA[evA];
                    }
                    if (bufferInA[evA+1]>=0) {
                        tempmta += 256*bufferInA[evA+1];
                    } else {
                        tempmta += 256*(256+bufferInA[evA+1]);
                    }                
                    iatA = tempmta - iatA;
                    photonsIat[nEvents]=iatA;
                    wAint[nEvents] = 0;
                    wBint[nEvents] = 1;
                    nEvents++;
                    iatA=tempmta; // iatA holds the last mta value now
                } else if (chNumber==15) {
                    iatA = iatA - 65535;
                }
                evA += 4;

                if (evA>=bufferInSizeA) { // read in the next chunk of data
                    bufferInSizeA = inputFile.read(bufferInA);
                    evA=0;
                }
            }
            return bufferInSizeA;
        } finally {
            return bufferInSizeA;
        }
    }
    
    public int readHeaderPicoHarpPT2(FileInputStream inputFile) {
             
        int aa;
        int ii,jj;
        byte[] bb4;
        bb4 = new byte[4];

        int NumberOfBoards;
        int NumberOfRouterChannels=4;
        int ImgHdrSize;
        
        try {
            
            for (ii=0; ii<85; ii++) { aa = inputFile.read(bb4); };

            aa = inputFile.read(bb4); //  NumberOfBoards
            NumberOfBoards = bb4[0];

            for (ii=0; ii<48; ii++) { aa = inputFile.read(bb4); };

            for (jj=0; jj<NumberOfBoards; jj++) {
                for (ii=0; ii<15; ii++) { aa = inputFile.read(bb4); }
            }

            for (jj=0; jj<NumberOfRouterChannels; jj++) { 
                for (ii=0; ii<6; ii++) { aa = inputFile.read(bb4); }
            }

            for (ii=0; ii<8; ii++) { aa = inputFile.read(bb4); }

            aa = inputFile.read(bb4); //  ImgHdrSize
            ImgHdrSize = byte4ToInt(bb4);

            for (jj=0; jj<ImgHdrSize; jj++) { aa = inputFile.read(bb4); }
            
            return 1;
        } finally {
            return 1;
        }
        
    }

    public int updateDataPicoHarpPT2chA(FileInputStream inputFile) {

        byte byteconst15 = 15;
        int chNumber, time4bits;
        int tempmta;
        
        nEvents=0;
        try {
            if (evA>=bufferInSizeA) { // read in the next chunk of data
                bufferInSizeA = inputFile.read(bufferInA);
                evA=0;
            }

            while ( nEvents<32768 && bufferInSizeA!=-1 ) {
                chNumber = (bufferInA[evA+3]>>4) & byteconst15; // 4 bits with channel number
                if ( chNumber==1) { 
                    if (bufferInA[evA]>=0) {
                        tempmta = bufferInA[evA];
                    } else {
                        tempmta = 256+bufferInA[evA];
                    }
                    if (bufferInA[evA+1]>=0) {
                        tempmta += 256*bufferInA[evA+1];
                    } else {
                        tempmta += 256*(256+bufferInA[evA+1]);
                    }                
                    if (bufferInA[evA+2]>=0) {
                        tempmta += 256*256*bufferInA[evA+2];
                    } else {
                        tempmta += 256*256*(256+bufferInA[evA+2]);
                    }
                    tempmta += (bufferInA[evA+3] & byteconst15)*256*256*256;
                    iatA = tempmta - iatA;
                    photonsIat[nEvents]=iatA;
                    wAint[nEvents] = 1;
                    nEvents++;
                    iatA=tempmta; // iatA holds the last mta value now
                } else if (chNumber==15) {
                    iatA = iatA - 2147483647; // an overlow occurs after 128*256*256*256-1
                }
                evA += 4;

                if (evA>=bufferInSizeA) { // read in the next chunk of data
                    bufferInSizeA = inputFile.read(bufferInA);
                    evA=0;
                }
            }
            return bufferInSizeA;
        } finally {
            return bufferInSizeA;
        }
    }

    public void updateIntensityTraceOneChInt() {
        
        // intensity trace is calculated by dividing the event block
        // into chunks that contain binSize number of events each
        // and dividing the sum of all event weights by the time 
        // it took to collect these events
        
        int ev;
        int ii=0; // counter within one bin
        int x=0, y=0;
        
        ev=0;
        while (ev<nEvents && intensityTracePosition<32768) { // for all events in this block
        // note that the length of the intensity trace is limited to 32768 points
            y += wAint[ev];         // intensity value for this bin
            x += photonsIat[ev];    // duration of the current bin
            ev++;
            ii++;
            if (ii==binSize || ev==nEvents) {
                intensityTraceYA[intensityTracePosition] = y / (double)x;
                intensityTraceX[intensityTracePosition] = macroTime + x/2.0;
                macroTime += x;
                x = 0;
                y = 0;
                intensityTracePosition++;
                ii=0;
            }
        }
        
    } 
   
    public void updateIntensityTraceTwoChInt() {
        
        // intensity trace is calculated by dividing the event block
        // into chunks that contain binSize number of events each
        // and dividing the sum of all event weights by the time 
        // it took to collect these events
        
        int ev;
        int ii=0; // counter within one bin
        int x=0, yA=0, yB=0;
        
        ev=0;
        while (ev<nEvents && intensityTracePosition<32768) { // for all events in this block
        // note that the length of the intensity trace is limited to 32768 points
            yA += wAint[ev];         // intensity value for this bin
            yB += wBint[ev];         // intensity value for this bin
            x += photonsIat[ev];    // duration of the current bin
            ev++;
            ii++;
            if (ii==binSize || ev==nEvents) {
                intensityTraceYA[intensityTracePosition] = yA / (double)x;
                intensityTraceYB[intensityTracePosition] = yB / (double)x;
                intensityTraceX[intensityTracePosition] = macroTime + x/2.0;
                macroTime += x;
                x = 0;
                yA = 0;
                yB = 0;
                intensityTracePosition++;
                ii=0;
            }
        }
        
    } 
   
    public void updateCorrAutoInt() {

        int ii; // index that runs over all timelag values within one cascade
        int jj; // index that runs over all timebase values (cascades)
        int jjoffset; // helps to account for twice as many points in the first cascade
        int iioffset; // helps to account for twice as many points in the first cascade

        int[] timebase;         // timebase values
        timebase = new int[nc];
        int[] ntlag;
        ntlag = new int[nc];    // array with the number of points in each cascade
        int[] cc;
        cc = new int[(nc+1)*np];// lagtime values in units of timebase
        
        int tt;                 // temporary variable for storing lag time
        int iat;                // iat and
        int w;               // weight for the most recent photon
        int ev;                 // counter over all events
        byte st=-128;                 // index of the stop event
        double cum10, cum01, cum11, n, s10;
        double tlast_old;
        double wtotal_old;


        // A. calculate lagtime values
        ntlag[0] = np*2; // twice as many points in the first cascade
        for (jj=1; jj<nc; jj++) { ntlag[jj]=np; }
        tt = 0;
        timebase[0]= t0;
        for (jj=1; jj<nc; jj++) { timebase[jj] = timebase[jj-1] * 2; }
        for (jj=0; jj<nc; jj++) {
            if (jj==0) {jjoffset=0;} else {jjoffset=(jj+1)*np;}
            for (ii=0; ii<ntlag[jj]; ii++) {
                tt += timebase[jj];
                cc[jjoffset+ii] = tt / timebase[jj];
        }
    }
    // end of A

    // B. calculate the correlation
    tlast_old = tlast[0];
    wtotal_old = wAtotal[0];
        
    ev = 0;
    while (ev<nEvents) { // for all events in this block
        iat = photonsIat[ev];
        w = wAint[ev];
        ev++;
        for (jj=0; jj<nc; jj++) { // for all cascades
            if (jj==0) {jjoffset=0;} else {jjoffset=(jj+1)*np;} // use for cc, sXX array
            if (timebase[jj]!=1) { // if time coarsening is necessary
                // bin 1 - the bin that has been populated and is ready for correlation
                // bin 2 - the bin that is being populated
                mta[jj] += iat; // mta value of photons in bin 2
                if ( mta[jj]>=timebase[jj] ) { // bin 2 becomes bin 1
                    ciatarray[128+jj*256+sp[jj]] = iatnew[jj];
                    sp[jj]++; // array pointer for bin 1
                    wAintArray[128+jj*256+sp[jj]] = wAintTemp[jj];
                    wAintTemp[jj] = w;

                    // calculation of normalization values
                    tlast[jj] += iatnew[jj];
                    wAtotal[jj] += wAintArray[128+jj*256+sp[jj]];
                    if ( tlast[jj]<cc[jjoffset+ntlag[jj]-1] ) {
                        for (ii=ntlag[jj]-1; ii>=0; ii--) {
                            if (tlast[jj]<cc[jjoffset+ii]) {s01a[jjoffset+ii]+=wAintArray[128+jj*256+sp[jj]];}
                        }
                    }
                    // correlation search and calculation
                    st = sp[jj]; // check the last entry
                    st--;
                    if ( ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] ) {
                        while ( ciatarray[128+jj*256+st]<cc[jjoffset] ) {
                            st--;
                            ciatarray[128+jj*256+st] += iatnew[jj];
                        }
                        // while within the limits of timelag values
                        while ( ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] && ciatarray[128+jj*256+st]>=cc[jjoffset] ) {
                            s11aa[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st];
                            st--;
                            ciatarray[128+jj*256+st] += iatnew[jj];
                        }
                    }

                    iatnew[jj] = mta[jj]/timebase[jj]; // iat value for bin 2
                    mta[jj] = mta[jj]-iatnew[jj]*timebase[jj]; // offset mta value for bin 2
                } else { // processing the current time bin
                    wAintTemp[jj] += w; // add photon weights
                }
            } else { // timebase=1; no need to coarsen time
                ciatarray[128+jj*256+sp[jj]] = iat; // iat value for bin 1
                sp[jj]++; // array pointer for bin 1
                wAintArray[128+jj*256+sp[jj]] = w;

                // calculation of normalization values
                tlast[jj] += iat;
                wAtotal[jj] += wAintArray[128+jj*256+sp[jj]];
                if ( tlast[jj]<cc[jjoffset+ntlag[jj]-1] ) {
                    for (ii=ntlag[jj]-1; ii>=0; ii--) {
                        if (tlast[jj]<cc[jjoffset+ii]) {s01a[jjoffset+ii]+=wAintArray[128+jj*256+sp[jj]];}
                    }
                }

                // correlation search and calculation
                st = sp[jj]; // check the last entry
                st--;
                if ( ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] ) {
                    while ( ciatarray[128+jj*256+st]<cc[jjoffset] ) {
                        st--;
                        ciatarray[128+jj*256+st] += iatnew[jj];
                    }
                    // while within the limits of timelag values
                    while ( ciatarray[128+jj*256+st]>=cc[jjoffset] && ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] ) {
                        s11aa[jjoffset+ciatarray[128+jj*256+st]-1] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st];
                        st--;
                        ciatarray[128+jj*256+st] += iat;
                    }
                }
            }
        } // end of cycle over timebase values          
    } // end of cycle over all events
    // end of B

    // C. normalize the correlation function
    for (jj=0; jj<nc; jj++) {
        if (jj==0) {jjoffset=0;} else {jjoffset=(jj+1)*np;}
        for (ii=0; ii<ntlag[jj]; ii++) {
            s10 = wAtotal[jj] - wAintArray[128+jj*256+sp[jj]]; // subtract last bin
            st = sp[jj];
            st--;
            while ( (ciatarray[128+jj*256+st]+1)<=cc[jjoffset+ii] ) {
                  s10 = s10 - wAintArray[128+jj*256+st];
                  st--;
            }
            n = tlast[jj]-cc[jjoffset+ii];
            // basically, true s10 is s10 - wtotal
            cum10 = s10/n;
            cum01 = (wAtotal[jj] - s01a[jjoffset+ii])/n;
            cum11 = s11aa[jjoffset+ii]/(n-1.0) - s10/n*(wAtotal[jj] - s01a[jjoffset+ii])/(n-1.0);
            g11[jjoffset+ii] =  cum11/cum10/cum01; //s11aa[jjoffset+ii];
        }
    }
    // end of C

//    countrateA = wAtotal[0] / tlast[0] / dt;

    //return (wAtotal[0]-wtotal_old)/(tlast[0]-tlast_old)/dt; // local countrate
//////////////////////////////////////////////////////////////////////////////
        
    }

    public void updateCorrCrossInt() {

        int ii; // index that runs over all timelag values within one cascade
        int jj; // index that runs over all timebase values (cascades)
        int jjoffset; // helps to account for twice as many points in the first cascade
        int iioffset; // helps to account for twice as many points in the first cascade

        int[] timebase;         // timebase values
        timebase = new int[nc];
        int[] ntlag;
        ntlag = new int[nc];    // array with the number of points in each cascade
        int[] cc;
        cc = new int[(nc+1)*np];// lagtime values in units of timebase
        
        int tt;                 // temporary variable for storing lag time
        int iat;                // iat and
        int wA, wB;               // weight for the most recent photon
        int ev;                 // counter over all events
        byte st=-128;                 // index of the stop event
        double cum10a, cum01a, cum10b, cum01b, cum11aa, cum11bb, cum11ab, cum11ba, n, s10a, s10b;
        double tlast_old;
        double wAtotal_old, wBtotal_old;


        // A. calculate lagtime values
        ntlag[0] = np*2; // twice as many points in the first cascade
        for (jj=1; jj<nc; jj++) { ntlag[jj]=np; }
        tt = 0;
        timebase[0]= t0;
        for (jj=1; jj<nc; jj++) { timebase[jj] = timebase[jj-1] * 2; }
        for (jj=0; jj<nc; jj++) {
            if (jj==0) {jjoffset=0;} else {jjoffset=(jj+1)*np;}
            for (ii=0; ii<ntlag[jj]; ii++) {
                tt += timebase[jj];
                cc[jjoffset+ii] = tt / timebase[jj];
        }
    }
    // end of A

    // B. calculate the correlation
    tlast_old = tlast[0];
    wAtotal_old = wAtotal[0];
    wBtotal_old = wBtotal[0];
        
    ev = 0;
    while (ev<nEvents) { // for all events in this block
        iat = photonsIat[ev];
        wA = wAint[ev];
        wB = wBint[ev];
        ev++;
        for (jj=0; jj<nc; jj++) { // for all cascades
            if (jj==0) {jjoffset=0;} else {jjoffset=(jj+1)*np;} // use for cc, sXX array
            if (timebase[jj]!=1) { // if time coarsening is necessary
                // bin 1 - the bin that has been populated and is ready for correlation
                // bin 2 - the bin that is being populated
                mta[jj] += iat; // mta value of photons in bin 2
                if ( mta[jj]>=timebase[jj] ) { // bin 2 becomes bin 1
                    ciatarray[128+jj*256+sp[jj]] = iatnew[jj];
                    sp[jj]++; // array pointer for bin 1
                    wAintArray[128+jj*256+sp[jj]] = wAintTemp[jj];
                    wBintArray[128+jj*256+sp[jj]] = wBintTemp[jj];
                    wAintTemp[jj] = wA;
                    wBintTemp[jj] = wB;

                    // calculation of normalization values
                    tlast[jj] += iatnew[jj];
                    wAtotal[jj] += wAintArray[128+jj*256+sp[jj]];
                    wBtotal[jj] += wBintArray[128+jj*256+sp[jj]];
                    if ( tlast[jj]<cc[jjoffset+ntlag[jj]-1] ) {
                        for (ii=ntlag[jj]-1; ii>=0; ii--) {
                            if (tlast[jj]<cc[jjoffset+ii]) {s01a[jjoffset+ii]+=wAintArray[128+jj*256+sp[jj]];}
                            if (tlast[jj]<cc[jjoffset+ii]) {s01b[jjoffset+ii]+=wBintArray[128+jj*256+sp[jj]];}
                        }
                    }
                    // correlation search and calculation
                    st = sp[jj]; // check the last entry
                    st--;
                    if ( ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] ) {
                        while ( ciatarray[128+jj*256+st]<cc[jjoffset] ) {
                            st--;
                            ciatarray[128+jj*256+st] += iatnew[jj];
                        }
                        // while within the limits of timelag values
                        while ( ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] && ciatarray[128+jj*256+st]>=cc[jjoffset] ) {
                            if (wAintArray[128+jj*256+st]>0) {
                                if (wAintArray[128+jj*256+sp[jj]]>0) {
                                    s11aa[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+st]*wAintArray[128+jj*256+sp[jj]];
                                }
                                if (wBintArray[128+jj*256+sp[jj]]>0) {
                                    s11ab[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+st]*wBintArray[128+jj*256+sp[jj]];
                                }
                            }
                            if (wBintArray[128+jj*256+st]>0) {
                                if (wBintArray[128+jj*256+sp[jj]]>0) {
                                    s11bb[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wBintArray[128+jj*256+st]*wBintArray[128+jj*256+sp[jj]];
                                }
                                if (wAintArray[128+jj*256+sp[jj]]>0) {
                                    s11ba[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wBintArray[128+jj*256+st]*wAintArray[128+jj*256+sp[jj]];
                                }
                            }
                            st--;
                            ciatarray[128+jj*256+st] += iatnew[jj];
                        }
                    }

                    iatnew[jj] = mta[jj]/timebase[jj]; // iat value for bin 2
                    mta[jj] = mta[jj]-iatnew[jj]*timebase[jj]; // offset mta value for bin 2
                } else { // processing the current time bin
                    wAintTemp[jj] += wA; // add photon weights
                    wBintTemp[jj] += wB; // add photon weights
                }
            } else { // timebase=1; no need to coarsen time
                ciatarray[128+jj*256+sp[jj]] = iat; // iat value for bin 1
                sp[jj]++; // array pointer for bin 1
                wAintArray[128+jj*256+sp[jj]] = wA;
                wBintArray[128+jj*256+sp[jj]] = wB;

                // calculation of normalization values
                tlast[jj] += iat;
                wAtotal[jj] += wAintArray[128+jj*256+sp[jj]];
                wBtotal[jj] += wBintArray[128+jj*256+sp[jj]];
                if ( tlast[jj]<cc[jjoffset+ntlag[jj]-1] ) {
                    for (ii=ntlag[jj]-1; ii>=0; ii--) {
                        if (tlast[jj]<cc[jjoffset+ii]) {
                            s01a[jjoffset+ii]+=wAintArray[128+jj*256+sp[jj]];
                            s01b[jjoffset+ii]+=wBintArray[128+jj*256+sp[jj]];
                        }
                    }
                }

                // correlation search and calculation
                st = sp[jj]; // check the last entry
                st--;
                if ( ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] ) {
                    while ( ciatarray[128+jj*256+st]<cc[jjoffset] ) {
                        st--;
                        ciatarray[128+jj*256+st] += iatnew[jj];
                    }
                    // while within the limits of timelag values
                    while ( ciatarray[128+jj*256+st]>=cc[jjoffset] && ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] ) {
                        if (wAintArray[128+jj*256+st]>0) {
                            if (wAintArray[128+jj*256+sp[jj]]>0) {
                                s11aa[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+st]*wAintArray[128+jj*256+sp[jj]];
                            }
                            if (wBintArray[128+jj*256+sp[jj]]>0) {
                                s11ab[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+st]*wBintArray[128+jj*256+sp[jj]];
                            }
                        }
                        if (wBintArray[128+jj*256+st]>0) {
                            if (wBintArray[128+jj*256+sp[jj]]>0) {
                                s11bb[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wBintArray[128+jj*256+st]*wBintArray[128+jj*256+sp[jj]];
                            }
                            if (wAintArray[128+jj*256+sp[jj]]>0) {
                                s11ba[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wBintArray[128+jj*256+st]*wAintArray[128+jj*256+sp[jj]];
                            }
                        }
                        st--;
                        ciatarray[128+jj*256+st] += iat;
                    }
                }
            }
        } // end of cycle over timebase values          
    } // end of cycle over all events
    // end of B

    // C. normalize the correlation function
    for (jj=0; jj<nc; jj++) {
        if (jj==0) {jjoffset=0;} else {jjoffset=(jj+1)*np;}
        for (ii=0; ii<ntlag[jj]; ii++) {
            s10a = wAtotal[jj] - wAintArray[128+jj*256+sp[jj]]; // subtract last bin
            s10b = wBtotal[jj] - wBintArray[128+jj*256+sp[jj]]; // subtract last bin
            st = sp[jj];
            st--;
            while ( (ciatarray[128+jj*256+st]+1)<=cc[jjoffset+ii] ) {
                  s10a = s10a - wAintArray[128+jj*256+st];
                  s10b = s10b - wBintArray[128+jj*256+st];
                  st--;
            }
            n = tlast[jj]-cc[jjoffset+ii];
            // basically, true s10 is s10 - wtotal
            cum10a = s10a/n;
            cum01a = (wAtotal[jj] - s01a[jjoffset+ii])/n;
            cum10b = s10b/n;
            cum01b = (wBtotal[jj] - s01b[jjoffset+ii])/n;
            cum11aa = s11aa[jjoffset+ii]/(n-1.0) - s10a/n*(wAtotal[jj] - s01a[jjoffset+ii])/(n-1.0);
            cum11bb = s11bb[jjoffset+ii]/(n-1.0) - s10b/n*(wBtotal[jj] - s01b[jjoffset+ii])/(n-1.0);
            cum11ab = s11ab[jjoffset+ii]/(n-1.0) - s10a/n*(wBtotal[jj] - s01b[jjoffset+ii])/(n-1.0);
            cum11ba = s11ba[jjoffset+ii]/(n-1.0) - s10b/n*(wAtotal[jj] - s01a[jjoffset+ii])/(n-1.0);
            gAA[jjoffset+ii] =  cum11aa/cum10a/cum01a; //s11aa[jjoffset+ii];
            gBB[jjoffset+ii] =  cum11bb/cum10b/cum01b; //s11aa[jjoffset+ii];
            gAB[jjoffset+ii] =  cum11ab/cum10a/cum01b; //s11aa[jjoffset+ii];
            gBA[jjoffset+ii] =  cum11ba/cum10b/cum01a; //s11aa[jjoffset+ii];
        }
    }
    // end of C

//    countrateA = wAtotal[0] / tlast[0] / dt;

    //return (wAtotal[0]-wtotal_old)/(tlast[0]-tlast_old)/dt; // local countrate
//////////////////////////////////////////////////////////////////////////////
        
    }

    public void updateCorrAutoHOmltInt() {

        int ii; // index that runs over all timelag values within one cascade
        int jj; // index that runs over all timebase values (cascades)
        int jjoffset; // helps to account for twice as many points in the first cascade
        int iioffset; // helps to account for twice as many points in the first cascade

        int[] timebase;         // timebase values
        timebase = new int[nc];
        int[] ntlag;
        ntlag = new int[nc];    // array with the number of points in each cascade
        int[] cc;
        cc = new int[(nc+1)*np];// lagtime values in units of timebase
        
        int tt;                 // temporary variable for storing lag time
        int iat;                // iat and
        int w;               // weight for the most recent photon
        int ev;                 // counter over all events
        byte st=-128;                 // index of the stop event
        double cum10, cum01, cum11, n, s10, s20, s30, s01aT, s02aT, s03aT;
        double cum20, cum02, cum30, cum03, cum12, cum21, cum13, cum31, cum22;
        double tlast_old;
        double wtotal_old;


        // A. calculate lagtime values
        ntlag[0] = np*2; // twice as many points in the first cascade
        for (jj=1; jj<nc; jj++) { ntlag[jj]=np; }
        tt = 0;
        timebase[0]= t0;
        for (jj=1; jj<nc; jj++) { timebase[jj] = timebase[jj-1] * 2; }
        for (jj=0; jj<nc; jj++) {
            if (jj==0) {jjoffset=0;} else {jjoffset=(jj+1)*np;}
            for (ii=0; ii<ntlag[jj]; ii++) {
                tt += timebase[jj];
                cc[jjoffset+ii] = tt / timebase[jj];
        }
    }
    // end of A

    // B. calculate the correlation
    tlast_old = tlast[0];
    wtotal_old = wAtotal[0];
        
    ev = 0;
    while (ev<nEvents) { // for all events in this block
        iat = photonsIat[ev];
        w = wAint[ev];
        ev++;
        for (jj=0; jj<nc; jj++) { // for all cascades
            if (jj==0) {jjoffset=0;} else {jjoffset=(jj+1)*np;} // use for cc, sXX array
            if (timebase[jj]!=1) { // if time coarsening is necessary
                // bin 1 - the bin that has been populated and is ready for correlation
                // bin 2 - the bin that is being populated
                mta[jj] += iat; // mta value of photons in bin 2
                if ( mta[jj]>=timebase[jj] ) { // bin 2 becomes bin 1
                    ciatarray[128+jj*256+sp[jj]] = iatnew[jj];
                    sp[jj]++; // array pointer for bin 1
                    wAintArray[128+jj*256+sp[jj]] = wAintTemp[jj];
                    wAintTemp[jj] = w;

                    // calculation of normalization values
                    tlast[jj] += iatnew[jj];
                    wAtotal[jj]  += wAintArray[128+jj*256+sp[jj]];
                    wAtotal2[jj] += wAintArray[128+jj*256+sp[jj]]
                                   *wAintArray[128+jj*256+sp[jj]];
                    wAtotal3[jj] += wAintArray[128+jj*256+sp[jj]]
                                   *wAintArray[128+jj*256+sp[jj]]
                                   *wAintArray[128+jj*256+sp[jj]];
                    if ( tlast[jj]<cc[jjoffset+ntlag[jj]-1] ) {
                        for (ii=ntlag[jj]-1; ii>=0; ii--) {
                            if (tlast[jj]<cc[jjoffset+ii]) {
                                s01a[jjoffset+ii]+=wAintArray[128+jj*256+sp[jj]];
                                s02a[jjoffset+ii]+=wAintArray[128+jj*256+sp[jj]]
                                                  *wAintArray[128+jj*256+sp[jj]];
                                s03a[jjoffset+ii]+=wAintArray[128+jj*256+sp[jj]]
                                                  *wAintArray[128+jj*256+sp[jj]]
                                                  *wAintArray[128+jj*256+sp[jj]];
                            }
                        }
                    }
                    // correlation search and calculation
                    st = sp[jj]; // check the last entry
                    st--;
                    if ( ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] ) {
                        while ( ciatarray[128+jj*256+st]<cc[jjoffset] ) {
                            st--;
                            ciatarray[128+jj*256+st] += iatnew[jj];
                        }
                        // while within the limits of timelag values
                        while ( ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] && ciatarray[128+jj*256+st]>=cc[jjoffset] ) {
                            s11aa[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st];
                            s12aa[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st]
                                                                                    *wAintArray[128+jj*256+sp[jj]];
                            s21aa[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st]
                                                                                                                  *wAintArray[128+jj*256+st];
                            s13aa[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st]
                                                                                    *wAintArray[128+jj*256+sp[jj]]
                                                                                    *wAintArray[128+jj*256+sp[jj]];
                            s31aa[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st]
                                                                                                                  *wAintArray[128+jj*256+st]
                                                                                                                  *wAintArray[128+jj*256+st];
                            s22aa[jjoffset+ciatarray[128+jj*256+st]-cc[jjoffset]] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st]
                                                                                    *wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st];
                            st--;
                            ciatarray[128+jj*256+st] += iatnew[jj];
                        }
                    }

                    iatnew[jj] = mta[jj]/timebase[jj]; // iat value for bin 2
                    mta[jj] = mta[jj]-iatnew[jj]*timebase[jj]; // offset mta value for bin 2
                } else { // processing the current time bin
                    wAintTemp[jj] += w; // add photon weights
                }
            } else { // timebase=1; no need to coarsen time
                ciatarray[128+jj*256+sp[jj]] = iat; // iat value for bin 1
                sp[jj]++; // array pointer for bin 1
                wAintArray[128+jj*256+sp[jj]] = w;

                // calculation of normalization values
                tlast[jj] += iat;
                wAtotal[jj]  += wAintArray[128+jj*256+sp[jj]];
                wAtotal2[jj] += wAintArray[128+jj*256+sp[jj]]
                               *wAintArray[128+jj*256+sp[jj]];
                wAtotal3[jj] += wAintArray[128+jj*256+sp[jj]]
                               *wAintArray[128+jj*256+sp[jj]]
                               *wAintArray[128+jj*256+sp[jj]];
                if ( tlast[jj]<cc[jjoffset+ntlag[jj]-1] ) {
                    for (ii=ntlag[jj]-1; ii>=0; ii--) {
                        if (tlast[jj]<cc[jjoffset+ii]) {
                            s01a[jjoffset+ii]+=wAintArray[128+jj*256+sp[jj]];
                            s02a[jjoffset+ii]+=wAintArray[128+jj*256+sp[jj]]
                                              *wAintArray[128+jj*256+sp[jj]];
                            s03a[jjoffset+ii]+=wAintArray[128+jj*256+sp[jj]]
                                              *wAintArray[128+jj*256+sp[jj]]
                                              *wAintArray[128+jj*256+sp[jj]];
                        }
                    }
                }

                // correlation search and calculation
                st = sp[jj]; // check the last entry
                st--;
                if ( ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] ) {
                    while ( ciatarray[128+jj*256+st]<cc[jjoffset] ) {
                        st--;
                        ciatarray[128+jj*256+st] += iatnew[jj];
                    }
                    // while within the limits of timelag values
                    while ( ciatarray[128+jj*256+st]>=cc[jjoffset] && ciatarray[128+jj*256+st]<=cc[jjoffset+ntlag[jj]-1] ) {
                        s11aa[jjoffset+ciatarray[128+jj*256+st]-1] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st];
                        s12aa[jjoffset+ciatarray[128+jj*256+st]-1] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st]
                                                                     *wAintArray[128+jj*256+sp[jj]];
                        s21aa[jjoffset+ciatarray[128+jj*256+st]-1] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st]
                                                                                                   *wAintArray[128+jj*256+st];
                        s13aa[jjoffset+ciatarray[128+jj*256+st]-1] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st]
                                                                     *wAintArray[128+jj*256+sp[jj]]
                                                                     *wAintArray[128+jj*256+sp[jj]];
                        s31aa[jjoffset+ciatarray[128+jj*256+st]-1] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st]
                                                                                                   *wAintArray[128+jj*256+st]
                                                                                                   *wAintArray[128+jj*256+st];
                        s22aa[jjoffset+ciatarray[128+jj*256+st]-1] += wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st]
                                                                     *wAintArray[128+jj*256+sp[jj]]*wAintArray[128+jj*256+st];
                        st--;
                        ciatarray[128+jj*256+st] += iat;
                    }
                }
            }
        } // end of cycle over timebase values          
    } // end of cycle over all events
    // end of B

    // C. normalize the correlation function
    for (jj=0; jj<nc; jj++) {
        if (jj==0) {jjoffset=0;} else {jjoffset=(jj+1)*np;}
        for (ii=0; ii<ntlag[jj]; ii++) {
            s10 = wAtotal[jj] - wAintArray[128+jj*256+sp[jj]]; // subtract last bin
            s20 = wAtotal2[jj] - wAintArray[128+jj*256+sp[jj]]
                                *wAintArray[128+jj*256+sp[jj]]; // subtract last bin
            s30 = wAtotal3[jj] - wAintArray[128+jj*256+sp[jj]]
                                *wAintArray[128+jj*256+sp[jj]]
                                *wAintArray[128+jj*256+sp[jj]]; // subtract last bin
            st = sp[jj];
            st--;
            while ( (ciatarray[128+jj*256+st]+1)<=cc[jjoffset+ii] ) {
                  s10 = s10 - wAintArray[128+jj*256+st];
                  s20 = s20 - wAintArray[128+jj*256+st]*wAintArray[128+jj*256+st];
                  s30 = s30 - wAintArray[128+jj*256+st]*wAintArray[128+jj*256+st]*wAintArray[128+jj*256+st];
                  st--;
            }
            n = tlast[jj]-cc[jjoffset+ii];
            // basically, true s01 is wtotal - s01
            s01aT = wAtotal[jj]  - s01a[jjoffset+ii];
            s02aT = wAtotal2[jj] - s02a[jjoffset+ii];
            s03aT = wAtotal3[jj] - s03a[jjoffset+ii];
            cum10 = s10/n;
	    cum20 = s20/(n-1.0) - s10/n*s10/(n-1.0);
	    cum30 = n/(n-1.0)*s30/(n-2.0)
	              - 3*s20/(n-1.0)*s10/(n-2.0)
	              + 2*s10/n*s10/(n-1.0)*s10/(n-2.0);
            cum01 = s01aT/n;
	    cum02 = s02aT/(n-1.0) - s01aT/n*s01aT/(n-1.0);
	    cum03 = n/(n-1.0)*s03aT/(n-2.0)
	              - 3*s02aT/(n-1.0)*s01aT/(n-2.0)
	              + 2*s01aT/n*s01aT/(n-1.0)*s01aT/(n-2.0);
            cum11 = s11aa[jjoffset+ii]/(n-1.0) - s10/n*s01aT/(n-1.0);
            cum21 = n/(n-1.0)*s21aa[jjoffset+ii]/(n-2.0)-2*s10/(n-1.0)*s11aa[jjoffset+ii]/(n-2.0)
                  - s20/(n-1.0)*s01aT/(n-2.0)
                  + 2*s10/n*s10/(n-1.0)*s01aT/(n-2.0);
            cum12 = n/(n-1.0)*s12aa[jjoffset+ii]/(n-2.0)-2*s01aT/(n-1.0)*s11aa[jjoffset+ii]/(n-2.0)
                  - s02aT/(n-1.0)*s10/(n-2.0) + 2*s01aT/n*s01aT/(n-1.0)*s10/(n-2.0);
            cum31 = n/(n-1.0)*(n+1.0)/(n-2.0)*s31aa[jjoffset+ii]/(n-3.0)
                  - (n+1)/(n-1.0)*s30/(n-2.0)*s01aT/(n-3.0)
                  - 3*s11aa[jjoffset+ii]/(n-2.0)*s20/(n-3.0)
                  - 3*(n+1)/(n-1.0)*s21aa[jjoffset+ii]/(n-2.0)*s10/(n-3.0)
                  + 6*s11aa[jjoffset+ii]/(n-1.0)*s10/(n-2.0)*s10/(n-3.0)
                  + 6*s20/(n-1.0)*s10/(n-2.0)*s01aT/(n-3.0)
                  - 6*s01aT/n*s10/(n-1.0)*s10/(n-2)*s10/(n-3.0);
            cum13 = n/(n-1.0)*(n+1.0)/(n-2.0)*s13aa[jjoffset+ii]/(n-3.0)
                  - (n+1)/(n-1.0)*s03aT/(n-2.0)*s10/(n-3.0)
                  - 3*s11aa[jjoffset+ii]/(n-2.0)*s02aT/(n-3.0)
                  - 3*(n+1)/(n-1.0)*s12aa[jjoffset+ii]/(n-2.0)*s01aT/(n-3.0)
                  + 6*s11aa[jjoffset+ii]/(n-1.0)*s01aT/(n-2.0)*s01aT/(n-3.0)
                  + 6*s02aT/(n-1.0)*s01aT/(n-2.0)*s10/(n-3.0)
                  - 6*s10/n*s01aT/(n-1.0)*s01aT/(n-2)*s01aT/(n-3.0);
            cum22 = n/(n-1.0)*(n+1)/(n-2.0)*s22aa[jjoffset+ii]/(n-3.0)
                  - 2*(n+1)/(n-1.0)*s21aa[jjoffset+ii]/(n-2.0)*s01aT/(n-3.0)
                  - 2*(n+1)/(n-1.0)*s12aa[jjoffset+ii]/(n-2.0)*s10/(n-3.0)
                  - s20/(n-2.0)*s02aT/(n-3.0)
                  - 2*s11aa[jjoffset+ii]/(n-2.0)*s11aa[jjoffset+ii]/(n-3.0)
                  + 8*s11aa[jjoffset+ii]/(n-1.0)*s10/(n-2.0)*s01aT/(n-3.0)
                  + 2*s02aT/(n-1.0)*s10/(n-2.0)*s10/(n-3.0)
                  + 2*s20/(n-1.0)*s01aT/(n-2.0)*s01aT/(n-3.0)
                  - 6*s10/n*s10/(n-1.0)*s01aT/(n-2.0)*s01aT/(n-3.0);


            g11[jjoffset+ii] = cum11/cum10/cum01; //s11aa[jjoffset+ii]; //
            g12[jjoffset+ii] = (cum12 - cum11)/cum10/(cum02 - cum01); //s12aa[jjoffset+ii];// 
            g21[jjoffset+ii] = (cum21 - cum11)/cum01/(cum20 - cum10); //s21aa[jjoffset+ii];// 
            g13[jjoffset+ii] = (cum13 - 3*cum12 + 2*cum11)/cum10/(2*cum01 - 3*cum02 + cum03);//s13aa[jjoffset+ii];//
            g31[jjoffset+ii] = (cum31 - 3*cum21 + 2*cum11)/cum01/(2*cum10 - 3*cum20 + cum30);//s31aa[jjoffset+ii];//
            g22[jjoffset+ii] = (cum11 - cum12 - cum21 + cum22)/(cum20-cum10)/(cum02-cum01); //s22aa[jjoffset+ii];//

        }
    }
    // end of C
        
    }
    
}

