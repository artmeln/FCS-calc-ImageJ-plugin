/*
 * Copyright (C) 2016 Artem Melnykov
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
    private int[] wAint;        // weights (channel A)
    private int[] wBint;        // weights (channel B)
    private int iatA;           // current inter-arrival time (ch. A)
    private int iatB;           // current inter-arrival time (ch. B)
    private int nEvents;      // number of events in the current photon array
    
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
    private double[] wBtotal;   // total number of photons for each cascade (channel B)
    private double[] s01a;      // right normalization constant (ch. A)
    private double[] s01b;      // right normalization constant (ch. B)
    private double[] s11aa;     // correlation function values (not normalized)
    private double[] s11bb;     // correlation function values (not normalized)
    private double[] s11ab;     // correlation function values (not normalized)
    private double[] s11ba;     // correlation function values (not normalized)
    private double countrateA;  // averagecountrate for all data (ch. A)
    private double countrateB;  // average countrate for all data (ch. B)
    
    // various time values
    private double dt;              // time resolution of the data (in seconds)
    private double bint;            // binning time for intensity traces
        
    public int initializeGmn(String correlationType) {
        int ff=0;
        try {
            if (correlationType=="auto") {
                evA = 32769;
                evB = 32769;
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

                return ff=1;
            } else if (correlationType=="high order") {
                evA = 32769;
                evB = 32769;
                bufferInA = new byte[32768];
                bufferInSizeA = 32768;
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
                tlast = new double[nc];
                wAtotal = new double[nc];
                s01a = new double[np*(nc+1)];
                s11aa = new double[np*(nc+1)];
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

}
