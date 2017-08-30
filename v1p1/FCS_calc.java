import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.frame.PlugInFrame;
import java.lang.*;

import ij.gui.Plot;
import java.awt.*;

import ij.util.Java2;
import javax.swing.JFileChooser;

import ij.io.*;
import java.io.*;

public class FCS_calc extends PlugInFrame {

	public  FCS_calc() {
		super(" FCS_calc");
	}

	public void run(String arg) {
		correlate();
		IJ.register(  FCS_calc .class);
	}

	public void correlate() {

		//aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
		// set up the correlator options 

		GenericDialog gd = new GenericDialog("Correlator settings");
		
		// correlator parameters
		gd.addNumericField("Time base",1,0);
		gd.addNumericField("N of cascades",20,0);
		gd.addNumericField("N of points per cascade",16,0);

		// data types supported
		String[] dataTypes;
		dataTypes = new String[4];
		dataTypes[0] = "Flex (8-bit)";
		dataTypes[1] = "Confocor2";
		dataTypes[2] = "Confocor3";
		dataTypes[3] = "PicoHarp pt3";
		gd.addChoice ("Data type", dataTypes, dataTypes[0]);
		
		// calculation types supported
		String[] calculationTypes;
		calculationTypes = new String[3];
		calculationTypes[0] = "auto";
		calculationTypes[1] = "cross";
		calculationTypes[2] = "autoHOmlt";
		gd.addChoice ("Calculation type", calculationTypes, calculationTypes[0]);

		// option for autoscaling the plots
		gd.addCheckbox("Autoscale", true);
		
		gd.showDialog();
		if (gd.wasCanceled()) { return; }

		// get values provided by the user
		double dblt0 = gd.getNextNumber();
		double  dblnc = gd.getNextNumber();
		double dblnp = gd.getNextNumber();
		int t0 = (int)dblt0;
		int nc = (int)dblnc;
		int np = (int)dblnp;
		String dataTypeSelected;
		dataTypeSelected = gd.getNextChoice();
		String calculationTypeSelected;
		calculationTypeSelected = gd.getNextChoice();
		boolean autoscale = gd.getNextBoolean(); 
		double minX=1e-6, maxX=0.1, minY=-0.01, maxY=1.2;
		double[] traceX, traceYA, traceYB;

		String chSelected = "";
		String legendAA="", legendBB="", legendAB="", legendBA="";
		int kk=0;
		int[] cB = {0, 255, 0, 0};
		int[] cG = {0, 0, 255, 0};
		int[] cR = {0, 0, 0, 255};
		FileInputStream inA = null;
		FileInputStream inB = null;
		File[] filesA = null;
		File[] filesB = null;
		//aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
		
		//bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
		// prepare a plot
		int ii, jj, jjoffset, tt, timebase, ntlag;
		double[] timeArray;
		double dt=1;
		//if ( dataTypeSelected=="Flex (8-bit)") { dt=16.6667e-9; }
		if ( dataTypeSelected=="Flex (8-bit)" ) { dt=12.5e-9; }
		if ( dataTypeSelected=="Confocor2" ) { dt=50e-9; }
		if ( dataTypeSelected=="Confocor3" ) { dt=50e-9; }
		if ( dataTypeSelected=="PicoHarp pt3" ) { dt=50e-9; }
		timeArray = new double[(nc+1)*np];
		tt = 0;
		timebase = t0;
		for (jj=0; jj<nc; jj++) {
			if (jj==0) { jjoffset=0; ntlag=2*np; } else { jjoffset=(jj+1)*np; ntlag=np; }
			for (ii=0; ii<ntlag; ii++) {
				tt += timebase;
				timeArray[jjoffset+ii] = tt*dt;
			}
			timebase = timebase * 2;
 		}
		Plot plotG11 = new Plot("Auto","Time","Correlation");
		plotG11.setLogScaleX();
		plotG11.setLimits(minX, maxX,  minY, maxY);
		plotG11.setLineWidth(1);
			
		Plot plotInt = new Plot("Intensity","Time","Countrate");
		Plot plotIntA = new Plot("Intensity (Ch. A)","Time","Countrate");
		Plot plotIntB = new Plot("Intensity (Ch. B)","Time","Countrate");

		Plot plotGaa = new Plot("Auto-gAA","Time","Correlation");
		plotGaa.setLogScaleX();
		plotGaa.setLimits(dt, timeArray[(nc+1)*np-1], minY, maxY);
		plotGaa.setLineWidth(1);
			
		Plot plotGbb = new Plot("Auto-gBB","Time","Correlation");
		plotGbb.setLogScaleX();
		plotGbb.setLimits(dt, timeArray[(nc+1)*np-1],  minY, maxY);
		plotGbb.setLineWidth(1);

		Plot plotGab = new Plot("Cross-gAB","Time","Correlation");
		plotGab.setLogScaleX();
		plotGab.setLimits(dt, timeArray[(nc+1)*np-1],  minY, maxY);
		plotGab.setLineWidth(1);

		Plot plotGba = new Plot("Cross-gBA","Time","Correlation");
		plotGba.setLogScaleX();
		plotGba.setLimits(dt, timeArray[(nc+1)*np-1],  minY, maxY);
		plotGba.setLineWidth(1);

		Plot plotG12 = new Plot("Auto-g12","Time","Correlation");
		plotG12.setLogScaleX();
		plotG12.setLimits(dt, timeArray[(nc+1)*np-1], minY, maxY);
		plotG12.setLineWidth(1);
			
		Plot plotG21 = new Plot("Auto-g21","Time","Correlation");
		plotG21.setLogScaleX();
		plotG21.setLimits(dt, timeArray[(nc+1)*np-1],  minY, maxY);
		plotG21.setLineWidth(1);

		Plot plotG13 = new Plot("Cross-g13","Time","Correlation");
		plotG13.setLogScaleX();
		plotG13.setLimits(dt, timeArray[(nc+1)*np-1],  minY, maxY);
		plotG13.setLineWidth(1);

		Plot plotG31 = new Plot("Cross-g31","Time","Correlation");
		plotG31.setLogScaleX();
		plotG31.setLimits(dt, timeArray[(nc+1)*np-1],  minY, maxY);
		plotG31.setLineWidth(1);

		Plot plotG22 = new Plot("Cross-g22","Time","Correlation");
		plotG22.setLogScaleX();
		plotG22.setLimits(dt, timeArray[(nc+1)*np-1],  minY, maxY);
		plotG22.setLineWidth(1);

		if (calculationTypeSelected=="auto") {
			PlotWindow winG11 = plotG11.show();
			winG11.setLocationAndSize(20,20,500,500);
			PlotWindow winInt = plotInt.show();
			winInt.setLocationAndSize(20,520,500,200);
		} else if (calculationTypeSelected=="cross") {
			PlotWindow winGaa = plotGaa.show();
			winGaa.setLocationAndSize(20,20,370,370);
			PlotWindow winGbb = plotGbb.show();
			winGbb.setLocationAndSize(20,400,370,370);
			PlotWindow winGab = plotGab.show();
			winGab.setLocationAndSize(400,20,370,370);
			PlotWindow winGba = plotGba.show();
			winGba.setLocationAndSize(400,400,370,370);
			PlotWindow winIntA = plotIntA.show();
			winIntA.setLocationAndSize(770,200,370,200);
			PlotWindow winIntB = plotIntB.show();
			winIntB.setLocationAndSize(770,400,370,200);
		} else if (calculationTypeSelected=="autoHOmlt") {
			PlotWindow winG11 = plotG11.show();
			winG11.setLocationAndSize(20,20,320,320);

			PlotWindow winG12 = plotG12.show();
			winG12.setLocationAndSize(350,20,320,320);

			PlotWindow winG21 = plotG21.show();
			winG21.setLocationAndSize(680,20,320,320);

			PlotWindow winG13 = plotG13.show();
			winG13.setLocationAndSize(20,350,320,320);

			PlotWindow winG31 = plotG31.show();
			winG31.setLocationAndSize(350,350,320,320);

			PlotWindow winG22 = plotG22.show();
			winG22.setLocationAndSize(680,350,320,320);

			PlotWindow winInt = plotInt.show();
			winInt.setLocationAndSize(1010,200,320,200);
		}
		//bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb

		//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		// set up the data stream
		if ( dataTypeSelected=="Flex (8-bit)" || dataTypeSelected == "Confocor2" || dataTypeSelected == "Confocor3" || dataTypeSelected == "PicoHarp pt3") { // open files
			Java2.setSystemLookAndFeel();
			JFileChooser fcA = new JFileChooser( OpenDialog.getLastDirectory() );
			fcA.setMultiSelectionEnabled(true);
			int returnVal = fcA.showOpenDialog(null);
			filesA = fcA.getSelectedFiles();
			File file0 = filesA[0];
			OpenDialog.setLastDirectory( file0.getAbsolutePath() );

			// select second channel data for cross-correlations
			if  ((dataTypeSelected=="Flex (8-bit)" || dataTypeSelected == "Confocor3") && calculationTypeSelected=="cross") {
				JFileChooser fcB = new JFileChooser( OpenDialog.getLastDirectory() );
				fcB.setMultiSelectionEnabled(true);
				fcB.showOpenDialog(null);
				filesB = fcB.getSelectedFiles();
				// check that the same number of files was selected
				if ( filesA.length != filesB.length ) { IJ.error("Same number of files must be selected for channels A and B!"); return; }
			}

			// select channel (A or B) in the Confocor2 data for auto-correlations
			if ( (dataTypeSelected=="Confocor2" && (calculationTypeSelected=="auto" || calculationTypeSelected=="autoHOmlt")) || 
			     (dataTypeSelected=="PicoHarp pt3" && (calculationTypeSelected=="auto" || calculationTypeSelected=="autoHOmlt")) ) {
				GenericDialog chSelector = new GenericDialog("Data channel");
				String[] chTypes;
				chTypes = new String[2];
				chTypes[0] = "A";
				chTypes[1] = "B";
				chSelector.addChoice ("Channel", chTypes, chTypes[1]);
				chSelector.showDialog();
				if (chSelector.wasCanceled()) { return; }
				chSelected = chSelector.getNextChoice();
			}

			double filesProgress=0, filesTotal;
			filesTotal = filesA.length;
			IJ.showProgress( filesProgress / filesTotal );	

			for (int ff=0; ff<filesA.length; ff++) { 
				// open data files: Channel A (or channel A and B)
				try {inA = new FileInputStream( filesA[ff].getAbsolutePath());}
				catch (Throwable e) {IJ.error("Unable to open"+filesA[ff].getAbsolutePath()); return;}
				// open data files: Channel B (if stored as a separate file)
				if  (dataTypeSelected=="Flex (8-bit)" && calculationTypeSelected=="cross") {
					try {inB = new FileInputStream( filesB[ff].getAbsolutePath());}
					catch (Throwable e) {IJ.error("Unable to open"+filesB[ff].getAbsolutePath()); return;}
				}

				// create the correlator and set the calculation parameters
				Gmn cr;
	        			cr = new Gmn();
				cr.setParams(t0, nc,np);
				
				if (dataTypeSelected=="Flex (8-bit)" && calculationTypeSelected=="auto") {
					cr.initializeGmn("auto");
					cr.updateDataFlexOneCh(inA);
					while (cr.getNofEvents()!=0) {  
						cr.updateCorrAutoInt();
						cr.updateIntensityTraceOneChInt();
						cr.updateDataFlexOneCh(inA);
            				}
					 legendAA = legendAA + filesA[ff].getName() + "\n"; 
				} else if (dataTypeSelected=="Confocor3" && calculationTypeSelected=="auto") {
					cr.initializeGmn("auto");
					cr.readHeaderConfocor3(inA);
					cr.updateDataConfocor3OneCh(inA);
					while (cr.getNofEvents()!=0) {  
						cr.updateCorrAutoInt();
						cr.updateIntensityTraceOneChInt();
						cr.updateDataConfocor3OneCh(inA);
            				}
					 legendAA = legendAA + filesA[ff].getName() + "\n"; 
				} else if (dataTypeSelected=="Confocor2" && calculationTypeSelected=="auto") {
					cr.initializeGmn("auto");
					if (chSelected=="A") {
						cr.readHeaderConfocor2(inA);
						cr.updateDataConfocor2chA(inA);
						while (cr.getNofEvents()!=0) {  
							cr.updateCorrAutoInt();
							cr.updateIntensityTraceOneChInt();
							cr.updateDataConfocor2chA(inA);
            					}
						legendAA = legendAA + filesA[ff].getName() + " - Ch.A\n"; 
					}					
					if (chSelected=="B") {
						cr.readHeaderConfocor2(inA);
						cr.updateDataConfocor2chB(inA);
						while (cr.getNofEvents()!=0) {  
							cr.updateCorrAutoInt();
							cr.updateIntensityTraceOneChInt();
							cr.updateDataConfocor2chB(inA);
            					}
						legendAA = legendAA + filesA[ff].getName() + " - Ch.B\n"; 
					}					
				} else if (dataTypeSelected=="PicoHarp pt3" && calculationTypeSelected=="auto") {
					cr.initializeGmn("auto");
					if (chSelected=="A") {
						cr.readHeaderPicoHarpPT3(inA);
						cr.updateDataPicoHarpPT3chA(inA);
						while (cr.getNofEvents()!=0) {  
							cr.updateCorrAutoInt();
							cr.updateIntensityTraceOneChInt();
							cr.updateDataPicoHarpPT3chA(inA);
            					}
						legendAA = legendAA + filesA[ff].getName() + " - Ch.A\n"; 
					}					
					if (chSelected=="B") {
						cr.readHeaderPicoHarpPT3(inA);
						cr.updateDataPicoHarpPT3chB(inA);
						while (cr.getNofEvents()!=0) {  
							cr.updateCorrAutoInt();
							cr.updateIntensityTraceOneChInt();
							cr.updateDataPicoHarpPT3chB(inA);
            					}
						legendAA = legendAA + filesA[ff].getName() + " - Ch.B\n"; 
					}					
				} else if (dataTypeSelected=="Flex (8-bit)" && calculationTypeSelected=="cross") {
					cr.initializeGmn("cross");
					cr.updateDataFlexTwoCh(inA, inB);
					while (cr.getNofEvents()!=0) {  
						cr.updateCorrCrossInt();
						cr.updateIntensityTraceTwoChInt();
						cr.updateDataFlexTwoCh(inA, inB);
            				}
					legendAA = legendAA + filesA[ff].getName() + "\n"; 
					legendBB = legendBB + filesB[ff].getName() + "\n"; 
					legendAB = legendAB + filesA[ff].getName() + " x " + filesB[ff].getName() + "\n"; 
					legendBA = legendBA + filesB[ff].getName() + " x " + filesA[ff].getName() + "\n"; 
				} else if (dataTypeSelected=="Confocor2" && calculationTypeSelected=="cross") {
					cr.initializeGmn("cross");
					cr.readHeaderConfocor2(inA);
					cr.updateDataConfocor2chAB(inA);
					while (cr.getNofEvents()!=0) {  
						cr.updateCorrCrossInt();
						cr.updateIntensityTraceTwoChInt();
						cr.updateDataConfocor2chAB(inA);
            				}
					legendAA = legendAA + filesA[ff].getName() + " - Ch.A\n"; 
					legendBB = legendBB + filesA[ff].getName() + " - Ch.B\n"; 
					legendAB = legendAB + filesA[ff].getName() + " - Ch.A x Ch.B\n"; 
					legendBA = legendBA + filesA[ff].getName() + " - Ch.B x Ch.A\n"; 
				} else if (dataTypeSelected=="PicoHarp pt3" && calculationTypeSelected=="cross") {
					cr.initializeGmn("cross");
					cr.readHeaderPicoHarpPT3(inA);
					cr.updateDataPicoHarpPT3chAB(inA);
					while (cr.getNofEvents()!=0) {  
						cr.updateCorrCrossInt();
						cr.updateIntensityTraceTwoChInt();
						cr.updateDataPicoHarpPT3chAB(inA);
            				}
					legendAA = legendAA + filesA[ff].getName() + " - Ch.A\n"; 
					legendBB = legendBB + filesA[ff].getName() + " - Ch.B\n"; 
					legendAB = legendAB + filesA[ff].getName() + " - Ch.A x Ch.B\n"; 
					legendBA = legendBA + filesA[ff].getName() + " - Ch.B x Ch.A\n"; 
				} else if (dataTypeSelected=="Flex (8-bit)" && calculationTypeSelected=="autoHOmlt") {
					cr.initializeGmn("autoHOmlt");
					cr.updateDataFlexOneCh(inA);
					while (cr.getNofEvents()!=0) {  
						cr.updateCorrAutoHOmltInt();
						cr.updateIntensityTraceOneChInt();
						cr.updateDataFlexOneCh(inA);
            				}
					 legendAA = legendAA + filesA[ff].getName() + "\n"; 
				} else if (dataTypeSelected=="Confocor3" && calculationTypeSelected=="autoHOmlt") {
					cr.initializeGmn("autoHOmlt");
					cr.readHeaderConfocor3(inA);
					cr.updateDataConfocor3OneCh(inA);
					while (cr.getNofEvents()!=0) {  
						cr.updateCorrAutoHOmltInt();
						cr.updateIntensityTraceOneChInt();
						cr.updateDataConfocor3OneCh(inA);
            				}
					 legendAA = legendAA + filesA[ff].getName() + "\n"; 
				} else if (dataTypeSelected=="Confocor2" && calculationTypeSelected=="autoHOmlt") {
					cr.initializeGmn("autoHOmlt");
					if (chSelected=="A") {
						cr.readHeaderConfocor2(inA);
						cr.updateDataConfocor2chA(inA);
						while (cr.getNofEvents()!=0) {  
							cr.updateCorrAutoHOmltInt();
							cr.updateIntensityTraceOneChInt();
							cr.updateDataConfocor2chA(inA);
            					}
						legendAA = legendAA + filesA[ff].getName() + " - Ch.A\n"; 
					}					
					if (chSelected=="B") {
						cr.readHeaderConfocor2(inA);
						cr.updateDataConfocor2chB(inA);
						while (cr.getNofEvents()!=0) {  
							cr.updateCorrAutoHOmltInt();
							cr.updateIntensityTraceOneChInt();
							cr.updateDataConfocor2chB(inA);
            					}
						legendAA = legendAA + filesA[ff].getName() + " - Ch.B\n"; 
					}					
				} else if (dataTypeSelected=="PicoHarp pt3" && calculationTypeSelected=="autoHOmlt") {
					cr.initializeGmn("autoHOmlt");
					if (chSelected=="A") {
						cr.readHeaderPicoHarpPT3(inA);
						cr.updateDataPicoHarpPT3chA(inA);
						while (cr.getNofEvents()!=0) {  
							cr.updateCorrAutoHOmltInt();
							cr.updateIntensityTraceOneChInt();
							cr.updateDataPicoHarpPT3chA(inA);
            					}
						legendAA = legendAA + filesA[ff].getName() + " - Ch.A\n"; 
					}					
					if (chSelected=="B") {
						cr.readHeaderPicoHarpPT3(inA);
						cr.updateDataPicoHarpPT3chB(inA);
						while (cr.getNofEvents()!=0) {  
							cr.updateCorrAutoHOmltInt();
							cr.updateIntensityTraceOneChInt();
							cr.updateDataPicoHarpPT3chB(inA);
            					}
						legendAA = legendAA + filesA[ff].getName() + " - Ch.B\n"; 
					}					
				} else {
					// not implemented
					return;
				}
		
				// plot correlation functions
				if (calculationTypeSelected=="auto") {
					Color clr = new Color( cB[kk], cG[kk], cR[kk]);
					if (kk>2) {kk=0;} else {kk++;}
					plotG11.setColor(clr);
					plotG11.addPoints(timeArray, cr.getG11(), Plot.LINE);
					if (autoscale  == true) {
						plotG11.setLimits(timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getG11()), cr.getMax(cr.getG11()));
					}
					plotG11.setColor("black");
					plotG11.addLegend(legendAA);
					plotG11.updateImage();
					plotInt.setColor(clr);
					traceX = cr.getIntensityTraceX();
					for (ii=0; ii<cr.getIntensityTracePosition(); ii++) { traceX[ii] = traceX[ii]*dt;}
					traceYA = cr.getIntensityTraceYA();
					for (ii=0; ii<cr.getIntensityTracePosition(); ii++) { traceYA[ii] = traceYA[ii]/dt;}
					plotInt.addPoints(traceX, traceYA, Plot.LINE);
					plotInt.setLimits(cr.getMin(traceX), cr.getMax(traceX), cr.getMin(traceYA), cr.getMax(traceYA));
					filesProgress = filesProgress+1;
					IJ.showProgress( filesProgress / filesTotal );
				} else if (calculationTypeSelected=="cross") {
					Color clr = new Color( cB[kk], cG[kk], cR[kk]);
					if (kk>2) {kk=0;} else {kk++;}
					plotGaa.setColor(clr);
					plotGaa.addPoints(timeArray, cr.getGaa(), Plot.LINE);
					if (autoscale  == true) {
						plotGaa.setLimits( timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getGaa()), cr.getMax(cr.getGaa()));
					}
					plotGaa.setColor("black");
					plotGaa.addLegend(legendAA);
					plotGaa.updateImage();
					plotGbb.setColor(clr);
					plotGbb.addPoints(timeArray, cr.getGbb(), Plot.LINE);
					if (autoscale  == true) {
						plotGbb.setLimits( timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getGbb()), cr.getMax(cr.getGbb()));
					}
					plotGbb.setColor("black");
					plotGbb.addLegend(legendBB);
					plotGbb.updateImage();
					plotGab.setColor(clr);
					plotGab.addPoints(timeArray, cr.getGab(), Plot.LINE);
					if (autoscale  == true) {
						plotGab.setLimits( timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getGab()), cr.getMax(cr.getGab()));
					}
					plotGab.setColor("black");
					plotGab.addLegend(legendAB);
					plotGab.updateImage();
					plotGba.setColor(clr);
					plotGba.addPoints(timeArray, cr.getGba(), Plot.LINE);
					if (autoscale  == true) {
						plotGba.setLimits( timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getGba()), cr.getMax(cr.getGba()));
					}
					plotGba.setColor("black");
					plotGba.addLegend(legendBA);
					plotGba.updateImage();

					plotIntA.setColor(clr);
					traceX = cr.getIntensityTraceX();
					for (ii=0; ii<cr.getIntensityTracePosition(); ii++) { traceX[ii] = traceX[ii]*dt;}
					traceYA = cr.getIntensityTraceYA();
					for (ii=0; ii<cr.getIntensityTracePosition(); ii++) { traceYA[ii] = traceYA[ii]/dt;}
					plotIntA.addPoints(traceX, traceYA, Plot.LINE);
					plotIntA.setLimits(cr.getMin(traceX), cr.getMax(traceX), cr.getMin(traceYA), cr.getMax(traceYA));
					
					plotIntB.setColor(clr);
					traceYB = cr.getIntensityTraceYB();
					for (ii=0; ii<cr.getIntensityTracePosition(); ii++) { traceYB[ii] = traceYB[ii]/dt;}
					plotIntB.addPoints(traceX, traceYB, Plot.LINE);
					plotIntB.setLimits(cr.getMin(traceX), cr.getMax(traceX), cr.getMin(traceYB), cr.getMax(traceYB));
					filesProgress = filesProgress+1;
					IJ.showProgress( filesProgress / filesTotal );
				} else if (calculationTypeSelected=="autoHOmlt") {
					Color clr = new Color( cB[kk], cG[kk], cR[kk]);
					if (kk>2) {kk=0;} else {kk++;}

					plotG11.setColor(clr);
					plotG11.addPoints(timeArray, cr.getG11(), Plot.LINE);
					if (autoscale  == true) {
						plotG11.setLimits(timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getG11()), cr.getMax(cr.getG11()));
					}
					plotG11.setColor("black");
					plotG11.addLegend(legendAA);
					plotG11.updateImage();
					plotInt.setColor(clr);

					plotG12.setColor(clr);
					plotG12.addPoints(timeArray, cr.getG12(), Plot.LINE);
					if (autoscale  == true) {
						plotG12.setLimits(timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getG12()), cr.getMax(cr.getG12()));
					}
					plotG12.setColor("black");
					plotG12.addLegend(legendAA);
					plotG12.updateImage();

					plotG21.setColor(clr);
					plotG21.addPoints(timeArray, cr.getG21(), Plot.LINE);
					if (autoscale  == true) {
						plotG21.setLimits(timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getG21()), cr.getMax(cr.getG21()));
					}
					plotG21.setColor("black");
					plotG21.addLegend(legendAA);
					plotG21.updateImage();

					plotG13.setColor(clr);
					plotG13.addPoints(timeArray, cr.getG13(), Plot.LINE);
					if (autoscale  == true) {
						plotG13.setLimits(timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getG13()), cr.getMax(cr.getG13()));
					}
					plotG13.setColor("black");
					plotG13.addLegend(legendAA);
					plotG13.updateImage();

					plotG31.setColor(clr);
					plotG31.addPoints(timeArray, cr.getG31(), Plot.LINE);
					if (autoscale  == true) {
						plotG31.setLimits(timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getG31()), cr.getMax(cr.getG31()));
					}
					plotG31.setColor("black");
					plotG31.addLegend(legendAA);
					plotG31.updateImage();

					plotG22.setColor(clr);
					plotG22.addPoints(timeArray, cr.getG22(), Plot.LINE);
					if (autoscale  == true) {
						plotG22.setLimits(timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getG22()), cr.getMax(cr.getG22()));
					}
					plotG22.setColor("black");
					plotG22.addLegend(legendAA);
					plotG22.updateImage();

					plotInt.setColor(clr);
					traceX = cr.getIntensityTraceX();
					for (ii=0; ii<cr.getIntensityTracePosition(); ii++) { traceX[ii] = traceX[ii]*dt;}
					traceYA = cr.getIntensityTraceYA();
					for (ii=0; ii<cr.getIntensityTracePosition(); ii++) { traceYA[ii] = traceYA[ii]/dt;}
					plotInt.addPoints(traceX, traceYA, Plot.LINE);
					plotInt.setLimits(cr.getMin(traceX), cr.getMax(traceX), cr.getMin(traceYA), cr.getMax(traceYA));
					filesProgress = filesProgress+1;
					IJ.showProgress( filesProgress / filesTotal );
				}
	
				try { inA.close(); }
				catch (Throwable e) {IJ.error("Stream does not exist"); return;}
				if (dataTypeSelected=="Flex (8-bit)" && calculationTypeSelected=="cross") {
					try { inB.close(); }
					catch (Throwable e) {IJ.error("Stream does not exist"); return;}
				}
				
			}
		} else {		// not implemented
			return;
		}
		
	}
}
