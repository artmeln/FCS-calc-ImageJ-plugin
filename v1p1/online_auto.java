import ij.*;
import ij.gui.*;
import ij.plugin.frame.PlugInFrame;

public class online_auto extends PlugInFrame {

	public  online_auto() {
		super(" online_auto");
	}

	public void run(String arg) {
		correlate();
		IJ.register( online_auto .class);
	}

	public void correlate() {

		// parameters of the correlation function calculation
		int t0 = 1;	// lowest time base value
		int nc = 12;	// number of cascades
		int np = 16;	// number of points in a cascade


		// provide your time resolution here		
		// *********************************
		double dt=10e-9; // (10 ns)
		// *********************************

		// populate the array of time values
		int ii, jj, jjoffset, tt, timebase, ntlag;
		double[] timeArray;
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

		// prepare the plot
		Plot plotG11 = new Plot("Auto","Time","Correlation");
		plotG11.setLogScaleX();
		plotG11.setLimits(1e-6, 0.1,  -0.01, 1.2);
		plotG11.setLineWidth(1);
		PlotWindow winG11 = plotG11.show();
		

		// create the correlator and set the calculation parameters
		Gmn cr;
		cr = new Gmn();
		cr.setParams(t0, nc,np);
		cr.initializeGmn("auto");


        	// initialize hardware and create a data buffer
		// ******** replace this block with custom code **************
	        byte[] datablock;
        	datablock = new byte[100000];
	        int datablocksize = 100000;
		int datablockposition=0;
		// ***********************************************************

		// run a measurement for 10 seconds
	        for (byte seconds=1; seconds<=5; seconds++) { // run for 5 seconds
            
            		// get data from the hardware
			// **** provide hardware specific calls here ****
			if (seconds==1) {
				for (ii=0; ii<datablocksize; ii++) { datablock[ii] = 5; } // fake data
			} else {
				for (ii=0; ii<datablocksize; ii++) { datablock[ii] = 7; } // fake data
			}
			// ***********************************************************
            	
			
			// calculate g11 for this data block
			while ( datablockposition != datablocksize ) {  // till end of block is reached

				// **** note that updateDataCustomHardwareOneCh must be defined by the user ****
				// convert a chunk of data to iat format
				datablockposition = cr.updateOnlineCustomHardwareOneCh(datablock, datablocksize); 
				// *****************************************************************************

				cr.updateIntensityTraceOneChInt(); // update intensity trace
				cr.updateCorrAutoInt(); // update correlation function

			}
			datablockposition = 0;
				
			IJ.wait(1000); // wait for 1 second			

			// update the plot
			plotG11 = new Plot("Auto","Time","Correlation");
			plotG11.setLogScaleX();
			plotG11.setLineWidth(1);
			plotG11.addPoints(timeArray, cr.getG11(), Plot.LINE);
			plotG11.setLimits(timeArray[0], timeArray[(nc+1)*np-1], cr.getMin(cr.getG11()), cr.getMax(cr.getG11()));
			winG11.drawPlot(plotG11);


	        }

		// plot intensity trace
		double[] traceX, traceYA, traceYB;
		Plot plotInt = new Plot("Intensity","Time","Countrate");
		PlotWindow winInt = plotInt.show();
		winInt.setLocationAndSize(20,370,550,100);
		traceX = cr.getIntensityTraceX();
		for (ii=0; ii<cr.getIntensityTracePosition(); ii++) { traceX[ii] = traceX[ii]*dt;}
		traceYA = cr.getIntensityTraceYA();
		for (ii=0; ii<cr.getIntensityTracePosition(); ii++) { traceYA[ii] = traceYA[ii]/dt;}
		plotInt.addPoints(traceX, traceYA, Plot.LINE);
		plotInt.setLimits(cr.getMin(traceX), cr.getMax(traceX), cr.getMin(traceYA), cr.getMax(traceYA));
		
	}
}
