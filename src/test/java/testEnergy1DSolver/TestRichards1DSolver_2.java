/*
  * GNU GPL v3 License
 *
 * Copyright 2016 Niccolo` Tubini
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

package testEnergy1DSolver;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.net.URISyntaxException;
import java.util.*;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;

import bufferWriter.Buffer1D;
import energy1D.*;
import monodimensionalProblemTimeDependent.ReadNetCDFRichardsGrid1D;
import monodimensionalProblemTimeDependent.ReadNetCDFRichardsOutput1D;
import monodimensionalProblemTimeDependent.WriteNetCDFRichards1D;

import org.junit.Test;

/**
 * Test the {@link TestRichards1DSolver_2} module.
 * 
 * This test consider an initial hydrostatic condition with free drainage boundary condition at the 
 * bottom and no rain.
 * 
 * @author Niccolo' Tubini
 */
public class TestRichards1DSolver_2 {

	@Test
	public void Test() throws Exception {


		String startDate = "2017-01-01 00:00";
		String endDate = "2017-01-15 00:00";
		int timeStepMinutes = 5;
		String fId = "ID";
		
		String pathTopBC ="resources/Input/TestAll_0.csv";
		String pathBottomBC ="resources/Input/TestAll_0.csv";
		String pathGrid =  "resources/Input/Clay_noPonding_hydrostatic_VG.nc";
		String pathOutput = "resources/Output/Output_Clay_noPonding_hydrostatic_free_drainage_VG.nc";
		String outputDescription = "\n"
				+ "Initial condition hydrostatic no ponding\n		"
				+ "Bottom free drainage\n		"
				+ "Clay parameters VG: ks=0.000023m/s, alpha= 5.88m-1, n=1.16, thetaR=0.07, thetaS=0.5, alphaStorativity= 0 1/Pa, betaStorativity= 0 1/Pa\n		"
				+ "Grid input file: " + pathGrid +"\n		"
				+ "TopBC input file: " + pathTopBC +"\n		"
				+ "BottomBC input file: " + pathBottomBC +"\n		"
				+ "DeltaT: 50s\n		"
				+ "Picard iteration: 1\n		"
			    + "Interface k: mean";
		
		
		OmsTimeSeriesIteratorReader topBCReader = getTimeseriesReader(pathTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);

		Buffer1D buffer = new Buffer1D();
		WriteNetCDFRichards1D writeNetCDF = new WriteNetCDFRichards1D();
		ReadNetCDFRichardsGrid1D readNetCDF = new ReadNetCDFRichardsGrid1D();
		
		Richards1DSolver R1DSolver = new Richards1DSolver();
		
		
		
		readNetCDF.richardsGridFilename = pathGrid;
		
		readNetCDF.read();
		
		
		R1DSolver.z = readNetCDF.z;
		R1DSolver.spaceDeltaZ = readNetCDF.spaceDelta;
		R1DSolver.psiIC = readNetCDF.psiIC;
		R1DSolver.deltaZ = readNetCDF.deltaZ;
		R1DSolver.ks = readNetCDF.Ks;
		R1DSolver.thetaS = readNetCDF.thetaS;
		R1DSolver.thetaR = readNetCDF.thetaR;
		R1DSolver.par1SWRC = readNetCDF.par1SWRC;
		R1DSolver.par2SWRC = readNetCDF.par2SWRC;
		R1DSolver.par3SWRC = readNetCDF.par3SWRC;
		R1DSolver.par4SWRC = readNetCDF.par4SWRC;
		R1DSolver.par5SWRC = readNetCDF.par5SWRC;
		R1DSolver.psiStar1 = readNetCDF.par6SWRC;
		R1DSolver.psiStar2 = readNetCDF.par7SWRC;
		R1DSolver.psiStar3 = readNetCDF.par8SWRC;
		R1DSolver.alphaSpecificStorage = readNetCDF.alphaSS;
		R1DSolver.betaSpecificStorage = readNetCDF.betaSS;
		R1DSolver.et = readNetCDF.et;
		R1DSolver.soilHydraulicModel = "van genuchten";
		R1DSolver.interfaceHydraulicCondType = "mean";
		R1DSolver.topBCType = "Top Neumann";
		R1DSolver.bottomBCType = "Bottom free drainage";
		R1DSolver.delta = 0;
		R1DSolver.tTimestep = 300;
		R1DSolver.timeDelta = 50;
		R1DSolver.newtonTolerance = Math.pow(10,-11);
		R1DSolver.dir = "resources/Output";
		R1DSolver.nestedNewton =1;
		//R1DSolver.picardIteration = 1;
		while( topBCReader.doProcess  ) {
		
			
			topBCReader.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = topBCReader.outData;
			R1DSolver.inTopBC= bCValueMap;


			bottomBCReader.nextRecord();
			bCValueMap = bottomBCReader.outData;
			R1DSolver.inBottomBC = bCValueMap;

			R1DSolver.inCurrentDate = topBCReader.tCurrent;
			
			R1DSolver.solve();

			
			buffer.inputDate = R1DSolver.inCurrentDate;
			buffer.inputSpatialCoordinate = readNetCDF.eta;
			buffer.inputDualSpatialCoordinate = readNetCDF.etaDual;
			buffer.inputVariable = R1DSolver.outputToBuffer;
			
			buffer.solve();
			
			writeNetCDF.fileName = pathOutput;
			writeNetCDF.briefDescritpion = outputDescription;
			writeNetCDF.myVariables = buffer.myVariable;
			writeNetCDF.mySpatialCoordinate = buffer.mySpatialCoordinate;
			writeNetCDF.myDualSpatialCoordinate = buffer.myDualSpatialCoordinate;			
			writeNetCDF.doProcess = topBCReader.doProcess;
			writeNetCDF.writeNetCDF();


		}

		topBCReader.close();
		bottomBCReader.close();
				
		/*
		 * ASSERT 
		 */
		ReadNetCDFRichardsOutput1D readTestData = new ReadNetCDFRichardsOutput1D();
		readTestData.richardsOutputFilename = "resources/Output/test2.nc";
		readTestData.read();
		
		ReadNetCDFRichardsOutput1D readSimData = new ReadNetCDFRichardsOutput1D();
		readSimData.richardsOutputFilename = pathOutput;
		readSimData.read();
		
		assertEquals(readSimData.runOff.length, readTestData.runOff.length);
		assertTrue("Runoff mismatch", Arrays.equals(readSimData.runOff,readTestData.runOff));
		
	}

	private OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = timeStepMinutes;
		reader.tEnd = endDate;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}
}
