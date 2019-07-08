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

import bufferWriter.EnergyBuffer1D;
import bufferWriter.RichardsBuffer1D;
import energy1D.*;
import monodimensionalProblemTimeDependent.ReadNetCDFEnergyGrid1D;
import monodimensionalProblemTimeDependent.ReadNetCDFRichardsGrid1D;
import monodimensionalProblemTimeDependent.ReadNetCDFRichardsOutput1D;
import monodimensionalProblemTimeDependent.WriteNetCDFEnergy1D;
import monodimensionalProblemTimeDependent.WriteNetCDFRichards1D;

import org.junit.Test;

/**
 * Test the {@link TestEnergy1DSolver2} module.
 * 
 * 
 * @author Niccolo' Tubini
 */
public class TestEnergy1DSolver2 {

	@Test
	public void Test() throws Exception {


		String startDate = "2017-09-01 10:00";
		String endDate = "2018-12-31 00:00";
		
		int timeStepMinutes = 5;
		String fId = "ID";
		
//		String pathTopBCRichards ="resources/Input/TestAll_2.csv";
//		String pathBottomBCRichards ="resources/Input/TestAll_0.csv";
//		String pathGridRichards =  "resources/Input/Clay_noPonding_hydrostatic_VG.nc";
//		String pathOutputRichards = "resources/Output/Output_Clay_noPonding_hydrostatic_VG_b.nc";
		
		String pathTopBCRichards = "C:\\Users\\Niccolo\\Desktop\\EGU2019_notebook\\data\\Timeseries\\Rain_formatted.csv";
		String pathBottomBCRichards = "C:\\Users\\Niccolo\\Desktop\\EGU2019_notebook\\data\\Timeseries\\Rain_formatted.csv";
		String pathGridRichards =  "C:\\Users\\Niccolo\\Desktop\\EGU2019_notebook\\data\\//EGU_Soil_1_Richards.nc"; //Test_Richards.nc"; 
		String pathOutputRichards = "C:\\Users\\Niccolo\\Desktop\\EGU2019_notebook\\data\\Output_UnTrim_Richards_2_Temp_375.nc";
		String outputDescriptionRichards = "\n"
				//+ "K independent on T"
				+ "K dependent on T Ronan\n"
				+ "theta dependent on T beta =766\n"
				+ "Initial condition hydrostatic no ponding\n		"
				+ "Bottom free drainage\n		"
				//+ "Clay parameters VG: ks=0.000023m/s, alpha= 5.88m-1, n=1.16, thetaR=0.07, thetaS=0.5, alphaStorativity= 0 1/Pa, betaStorativity= 0 1/Pa\n		"
				//+ "Sand parameters VG: ks=0.003697m/s, alpha= 1.47m-1, n=1.7, thetaR=0.02, thetaS=0.38, alphaStorativity= 0 1/Pa, betaStorativity= 0 1/Pa\n		"
				+ "Grid input file: " + pathGridRichards +"\n		"
				+ "TopBC input file: " + pathTopBCRichards +"\n		"
				+ "BottomBC input file: " + pathBottomBCRichards +"\n		"
				+ "DeltaT: 300s\n		"
				+ "Picard iteration: 1\n		"
			    + "Interface k: max";
		
		String pathTopBCEnergy = "C:\\Users\\Niccolo\\Desktop\\EGU2019_notebook\\data\\Timeseries\\Temperature_5min_formatted.csv";
		String pathShortWave = "C:\\Users\\Niccolo\\Desktop\\EGU2019_notebook\\data\\Timeseries\\ShortWave_5min.csv";
		String pathWindVelocity = "C:\\Users\\Niccolo\\Desktop\\EGU2019_notebook\\data\\Timeseries\\Temperature_5min_formatted.csv";
		String pathBottomBCEnergy = "C:\\Users\\Niccolo\\Desktop\\EGU2019_notebook\\data\\Timeseries\\TemperatureBottom_5min_formatted.csv";
		String pathGridEnergy =  "C:\\Users\\Niccolo\\Desktop\\EGU2019_notebook\\data\\//EGU_Soil_1_Energy.nc"; //Test_Energy.nc"; //EGU_Soil_1_Energy.nc
		String pathOutputEnergy = "C:\\Users\\Niccolo\\Desktop\\EGU2019_notebook\\data\\Output_UnTrim_Energy_2_Temp_375.nc";
		String outputDescriptionEnergy = "\n"

				+ "Grid input file: " + pathGridEnergy +"\n		"
				+ "TopBC input file: " + pathTopBCEnergy +"\n		"
				+ "BottomBC input file: " + pathBottomBCEnergy +"\n		"
				+ "DeltaT: 300s\n		"
			    + "Interface k: max";
		
		
		OmsTimeSeriesIteratorReader topBCReaderRichards = getTimeseriesReader(pathTopBCRichards, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReaderRichards = getTimeseriesReader(pathBottomBCRichards, fId, startDate, endDate, timeStepMinutes);

		OmsTimeSeriesIteratorReader topBCReaderEnergy = getTimeseriesReader(pathTopBCEnergy, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReaderEnergy = getTimeseriesReader(pathBottomBCEnergy, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader shortWaveReader = getTimeseriesReader(pathShortWave, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader windReader = getTimeseriesReader(pathWindVelocity, fId, startDate, endDate, timeStepMinutes);

		RichardsBuffer1D bufferRichards = new RichardsBuffer1D();
		EnergyBuffer1D bufferEnergy = new EnergyBuffer1D();

		WriteNetCDFRichards1D writeNetCDFRichards = new WriteNetCDFRichards1D();
		WriteNetCDFEnergy1D writeNetCDFEnergy = new WriteNetCDFEnergy1D();

		ReadNetCDFRichardsGrid1D readNetCDFRichards = new ReadNetCDFRichardsGrid1D();
		ReadNetCDFEnergyGrid1D readNetCDFEnergy = new ReadNetCDFEnergyGrid1D();
		
		//Richards1DSolver R1DSolver = new Richards1DSolver();
		CallEnergySolver2 solver = new CallEnergySolver2();
		
		
		readNetCDFRichards.richardsGridFilename = pathGridRichards;
		readNetCDFRichards.read();
		
		readNetCDFEnergy.energyGridFilename = pathGridEnergy;
		readNetCDFEnergy.read();
		
		
		solver.z = readNetCDFRichards.z;
		solver.spaceDeltaZ = readNetCDFRichards.spaceDelta;
		solver.psiIC = readNetCDFRichards.psiIC;
		solver.temperatureIC = readNetCDFEnergy.temperatureIC;
		solver.deltaZ = readNetCDFRichards.deltaZ;
		solver.ks = readNetCDFRichards.Ks;
		solver.thetaS = readNetCDFRichards.thetaS;
		solver.thetaR = readNetCDFRichards.thetaR;
		solver.par1SWRC = readNetCDFRichards.par1SWRC;
		solver.par2SWRC = readNetCDFRichards.par2SWRC;
		solver.par3SWRC = readNetCDFRichards.par3SWRC;
		solver.par4SWRC = readNetCDFRichards.par4SWRC;
		solver.par5SWRC = readNetCDFRichards.par5SWRC;
		solver.psiStar1 = readNetCDFRichards.par6SWRC;
		solver.psiStar2 = readNetCDFRichards.par7SWRC;
		solver.psiStar3 = readNetCDFRichards.par8SWRC;
		solver.alphaSpecificStorage = readNetCDFRichards.alphaSS;
		solver.betaSpecificStorage = readNetCDFRichards.betaSS;
		solver.clayFraction = readNetCDFEnergy.clayFraction;
		solver.sandFraction = readNetCDFEnergy.sandFraction;
		solver.grainThermalConductivity = readNetCDFEnergy.grainThermalConductivity;
		solver.grainDensity = readNetCDFEnergy.grainDensity;
		solver.par1 = readNetCDFEnergy.par1;
		solver.par2 = readNetCDFEnergy.par2;
		solver.par3 = readNetCDFEnergy.par3;
		solver.par4 = readNetCDFEnergy.par4;
		solver.par5 = readNetCDFEnergy.par5;
		solver.par6 = readNetCDFEnergy.par6;
		solver.et = readNetCDFRichards.et;
		solver.soilHydraulicModel = "VanGenuchtenB";
		solver.typeUHCModel = "MualemVanGenuchten";
		solver.typeUHCTemperatureModel = "Ronan1998";  //"notemperature"; //    "notemperature"; //
		solver.interfaceHydraulicCondType = "max";
		solver.soilHeatCapacityModel = "Farouki";
		solver.soilThermalConductivityModel = "Johansen";
		solver.interfaceThermalCondType = "max";
		solver.richardsTopBCType = "Top Neumann";
		solver.richardsBottomBCType = "Bottom free drainage"; //"Bottom Dirichlet"; //
		solver.energyTopBCType = "Top Dirichlet";
		solver.energyBottomBCType = "Bottom Dirichlet";
		solver.delta = 0;
		solver.tTimestep = 300;
		solver.timeDelta = 300;
		solver.newtonTolerance = Math.pow(10,-11);
		solver.dir = "resources/Output";
		solver.nestedNewton =1;
		//R1DSolver.picardIteration = 1;
		while( topBCReaderRichards.doProcess  ) {
		
			
			topBCReaderRichards.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = topBCReaderRichards.outData;
			solver.inRichardsTopBC= bCValueMap;


			bottomBCReaderRichards.nextRecord();
			bCValueMap = bottomBCReaderRichards.outData;
			solver.inRichardsBottomBC = bCValueMap;
			
			topBCReaderEnergy.nextRecord();	
			bCValueMap = topBCReaderEnergy.outData;
			solver.inTemperatureTopBC= bCValueMap;
			
			bottomBCReaderEnergy.nextRecord();	
			bCValueMap = bottomBCReaderEnergy.outData;
			solver.inTemperatureBottomBC= bCValueMap;


			shortWaveReader.nextRecord();
			bCValueMap = shortWaveReader.outData;
			solver.inShortWave = bCValueMap;
			
			windReader.nextRecord();
			bCValueMap = windReader.outData;
			solver.inWind = bCValueMap;

			solver.inCurrentDate = topBCReaderRichards.tCurrent;
			
			solver.solve();

			
			bufferRichards.inputDate = solver.inCurrentDate;
			bufferRichards.inputSpatialCoordinate = readNetCDFRichards.eta;
			bufferRichards.inputDualSpatialCoordinate = readNetCDFRichards.etaDual;
			bufferRichards.inputVariable = solver.richardsOutputToBuffer;
			bufferRichards.solve();
			
			bufferEnergy.inputDate = solver.inCurrentDate;
			bufferEnergy.inputSpatialCoordinate = readNetCDFEnergy.eta;
			bufferEnergy.inputDualSpatialCoordinate = readNetCDFEnergy.etaDual;
			bufferEnergy.inputVariable = solver.energyOutputToBuffer;
			bufferEnergy.solve();
			
			writeNetCDFRichards.fileName = pathOutputRichards;
			writeNetCDFRichards.briefDescritpion = outputDescriptionRichards;
			writeNetCDFRichards.myVariables = bufferRichards.myVariable;
			writeNetCDFRichards.mySpatialCoordinate = bufferRichards.mySpatialCoordinate;
			writeNetCDFRichards.myDualSpatialCoordinate = bufferRichards.myDualSpatialCoordinate;			
			writeNetCDFRichards.doProcess = topBCReaderRichards.doProcess;
			writeNetCDFRichards.writeNetCDF();
			
			writeNetCDFEnergy.fileName = pathOutputEnergy;
			writeNetCDFEnergy.briefDescritpion = outputDescriptionEnergy;
			writeNetCDFEnergy.myVariables = bufferEnergy.myVariable;
			writeNetCDFEnergy.mySpatialCoordinate = bufferEnergy.mySpatialCoordinate;
			writeNetCDFEnergy.myDualSpatialCoordinate = bufferEnergy.myDualSpatialCoordinate;			
			writeNetCDFEnergy.doProcess = topBCReaderRichards.doProcess;
			writeNetCDFEnergy.writeNetCDF();


		}

		topBCReaderRichards.close();
		bottomBCReaderRichards.close();
		topBCReaderEnergy.close();
		bottomBCReaderEnergy.close();
		shortWaveReader.close();
		windReader.close();

				
		/*
		 * ASSERT 
		 */
//		ReadNetCDFRichardsOutput1D readTestData = new ReadNetCDFRichardsOutput1D();
//		readTestData.richardsOutputFilename = "resources/Output/test3.nc";
//		readTestData.read();
//		
//		ReadNetCDFRichardsOutput1D readSimData = new ReadNetCDFRichardsOutput1D();
//		readSimData.richardsOutputFilename = pathOutputRichards;
//		readSimData.read();
//		
//		assertEquals(readSimData.runOff.length, readTestData.runOff.length);
//		assertTrue("Runoff mismatch", Arrays.equals(readSimData.runOff,readTestData.runOff));
//		
//		assertEquals(readSimData.errorVolume.length, readTestData.errorVolume.length);
//		assertTrue("Runoff mismatch", Arrays.equals(readSimData.errorVolume,readTestData.errorVolume));
//	
		
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
