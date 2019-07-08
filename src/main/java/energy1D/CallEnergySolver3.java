/*
 * GNU GPL v3 License
 *
 * Copyright 2019 Niccolo` Tubini
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

package energy1D;

import java.util.ArrayList;
import java.util.HashMap;

import Richards1DSolver.Richards1DSolver3;
import oms3.annotations.*;
import physicalquantities.Variables;
import soilparameters.SoilParameters;



@Description("Solve the energy equation for the 1D domain coupled with Richards' equation.")
@Documentation("")
@Author(name = "Niccolo' Tubini, and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Hydrology, Energy, Infiltration")
@Bibliography("Casulli (2005)")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class CallEnergySolver3 {

	// SOIL PARAMETERS
	@Description("The hydraulic conductivity at saturation")
	@In 
	@Unit ("m/s")
	public double[] ks;

	@Description("Saturated water content")
	@In 
	@Unit ("-")
	public double[] thetaS;

	@Description("Residual water content")
	@In 
	@Unit ("-")
	public double[] thetaR;

	@Description("First parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par1SWRC;

	@Description("Second parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par2SWRC;

	@Description("Third parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par3SWRC;

	@Description("Fourth parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par4SWRC;

	@Description("Fifth parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par5SWRC;

	@Description("Critical value of psi for which the moisture capacity is null")
	@In 
	@Unit ("m")
	public double[] psiStar1;

	@Description("Critical value of psi for which the moisture capacity is null")
	@In 
	@Unit ("m")
	public double[] psiStar2;

	@Description("Critical value of psi for which the moisture capacity is null")
	@In 
	@Unit ("m")
	public double[] psiStar3;

	@Description("Aquitard compressibility")
	@In 
	@Unit ("1/Pa")
	public double[] alphaSpecificStorage;

	@Description("Water compressibility")
	@In 
	@Unit ("1/Pa")
	public double[] betaSpecificStorage;

	@Description("It is possibile to chose between 3 different models to compute "
			+ "the soil hydraulic properties: Van Genuchten; Brooks and Corey; Kosugi unimodal")
	@In 
	public String soilHydraulicModel;

	@Description("It is possible to choose among these models:"
			+ "notemperature, ....")
	@In 
	public String typeUHCTemperatureModel;

	@Description("It is possible to choose among these models:"
			+ "Mualem Van Genuchten, Mualem Brooks Corey, ....")
	@In 
	public String typeUHCModel;

	@Description("Hydraulic conductivity at control volume interface can be evaluated as"
			+ " the average of kappas[i] and kappas[i+1]"
			+ " the maximum between kappas[i] and kappas[i+1]"
			+ " the minimum between kappas[i] and kappas[i+1]"
			+ " a weighted average of kappas[i] and kappas[i+1] where weights are dx[i] and dx[i+1]")
	@In
	public String interfaceHydraulicCondType;
	
	@Description("Clay fraction")
	@In 
	@Unit ("-")
	public double[] clayFraction;

	@Description("Sand fraction")
	@In 
	@Unit ("-")
	public double[] sandFraction;
	
	@Description("Thermal conductivity of grain soils")
	@In 
	@Unit ("W/mK")
	public double[] grainThermalConductivity;
	
	@Description("Density of grain soils")
	@In 
	@Unit ("kg/m^3")
	public double[] grainDensity;
	
	@Description("Parameter 1 of thermal conductivity. Not used")
	@In 
	@Unit ("")
	public double[] par1;
	
	@Description("Parameter 2 of thermal conductivity. Not used")
	@In 
	@Unit ("")
	public double[] par2;
	
	@Description("Parameter 3 of thermal conductivity. Not used")
	@In 
	@Unit ("")
	public double[] par3;
	
	@Description("Parameter 4 of thermal conductivity. Not used")
	@In 
	@Unit ("")
	public double[] par4;
	
	@Description("Parameter 5 of thermal conductivity. Not used")
	@In 
	@Unit ("")
	public double[] par5;
	
	@Description("Parameter 6 of thermal conductivity. Not used")
	@In 
	@Unit ("")
	public double[] par6;
	
	@Description("Soil heat capacity model")
	@In
	public String soilHeatCapacityModel;
	
	@Description("Thermal conductivity model: up to now Johansen")
	@In
	public String soilThermalConductivityModel;
	
	@Description("Thermal conductivity at control volume interface can be evaluated as"
			+ " the average of lambdas[i] and lambdas[i+1]"
			+ " the maximum between lambdas[i] and lambdas[i+1]"
			+ " the minimum between lambdas[i] and lambdas[i+1]"
			+ " a weighted average of lambdas[i] and lambdas[i+1] where weights are dx[i] and dx[i+1]")
	@In
	public String interfaceThermalCondType;

	@Description("Number of Picard iteration to update the diffusive flux matrix")
	@In
	public int picardIteration=1;
	/////////////////////////////////////////////

	// ENERGY EQUATION

	
	@Description("Initial condition for temperature read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] temperatureIC;
	


	//////////////////////////////////////////////
	@Description("Coefficient to simulate ET by making use of Casulli's formula")
	@In
	@Unit("1/s")
	public double[] et;

	@Description("Initial condition for water head read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] psiIC;

	@Description("z coordinate read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] z;

	@Description("Space delta to compute gradients read from grid NetCDF file")
	@In 
	@Unit("m")
	public double[] spaceDeltaZ;

	@Description("Length of control volumes read from grid NetCDF file")
	@In 
	@Unit("m")
	public double[] deltaZ;

	@Description("Time amount at every time-loop")
	@In
	@Unit ("s")
	public double tTimestep;

	@Description("Time step of integration")
	@In
	@Unit ("s")
	public double timeDelta;

	@Description("Tolerance for Newton iteration")
	@In
	public double newtonTolerance;

	@Description("Control parameter for nested Newton algorithm:"
			+"0 --> simple Newton method"
			+"1 --> nested Newton method")
	@In
	public int nestedNewton; 

	@Description("Slope of the soil")
	@In 
	@Unit ("°")
	public double delta;

	// BOUNDARY CONDITIONS for RICHARDS

	@Description("The HashMap with the time series of the boundary condition at the top of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inRichardsTopBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann")
	@In 
	public String richardsTopBCType;

	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inRichardsBottomBC;

	@Description("It is possibile to chose among 3 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann"
			+ "- Impervious boundary condition --> Bottom Impervious")
	@In 
	public String richardsBottomBCType;

	// BOUNDARY CONDITIONS for ENERGY

	@Description("The HashMap with the time series of the boundary condition at the top of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inTemperatureTopBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann")
	@In 
	public String energyTopBCType;

	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inTemperatureBottomBC;

	@Description("It is possibile to chose among 3 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann"
			+ "- Impervious boundary condition --> Bottom Impervious")
	@In 
	public String energyBottomBCType;

	@Description("The HashMap with the time series of the short wave radiation")
	@In
	@Unit ("W/m^2")
	public HashMap<Integer, double[]> inShortWave;
	
	@Description("The HashMap with the time series of the wind velocity")
	@In
	@Unit ("m/s")
	public HashMap<Integer, double[]> inWind;
	
	
	@Description("The current date of the simulation.")
	@In
	public String inCurrentDate;

	@Description("Path of output files")
	@In
	public String dir;

	@Description("ArrayList of variable to be stored in the buffer writer for Richards' equation")
	@Out
	public ArrayList<double[]> richardsOutputToBuffer;

	@Description("ArrayList of variable to be stored in the buffer writer for energy equation")
	@Out
	public ArrayList<double[]> energyOutputToBuffer;

	//////////////////////////////////////////
	//////////////////////////////////////////

	@Description("Maximun number of Newton iterations")
	final int MAXITER_NEWT = 50;

	@Description("Richards top boundary condition according with topBCType")
	@Unit ("")
	double richardsTopBC;

	@Description("Richards bottom boundary condition according with bottomBCType")
	@Unit ("")
	double richardsBottomBC;

	@Description("Temperature top boundary condition according with topBCType")
	@Unit ("")
	double temperatureTopBC;

	@Description("Tempurature bottom boundary condition according with bottomBCType")
	@Unit ("")
	double temperatureBottomBC;
	
	@Description("Short wave radiation")
	@Unit ("")
	double shortWave;
	
	@Description("Wind velocity")
	@Unit ("")
	double wind;

	@Description("Number of control volume for domain discetrization")
	@Unit (" ")
	int NUM_CONTROL_VOLUMES; 

	@Description("It is needed to iterate on the date")
	int step;

	///////////////////////////////

	Energy1DSolver energySolver;
	Richards1DSolver3 richardsSolver;

	@Execute
	public void solve() {

		//System.out.println("RICHARDS 1D "+inCurrentDate);

		if(step==0){
//			energyTopBCType = "Top Dirichlet";
//			energyBottomBCType = "Bottom Dirichlet";
//			interfaceThermalCondType = "mean";
//			soilHeatCapacityModel = "Farouki";
//			soilThermalConductivityModel = "Johansen";
//			clayFraction = new double[thetaS.length];
//			sandFraction = new double[thetaS.length];
//			temperatureIC= new double[psiIC.length];
//			for(int i=0; i<thetaS.length; i++) {
//				clayFraction[i] = 0.5;
//				sandFraction[i] = 0.5;
//			}
//			for(int i=0; i<psiIC.length; i++) {
//				temperatureIC[i] = 280;
//			}
			
			NUM_CONTROL_VOLUMES = psiIC.length;
			Variables variables = Variables.getInstance(psiIC.clone(), temperatureIC.clone());
			SoilParameters soilParameters = SoilParameters.getInstance(alphaSpecificStorage, betaSpecificStorage, ks, 
					par1SWRC, par2SWRC, par3SWRC, par4SWRC, par5SWRC, psiStar1, psiStar2, psiStar3, thetaR, thetaS, clayFraction, sandFraction);

			richardsSolver = new Richards1DSolver3( soilHydraulicModel, typeUHCModel, typeUHCTemperatureModel, richardsTopBCType,
					richardsBottomBCType, interfaceHydraulicCondType, NUM_CONTROL_VOLUMES, nestedNewton, newtonTolerance,
					MAXITER_NEWT, picardIteration, tTimestep, deltaZ, spaceDeltaZ, psiIC);
			
			energySolver = new Energy1DSolver( soilHeatCapacityModel, soilThermalConductivityModel, energyTopBCType,
					energyBottomBCType, interfaceThermalCondType, NUM_CONTROL_VOLUMES, picardIteration,
					tTimestep, deltaZ, spaceDeltaZ, temperatureIC );

			richardsOutputToBuffer= new ArrayList<double[]>();
			energyOutputToBuffer= new ArrayList<double[]>();


		} // close step==0


		//time = time + tTimestep;

		richardsTopBC = 0.0;
		richardsTopBC = (inRichardsTopBC.get(1)[0]/1000)/tTimestep;

		richardsBottomBC = 0.0;
		if(inRichardsBottomBC != null)
			richardsBottomBC = inRichardsBottomBC.get(1)[0];
		if(richardsBottomBCType.equalsIgnoreCase("Bottom Neumann") || richardsBottomBCType.equalsIgnoreCase("BottomNeumann")) {
			richardsBottomBC = richardsBottomBC/tTimestep;
		}

		//temperatureTopBC = 288.0;
		temperatureTopBC = inTemperatureTopBC.get(1)[0]+273.15;
		
		//temperatureBottomBC = 3/100;
		temperatureBottomBC = inTemperatureBottomBC.get(1)[0];
		
		//shortWave = 0.0;
		shortWave = inShortWave.get(1)[0]*1000/tTimestep;
		
		wind = 0.0;
		//wind = inWind.get(1)[0];
		
		richardsOutputToBuffer.clear();
		energyOutputToBuffer.clear();
		
		double sumTimeDelta = 0;

		//double volume = 0.0;
		//double volumeNew = 0.0;
		while(sumTimeDelta < tTimestep) {

			if(sumTimeDelta + timeDelta>tTimestep) {
				timeDelta = tTimestep - sumTimeDelta;
			}
			sumTimeDelta = sumTimeDelta + timeDelta;
			
			richardsSolver.solve(richardsTopBC, richardsBottomBC, inCurrentDate, timeDelta);
			energySolver.solve(temperatureTopBC, temperatureBottomBC, richardsTopBC, shortWave, wind, inCurrentDate, timeDelta);

		}
		richardsOutputToBuffer.add(Variables.waterSuctions);
		richardsOutputToBuffer.add(Variables.thetas);
		richardsOutputToBuffer.add(psiIC);
		richardsOutputToBuffer.add(Variables.darcyVelocities);
		richardsOutputToBuffer.add(Variables.darcyVelocitiesCapillary);
		richardsOutputToBuffer.add(Variables.darcyVelocitiesGravity);
		richardsOutputToBuffer.add(Variables.poreVelocities);
		richardsOutputToBuffer.add(Variables.celerities);
		richardsOutputToBuffer.add(Variables.kinematicRatio);
		richardsOutputToBuffer.add(new double[] {Variables.errorVolume});
		richardsOutputToBuffer.add(new double[] {richardsTopBC*tTimestep*1000}); // I want to have rainfall height instead of water flux
		if(richardsBottomBCType.equalsIgnoreCase("Bottom Neumann") || richardsBottomBCType.equalsIgnoreCase("BottomNeumann")) {
			richardsBottomBC = richardsBottomBC*tTimestep;
		}
		richardsOutputToBuffer.add(new double[] {richardsBottomBC});
		richardsOutputToBuffer.add(new double[] {richardsTopBC+Variables.darcyVelocities[NUM_CONTROL_VOLUMES-1]}); // surface runoff
		
		energyOutputToBuffer.add(Variables.temperatures);
		energyOutputToBuffer.add(Variables.internalEnergies);
		energyOutputToBuffer.add(temperatureIC);
		energyOutputToBuffer.add(Variables.energyFluxes);
		energyOutputToBuffer.add(Variables.diffusionFluxes);
		energyOutputToBuffer.add(Variables.advectionFluxes);
		energyOutputToBuffer.add(Variables.pecletEnergy);
		energyOutputToBuffer.add(new double[] {Variables.errorEnergy});
		energyOutputToBuffer.add(new double[] {temperatureTopBC});
		energyOutputToBuffer.add(new double[] {shortWave}); //{shortWave}
		energyOutputToBuffer.add(new double[] {0.0}); //{windVelocity}
		energyOutputToBuffer.add(new double[] {temperatureBottomBC});
		
		
		step++;
	} //// MAIN CYCLE END ////

}  /// CLOSE Richards1d ///



