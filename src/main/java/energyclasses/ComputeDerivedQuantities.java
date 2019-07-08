package energyclasses;
/*
 * GNU GPL v3 License
 *
 * Copyright 2017  Niccolo` Tubini
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

import interfaceconductivity.*;
import physicalquantities.Variables;
import soilparameters.SoilParameters;
import soilheatcapacity.SoilHeatCapacity;
import soilthermalconductivity.SoilThermalConductivity;
//import unsaturatedhydraulicconductivity.UnsaturatedHydraulicConductivity;

/**
 * This class computes derived quantities for Richards problem such as: 
 * 	- velocities at interfaces;
 *  - adimensional water content;
 *  - total water depth;
 *  - water volume
 *  
 *  @author Niccolo' Tubini
 */


public class ComputeDerivedQuantities {

	int NUM_CONTROL_VOLUMES;
	double kM;
	double kP;
	double bottomBC;
	double k_b;
	double internalEnergy;


	double[] spaceDelta;
	double[] dx;

	String bottomBCType;

	SoilHeatCapacity soilHeatCapacity;
	SoilThermalConductivity soilThermalConductivity;
	TotalDepth totalDepth;
	InterfaceConductivity interfaceThermalConductivity;
	Variables variables;
	SoilParameters soilParameters;


	/**
	 * @param NUM_CONTROL_VOLUMES number of control volumes
	 * @param dx vector containing control volume length
	 * @param spaceDelta vector containing distances between control volumes centroids
	 * @param par1SWRC vector containing first parameter of SWRC model
	 * @param par2SWRC vector containing second parameter of SWRC model
	 * @param thetaR vector containing thetaR
	 * @param thetaS vector containing thetaS
	 * @param soilPar is the class to compute the soil hydraulic properties
	 * @param totalDepth is the class to compute the total water depth
	 */
	public ComputeDerivedQuantities(int NUM_CONTROL_VOLUMES, double[] dx, double[] spaceDelta, SoilHeatCapacity soilHeatCapacity, 
			SoilThermalConductivity soilThermalConductivity, TotalDepth totalDepth, InterfaceConductivity interfaceThermalConductivity, String bottomBCType){

		this.NUM_CONTROL_VOLUMES = NUM_CONTROL_VOLUMES;
		this.spaceDelta = spaceDelta;
		this.dx = dx;
		this.soilHeatCapacity = soilHeatCapacity;
		System.out.println(this.soilHeatCapacity.toString());
		this.soilThermalConductivity = soilThermalConductivity;
		this.totalDepth = totalDepth;
		this.interfaceThermalConductivity = interfaceThermalConductivity;
		this.bottomBCType = bottomBCType; 
		//this.variables = Variables.getInstance();
		//this.soilParameters = SoilParameters.getInstance();
		
		//this.thetas = new double[NUM_CONTROL_VOLUMES];
		//this.volumes = new double[NUM_CONTROL_VOLUMES];
		//this.kappas = new double[NUM_CONTROL_VOLUMES];
		//this.velocities = new double[NUM_CONTROL_VOLUMES+1];
	}



	/**
	 * This method is used to update the following variables at each time step.
	 *
	 * @param bottomBC water head at the bottom
	 */
	public void setComputeDerivedQuantities(double bottomBC) {

		this.bottomBC = bottomBC;

	}



	/**
	 * This method computes the heat capacity at time level n
	 * 
	 */
	public void computeHeatCapacitiesOld() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			
			if(i == NUM_CONTROL_VOLUMES-1) {
				Variables.heatCapacitiesOld[i] = Variables.thetasOld[i] * 1000 * 4188;
			} else {
				Variables.heatCapacitiesOld[i] = soilHeatCapacity.cT(Variables.thetasOld[i], i);
			}
		}

	}	

	
	
	/**
	 * This method computes the heat capacity at time level n+1
	 * 
	 */
	public void computeHeatCapacities() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {

			if(i == NUM_CONTROL_VOLUMES-1) {
				Variables.heatCapacities[i] = Variables.thetas[i] * 1000 * 4188;
			} else {
				Variables.heatCapacities[i] = soilHeatCapacity.cT(Variables.thetas[i], i);
			}
		}

	}

	
	/**
	 * This method computes the thermal conductivity
	 * 
	 */
	public void computeLambdas() {
		
		//Variables.kappaBottom = unsaturatedHydraulicConductivity.hydraulicConductivity(bottomBC, 0);
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			
			if(i==NUM_CONTROL_VOLUMES-1) {
				Variables.lambdas[i] = soilThermalConductivity.thermalConductivity(i-1); // devo mettere waterSuctions[i-1]
			} else {
				Variables.lambdas[i] = soilThermalConductivity.thermalConductivity(i);
			}

		}
		//return kappas;	

	}



	/**
	 * This method computes the internal energy for each control volumes at time level n
	 * Note that for the 1D case the water volume is a length (volume per unit area)
	 * 
	 */
	public void computeInternalEnergyOld() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			
			if(i == NUM_CONTROL_VOLUMES-1) {
				Variables.internalEnergiesOld[i] = 1000.0*4188.0*(Variables.temperatures[i]-273.15)*Variables.thetasOld[i];
			} else {
				Variables.internalEnergiesOld[i] = Variables.heatCapacitiesOld[i]*(Variables.temperatures[i]-273.15)*this.dx[i];
			}
		}
		//return volumes;	

	}	


	
	/**
	 * This method computes the internal energy for each control volumes at time level n+1
	 * Note that for the 1D case the water volume is a length (volume per unit area)
	 * 
	 */
	public void computeInternalEnergyNew() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			
			if(i == NUM_CONTROL_VOLUMES-1) {
				Variables.internalEnergies[i] = 1000.0*4188.0*(Variables.temperatures[i]-273.15)*Variables.thetas[i];
			} else {
				Variables.internalEnergies[i] = Variables.heatCapacities[i]*(Variables.temperatures[i]-273.15)*this.dx[i];
			}
		}
		//return volumes;	

	}
	
	
	
	/**
	 * Compute the total internal energy at time level n
	 * @return internalEnergy
	 */
	public double computeTotalInternalEnergyOld() {

		internalEnergy = 0.0;
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			internalEnergy += Variables.internalEnergiesOld[i];
		}
		return internalEnergy;	

	}	
	
	
	

	/**
	 * Compute the total internal energy at time level n+1
	 * @return internalEnergy
	 */
	public double computeTotalInternalEnergy() {

		internalEnergy = 0.0;
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			internalEnergy += Variables.internalEnergies[i];
		}
		return internalEnergy;	

	}	



	/**
	 * FIXME
	 * This method computes energy flux at each control volume interface
	 * it is the sum of the diffusion flux and the convection flux (upwind)
	 * 
	 */
	public void computeEnergyFluxes() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {
				
				if (this.bottomBCType.equalsIgnoreCase("Bottom Neumann") || this.bottomBCType.equalsIgnoreCase("BottomNeumann")) {
					Variables.energyFluxes[i] = bottomBC;
				} else {
					kM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i],dx[i],dx[i]);
					Variables.energyFluxes[i] =  -kM * (Variables.temperatures[i]-bottomBC)/spaceDelta[i] + 1000*4188*( + 0.5*Variables.darcyVelocities[i]*(Variables.temperatures[i]-273.15+bottomBC-273.15) 
							- 0.5*Math.abs(Variables.darcyVelocities[i])*(Variables.temperatures[i]-273.15-bottomBC+273.15) );;
//					Variables.energyFluxes[i] =  -kM * (Variables.temperatures[i]-bottomBC)/spaceDelta[i] + 1000*4188*( + 0.5*Variables.darcyVelocities[i]*(Variables.temperatures[i]+bottomBC) 
//							- 0.5*Math.abs(Variables.darcyVelocities[i])*(Variables.temperatures[i]-bottomBC) );;
				}

			} else if(i == NUM_CONTROL_VOLUMES-1) {
				kP = Variables.lambdas[i-1];
				Variables.energyFluxes[i] =  -kP * (Variables.temperatures[i]-Variables.temperatures[i-1])/spaceDelta[i] + 1000*4188*( + 0.5*Variables.darcyVelocities[i]*(Variables.temperatures[i]-273.15+Variables.temperatures[i-1]-273.15) 
						- 0.5*Math.abs(Variables.darcyVelocities[i])*(Variables.temperatures[i]-273.15-Variables.temperatures[i-1]+273.15) );
//				Variables.energyFluxes[i] =  -kP * (Variables.temperatures[i]-Variables.temperatures[i-1])/spaceDelta[i] + 1000*4188*( + 0.5*Variables.darcyVelocities[i]*(Variables.temperatures[i]+Variables.temperatures[i-1]) 
//						- 0.5*Math.abs(Variables.darcyVelocities[i])*(Variables.temperatures[i]-Variables.temperatures[i-1]) );

			} else {
				kM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i-1],dx[i],dx[i-1]);
				Variables.energyFluxes[i] =  -kM * (Variables.temperatures[i]-Variables.temperatures[i-1])/spaceDelta[i] + 1000*4188*( + 0.5*Variables.darcyVelocities[i]*(Variables.temperatures[i]-273.15+Variables.temperatures[i-1]-273.15) 
						- 0.5*Math.abs(Variables.darcyVelocities[i])*(Variables.temperatures[i]-273.15-Variables.temperatures[i-1]+273.15) );
//				Variables.energyFluxes[i] =  -kM * (Variables.temperatures[i]-Variables.temperatures[i-1])/spaceDelta[i] + 1000*4188*( + 0.5*Variables.darcyVelocities[i]*(Variables.temperatures[i]+Variables.temperatures[i-1]) 
//						- 0.5*Math.abs(Variables.darcyVelocities[i])*(Variables.temperatures[i]-Variables.temperatures[i-1]) );

			}
		}

	}
	
	
	/**
	 * FIXME
	 * This method computes diffusion flux at each control volume interface
	 * 
	 */
	public void computeDiffusionFluxes() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if( i == 0 ) {
				
				if (this.bottomBCType.equalsIgnoreCase("Bottom Neumann") || this.bottomBCType.equalsIgnoreCase("BottomNeumann")) {
					Variables.diffusionFluxes[i] = bottomBC;
				} else {
					kM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i],dx[i],dx[i]);
					Variables.diffusionFluxes[i] =  -kM * (Variables.temperatures[i]-bottomBC)/spaceDelta[i];
				}

			} else if(i == NUM_CONTROL_VOLUMES-1) {
				kP = Variables.lambdas[i-1];
				Variables.diffusionFluxes[i] =  -kP * (Variables.temperatures[i]-Variables.temperatures[i-1])/spaceDelta[i];

			} else {
				kM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i-1],dx[i],dx[i-1]);
				Variables.diffusionFluxes[i] =  -kM * (Variables.temperatures[i]-Variables.temperatures[i-1])/spaceDelta[i];

			}
		}

	}
	
	
	/**
	 * FIXME
	 * This method computes convection flux at each control volume interface
	 * upwind formula
	 * 
	 */
	public void computeAdvectionFluxes() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			
			if( i == 0 ) {
				
				if (this.bottomBCType.equalsIgnoreCase("Bottom Neumann") || this.bottomBCType.equalsIgnoreCase("BottomNeumann")) {
					Variables.advectionFluxes[i] = Double.NaN;
				} else {
					Variables.advectionFluxes[i] =  1000*4188*( + 0.5*Variables.darcyVelocities[i]*(Variables.temperatures[i]-273.15+bottomBC-273.15) 
							- 0.5*Math.abs(Variables.darcyVelocities[i])*(Variables.temperatures[i]-273.15-bottomBC+273.15) );
//					Variables.advectionFluxes[i] =  1000*4188*( + 0.5*Variables.darcyVelocities[i]*(Variables.temperatures[i]+bottomBC) 
//							- 0.5*Math.abs(Variables.darcyVelocities[i])*(Variables.temperatures[i]-bottomBC) );
				}

			} else {
				Variables.advectionFluxes[i] =  1000*4188*( 0.5*Variables.darcyVelocities[i]*(Variables.temperatures[i]-273.15 + Variables.temperatures[i-1]-273.15) 
																		- 0.5*Math.abs(Variables.darcyVelocities[i])*(Variables.temperatures[i]-273.15 - Variables.temperatures[i-1]+273.15) );
//				Variables.advectionFluxes[i] =  1000*4188*( 0.5*Variables.darcyVelocities[i]*(Variables.temperatures[i] + Variables.temperatures[i-1]) 
//				- 0.5*Math.abs(Variables.darcyVelocities[i])*(Variables.temperatures[i] - Variables.temperatures[i-1]) );						

			}
		}

	}
	
	
	
	/**
	 * FIXME
	 * This method computes peclet number at each control volume interface
	 * 
	 */
	public void computePeclet() {

		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
							
			if (Variables.advectionFluxes[i] == 0) {
				
				Variables.pecletEnergy[i] = Double.NaN;
				
			} else {
				
				Variables.pecletEnergy[i] =  Variables.diffusionFluxes[i]/Variables.advectionFluxes[i];

			}
		}

	}




}






