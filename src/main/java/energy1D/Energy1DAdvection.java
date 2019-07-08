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

import energyboundaryconditions.BoundaryCondition;
import energyboundaryconditions.SimpleBoundaryConditionFactory;
import energyclasses.*;
import interfaceconductivity.*;
import oms3.annotations.*;
import physicalquantities.Variables;
import soilheatcapacity.*;
import soilthermalconductivity.*;



@Description("Solve the energy equation for the 1D domain.")
@Documentation("")
@Author(name = "Niccolo' Tubini, and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Hydrology, energy, soil")
//@Bibliography("Casulli (2010)")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class Energy1DAdvection {


	//	@Description("Number of Picard iteration to update the diffusive flux matrix")
	public int picardIteration=1;

	//	@Description("Initial condition for water head read from grid NetCDF file")
	//	@In
	//	@Unit("m")
	public double[] temperatureIC;

	@Description("Time amount at every time-loop")
	@Unit ("s")
	public double tTimestep;
	//
	@Description("Time step of integration")
	@Unit ("s")
	public double timeDelta;
	//
	@Description("Tolerance for Newton iteration")
	//	@In
	public double newtonTolerance;
	//
	//	@Description("Control parameter for nested Newton algorithm:"
	//			+"0 --> simple Newton method"
	//			+"1 --> nested Newton method")
	//	@In
	public int nestedNewton; 
	//
	//	@Description("Slope of the soil")
	//	@In 
	//	@Unit ("°")
	public double delta;
	//
	//	// BOUNDARY CONDITIONS

	//	@Description("It is possibile to chose between 2 different kind "
	//			+ "of boundary condition at the top of the domain: "
	//			+ "- Dirichlet boundary condition --> Top Dirichlet"
	//			+ "- Neumann boundary condition --> Top Neumann")
	//	@In 
	public String topBCType;

	//	@Description("It is possibile to chose among 3 different kind "
	//			+ "of boundary condition at the bottom of the domain: "
	//			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
	//			+ "- Neumann boundary condition --> Bottom Neumann"
	//			+ "- Impervious boundary condition --> Bottom Impervious")
	//	@In 
	public String bottomBCType;

	//////////////////////////////////////////
	//////////////////////////////////////////


	@Description("Top boundary condition according with topBCType")
	@Unit ("")
	double topBC;

	@Description("Bottom boundary condition according with bottomBCType")
	@Unit ("")
	double bottomBC;

	@Description("Rainfall at the top of the domain")
	@Unit ("")
	double rain;

	@Description("Total internal energy at time level n")
	@Unit ("-")
	double internalEnergy;

	@Description("Total internal energy at time level n+1")
	@Unit ("-")
	double internalEnergyNew;

	@Description("Volume error between time levels n+1 and n")
	@Unit ("-")
	double errorEnergy;

	@Description("Vector collects the lower diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] lowerDiagonal;

	@Description("Vector collects the main diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] mainDiagonal;

	@Description("Vector collects the upper diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] upperDiagonal;

	@Description("Right hand side vector of the scalar equation to solve")
	@Unit ("-")
	double[] rhss;

	@Description("Thermal conductivity at the cell interface i+1/2")
	@Unit ("m/s")
	double lambdaP;

	@Description("Thermal conductivity at the cell interface i-1/2")
	@Unit ("m/s")
	double lambdaM;

	@Description("Number of control volume for domain discetrization")
	@Unit (" ")
	int NUM_CONTROL_VOLUMES; 

	@Description("Control variable for the integration time loop ")
	@Unit ("s")
	public double sumTimeDelta = 0.0;

	@Description("Space step")
	@Unit ("m")
	double[] spaceDelta;

	@Description("Vector containing the length of each control volume")
	double[] dx;

	@Description("Object to perform the Thomas algortithm")
	Thomas thomasAlg;

	@Description("Object dealing with soil heat capacity parametrization")
	SoilHeatCapacity soilHeatCapacity;
	SoilHeatCapacityFactory soilHeatCapacityFactory;;

	@Description("Object dealing with thermal conductivity")
	SoilThermalConductivity soilThermalConductivity;
	SimpleSoilThermalConductivityFactory soilThermalConductivityFactory;

	@Description("Object to compute total water depth")
	TotalDepth totalDepth;

	@Description("This object compute the diagonal and right hand side entries for the uppermost cell accordingly with the prescribed top boundary condition.")
	BoundaryCondition topBoundaryCondition;

	@Description("This object compute the diagonal and right hand side entries for the lowermost cell accordingly with the prescribed bottom boundary condition.")
	BoundaryCondition bottomBoundaryCondition;

	@Description("This object compute hydraulic quantities.")
	ComputeDerivedQuantities compute;

	@Description("This object compute the interface hydraulic conductivity accordingly with the prescribed method.")
	InterfaceConductivity interfaceThermalConductivity;

	double cpAir = 1003.4; 
	double epsilon = 0.0;
	double epsilonAir = 0.0;
	double densityAir = 1.2;
	double densityWater = 1000;
	double cwWater = 4188;
	double ra = 1.0;

	double[] uDiagonal;
	double[] mDiagonal;
	double[] lDiagonal;
	double[] r;
	double[] tmp;
	///////////////////////////////



	public Energy1DAdvection(String soilHeatCapacityModel, String soilThermalConductivityModel, String topBCType,
			String bottomBCType, String interfaceHydraulicCondType, int NUM_CONTROL_VOLUMES, int picardIteration,
			double tTimestep, double[] deltaZ, double[] spaceDeltaZ, double[] temperatureIC) {

		soilHeatCapacityFactory = new SoilHeatCapacityFactory();
		soilHeatCapacity = soilHeatCapacityFactory.createSoilHeatCapacity(soilHeatCapacityModel);

		soilThermalConductivityFactory = new SimpleSoilThermalConductivityFactory();
		soilThermalConductivity = soilThermalConductivityFactory.createThermalSoilParameterization(soilThermalConductivityModel);

		totalDepth = new TotalDepth();

		SimpleBoundaryConditionFactory boundCondFactory = new SimpleBoundaryConditionFactory();
		//this.psiIC = psiIC;
		this.topBCType = topBCType;
		this.bottomBCType = bottomBCType;
		topBoundaryCondition = boundCondFactory.createBoundaryCondition(topBCType);		
		bottomBoundaryCondition = boundCondFactory.createBoundaryCondition(bottomBCType);	

		SimpleInterfaceConductivityFactory interfaceConductivityFactory = new SimpleInterfaceConductivityFactory();
		interfaceThermalConductivity = interfaceConductivityFactory.createInterfaceConductivity(interfaceHydraulicCondType);


		dx			  = new double[NUM_CONTROL_VOLUMES];
		this.spaceDelta = spaceDeltaZ;
		for(int i = 0; i < NUM_CONTROL_VOLUMES-1; i++) {
			dx[i] = deltaZ[i];
		}
		this.dx[NUM_CONTROL_VOLUMES-1] = deltaZ[deltaZ.length-1];
		this.picardIteration = picardIteration;
		this.tTimestep = tTimestep;
		this.NUM_CONTROL_VOLUMES = NUM_CONTROL_VOLUMES;

		compute = new ComputeDerivedQuantities(NUM_CONTROL_VOLUMES, dx, spaceDelta, soilHeatCapacity, soilThermalConductivity , totalDepth, interfaceThermalConductivity, bottomBCType);

		thomasAlg = new Thomas();

		// conversion from degree to radiant of slope angle
		delta = delta*Math.PI/180;


		lowerDiagonal = new double[NUM_CONTROL_VOLUMES];
		mainDiagonal  = new double[NUM_CONTROL_VOLUMES];
		upperDiagonal = new double[NUM_CONTROL_VOLUMES];
		rhss 		  = new double[NUM_CONTROL_VOLUMES];
		uDiagonal = new double[NUM_CONTROL_VOLUMES-1];
		mDiagonal = new double[NUM_CONTROL_VOLUMES-1];
		lDiagonal = new double[NUM_CONTROL_VOLUMES-1];
		r = new double[NUM_CONTROL_VOLUMES-1];
		tmp = new double[NUM_CONTROL_VOLUMES-1];
		lambdaP       = 0.0;
		lambdaM	      = 0.0;

	}


	/**
	 * 
	 * @param topBC
	 * @param bottomBC
	 * @param inCurrentDate
	 * @param timeDelta
	 */
	public void solve(double topBC, double bottomBC, double rain, double  shortWave, double wind,  String inCurrentDate, double timeDelta) {

		System.out.println("ENERGY 1D ADVECTION"+inCurrentDate);


		this.timeDelta = timeDelta;

		/* COMPUTE HEAT CAPACITIES AT TIME LEVEL n and n+1
		 * THERMAL CONDUCTIVITY AT TIME LEVEL n and 
		 * INTERNAL ENERGY AT TIME LEVEL n
		 */
		compute.setComputeDerivedQuantities(bottomBC);
		compute.computeHeatCapacitiesOld();
		compute.computeHeatCapacities();
		//compute.computeLambdas();
		compute.computeInternalEnergyOld();

		internalEnergy = compute.computeTotalInternalEnergyOld();


		/* COEFFICIENT MATRIX IS BUILD BY THREE VECTORS COLLECTING ELEMENTS OF THE THREE DIAGONAL:
				   	 a lower diagonal T_(i+1)
				   	 b main diagonal  T_i
				   	 c upper diagonal T_(i-1)
				   	 RIGHT HAND SIDE */
		if(Variables.thetas[NUM_CONTROL_VOLUMES-1]>0.0) {
			for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
				if( i == 0 ) {

					//						lambdaP = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i+1],dx[i],dx[i+1]);
					//						if (bottomBCType.equalsIgnoreCase("Bottom Neumann") || bottomBCType.equalsIgnoreCase("BottomNeumann")) {
					//							
					//							lambdaM = -999;
					//						
					//						} else {
					//
					//							lambdaM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i],dx[i],dx[i]);
					//						}

					lowerDiagonal[i] =  0;

					mainDiagonal[i] = Variables.heatCapacities[i]*dx[i]	+ timeDelta*densityWater*cwWater*( 0.5*Variables.darcyVelocities[i+1] + 0.5*Math.abs(Variables.darcyVelocities[i+1]) - (0.5*Variables.darcyVelocities[i] - 0.5*Math.abs(Variables.darcyVelocities[i])) );

					upperDiagonal[i] = + timeDelta*densityWater*cwWater*( 0.5*Variables.darcyVelocities[i+1] - 0.5*Math.abs(Variables.darcyVelocities[i+1]) );

					rhss[i] = Variables.internalEnergiesOld[i] + Variables.heatCapacities[i]*dx[i]*273.15 - timeDelta*densityWater*cwWater*( -0.5*Variables.darcyVelocities[i]*bottomBC - 0.5*Math.abs(Variables.darcyVelocities[i])*bottomBC )
							+ timeDelta*densityWater*cwWater*( Variables.darcyVelocities[i+1]*273.15 - Variables.darcyVelocities[i]*273.15 );

				} else if(i == NUM_CONTROL_VOLUMES -1) {
					/* FIXME: sistemare le classi per le boundary al momento solo Dirichlet per il flusso diffusivo meno,
					 * mancano i flussi alla superficie (quelli commentati)
					 * 
					 */
					//lambdaP = 0;	
					//lambdaM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i-1],dx[i],dx[i-1]);
					lowerDiagonal[i] = + timeDelta*densityWater*cwWater*( -0.5*Variables.darcyVelocities[i] - 0.5*Math.abs(Variables.darcyVelocities[i]) ); 

					mainDiagonal[i] = Variables.heatCapacities[i] + timeDelta*densityWater*cwWater*( +0.5*(-rain) + 0.5*Math.abs(-rain) -0.5*Variables.darcyVelocities[i] + 0.5*Math.abs(Variables.darcyVelocities[i]) );

					upperDiagonal[i] = 0;//topBoundaryCondition.upperDiagonal(-999, lambdaP, lambdaM,  -999, spaceDelta[i], timeDelta, delta);

					rhss[i] = Variables.internalEnergiesOld[i] + Variables.heatCapacities[i]*273.15
							- timeDelta*densityWater*cwWater*( 0.5*(-rain)*(topBC) - 0.5*Math.abs((-rain))*(topBC) )
							+ timeDelta*densityWater*cwWater*( -rain*273.15 - Variables.darcyVelocities[i]*273.15 ) + timeDelta*shortWave;;
					//flussi di energia alla superficie
					// + timeDelta*shortWave + timeDelta*epsilon*epsilonAir*5.67*Math.pow(10, -8)*Math.pow(topBC, 4) - timeDelta*epsilon*5.67*Math.pow(10, -8)*Math.pow(Variables.temperatures[i], 4)
					//+ timeDelta*densityAir*cpAir*wind/ra*topBC;

				} else {
					//lambdaP = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i+1],dx[i],dx[i+1]);
					//lambdaM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i-1],dx[i],dx[i-1]);

					lowerDiagonal[i] = + timeDelta*densityWater*cwWater*( -0.5*Variables.darcyVelocities[i] - 0.5*Math.abs(Variables.darcyVelocities[i]) );

					mainDiagonal[i] = Variables.heatCapacities[i]*dx[i] + timeDelta*densityWater*cwWater*( 0.5*Variables.darcyVelocities[i+1] + 0.5*Math.abs(Variables.darcyVelocities[i+1]) - 0.5*Variables.darcyVelocities[i] + 0.5*Math.abs(Variables.darcyVelocities[i]) );

					upperDiagonal[i] = + timeDelta*densityWater*cwWater*( 0.5*Variables.darcyVelocities[i+1] - 0.5*Math.abs(Variables.darcyVelocities[i+1]) );

					rhss[i] = Variables.internalEnergiesOld[i] + Variables.heatCapacities[i]*dx[i]*273.15 + timeDelta*densityWater*cwWater*( Variables.darcyVelocities[i+1]*273.15 - Variables.darcyVelocities[i]*273.15 ); ; 

				}

			}

			//				System.out.println("");

			//// THOMAS ALGORITHM ////
			thomasAlg.set(upperDiagonal,mainDiagonal,lowerDiagonal,rhss);
			Variables.temperatures = thomasAlg.solver();
		} else {
			for(int i = 0; i < NUM_CONTROL_VOLUMES-1; i++) {
				if( i == 0 ) {

//					lambdaP = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i+1],dx[i],dx[i+1]);
//					if (bottomBCType.equalsIgnoreCase("Bottom Neumann") || bottomBCType.equalsIgnoreCase("BottomNeumann")) {
//
//						lambdaM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i],dx[i],dx[i]);
//
//					} else {
//
//						lambdaM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i],dx[i],dx[i]);
//					}

					lDiagonal[i] = 0.0;// bottomBoundaryCondition.lowerDiagonal(-999, lambdaP, lambdaM, spaceDelta[i+1], spaceDelta[i], timeDelta, delta);

					mDiagonal[i] = Variables.heatCapacities[i]*dx[i] + timeDelta*densityWater*cwWater*( 0.5*Variables.darcyVelocities[i+1] + 0.5*Math.abs(Variables.darcyVelocities[i+1]) - (0.5*Variables.darcyVelocities[i] - 0.5*Math.abs(Variables.darcyVelocities[i])) );

					uDiagonal[i] = + timeDelta*densityWater*cwWater*( 0.5*Variables.darcyVelocities[i+1] - 0.5*Math.abs(Variables.darcyVelocities[i+1]) );

					r[i] = Variables.internalEnergiesOld[i] +  Variables.heatCapacities[i]*dx[i]*273.15 - timeDelta*densityWater*cwWater*( -0.5*Variables.darcyVelocities[i]*(bottomBC) - 0.5*Math.abs(Variables.darcyVelocities[i])*(bottomBC) )
							+ timeDelta*densityWater*cwWater*( Variables.darcyVelocities[i+1]*273.15 - Variables.darcyVelocities[i]*273.15 );


					//FIXME nel caso di neumann
					//- timeDelta*densityWater*cwWater*( -0.5*Variables.darcyVelocities[i]*Variables.temperatures[i] - 0.5*Math.abs(Variables.darcyVelocities[i])*Variables.temperatures[i] );

				} else if(i == NUM_CONTROL_VOLUMES -2) {
					/* FIXME: sistemare le classi per le boundary al momento solo Dirichlet per il flusso diffusivo meno,
					 * mancano i flussi alla superficie (quelli commentati)
					 * 
					 */
//					lambdaP = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i],dx[i],dx[i]);	
//					lambdaM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i-1],dx[i],dx[i-1]);
					lDiagonal[i] = timeDelta*densityWater*cwWater*( -0.5*Variables.darcyVelocities[i] - 0.5*Math.abs(Variables.darcyVelocities[i]) );

					mDiagonal[i] = Variables.heatCapacities[i]*dx[i] + timeDelta*densityWater*cwWater*( +0.5*(Variables.darcyVelocities[i+1]) + 0.5*Math.abs(Variables.darcyVelocities[i+1]) - 0.5*Variables.darcyVelocities[i] + 0.5*Math.abs(Variables.darcyVelocities[i]) );

					uDiagonal[i] = 0.0;//topBoundaryCondition.upperDiagonal(topBC, lambdaP, lambdaM,  -999, spaceDelta[i], timeDelta, delta);

					r[i] = Variables.internalEnergiesOld[i] + Variables.heatCapacities[i]*dx[i]*273.15 - timeDelta*densityWater*cwWater*( 0.5*(Variables.darcyVelocities[i+1])*(topBC) - 0.5*Math.abs((Variables.darcyVelocities[i+1]))*(topBC) )
							+ timeDelta*densityWater*cwWater*( Variables.darcyVelocities[i+1]*273.15 - Variables.darcyVelocities[i]*273.15 );
					//flussi di energia alla superficie
					//+ timeDelta*shortWave + timeDelta*epsilon*epsilonAir*5.67*Math.pow(10, -8)*Math.pow(topBC, 4)  + timeDelta*epsilon*5.67*Math.pow(10, -8)*Math.pow(Variables.temperatures[i], 4)
					//+ timeDelta*densityAir*cpAir*wind/ra*topBC;

				} else {
//					lambdaP = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i+1],dx[i],dx[i+1]);
//					lambdaM = interfaceThermalConductivity.compute(Variables.lambdas[i],Variables.lambdas[i-1],dx[i],dx[i-1]);

					lDiagonal[i] = + timeDelta*densityWater*cwWater*( -0.5*Variables.darcyVelocities[i] - 0.5*Math.abs(Variables.darcyVelocities[i]) );

					mDiagonal[i] = Variables.heatCapacities[i]*dx[i] + timeDelta*densityWater*cwWater*( 0.5*Variables.darcyVelocities[i+1] + 0.5*Math.abs(Variables.darcyVelocities[i+1]) - 0.5*Variables.darcyVelocities[i] + 0.5*Math.abs(Variables.darcyVelocities[i]) );

					uDiagonal[i] = timeDelta*densityWater*cwWater*( 0.5*Variables.darcyVelocities[i+1] - 0.5*Math.abs(Variables.darcyVelocities[i+1]) );

					r[i] = Variables.internalEnergiesOld[i] + Variables.heatCapacities[i]*dx[i]*273.15 + timeDelta*densityWater*cwWater*( Variables.darcyVelocities[i+1]*273.15 - Variables.darcyVelocities[i]*273.15 ); 

				}

			}

			//				System.out.println("");


			//// THOMAS ALGORITHM ////
			thomasAlg.set(uDiagonal,mDiagonal,lDiagonal,r);
			tmp = thomasAlg.solver();
			for(int i = 0; i < NUM_CONTROL_VOLUMES-1; i++) {
				Variables.temperatures[i] = tmp[i];
			}
			Variables.temperatures[NUM_CONTROL_VOLUMES-1] = tmp[NUM_CONTROL_VOLUMES-2];
		}
		//				//// THOMAS ALGORITHM ////
		//					if(rhss[NUM_CONTROL_VOLUMES -1] == 0) {
		//						double tmp = Variables.temperatures[NUM_CONTROL_VOLUMES -1];
		//						for(int i=0; i<NUM_CONTROL_VOLUMES -1;i++) {
		//							uD[i] = upperDiagonal[i];
		//							mD[i] = mainDiagonal[i];
		//							lD[i] = lowerDiagonal[i];
		//							r[i] = rhss[i];
		//						}
		//						thomasAlg.set(uD,mD,lD,r);
		//						t = thomasAlg.solver();
		//						for(int i=0; i<NUM_CONTROL_VOLUMES-1;i++) {
		//							Variables.temperatures[i] = t[i];
		//						}
		//						Variables.temperatures[NUM_CONTROL_VOLUMES -1] = tmp;
		//					} else {
		//						thomasAlg.set(upperDiagonal,mainDiagonal,lowerDiagonal,rhss);
		//						Variables.temperatures = thomasAlg.solver();
		//					}



		/* COMPUTE DERIVED QUANTITIES
		 * INTERNAL ENERGIES at time level n+1
		 * 
		 * FIXME: mancano i flussi di calore e il numero di Peclet 
		 */ 
		compute.computeInternalEnergyNew();
		compute.computeEnergyFluxes();
		//compute.computeDiffusionFluxes();
		compute.computeAdvectionFluxes();
		//compute.computePeclet();
		internalEnergyNew = compute.computeTotalInternalEnergy();

		//FIXME: dovrebbe essere OK
		if(Variables.thetas[NUM_CONTROL_VOLUMES-1]>0.0) {
			Variables.errorEnergy = internalEnergyNew - internalEnergy - timeDelta*( -densityWater*cwWater*( 0.5*(-rain)*(topBC-273.15 + Variables.temperatures[NUM_CONTROL_VOLUMES -1]-273.15) - 0.5*Math.abs(-rain)*(topBC-273.15 - Variables.temperatures[NUM_CONTROL_VOLUMES -1]+273.15) )
					- densityWater*cwWater*( -0.5*Variables.darcyVelocities[0]*(Variables.temperatures[0]-273.15+bottomBC-273.15) + 0.5*Math.abs(Variables.darcyVelocities[0])*(Variables.temperatures[0]-273.15-bottomBC+273.15) ) );
			//+ timeDelta*densityWater*cwWater*( rain*topBC ) + timeDelta*shortWave + timeDelta*epsilon*epsilonAir*5.67*Math.pow(10, -8)*Math.pow(topBC, 4) - timeDelta*epsilon*5.67*Math.pow(10, -8)*Math.pow(Variables.temperatures[NUM_CONTROL_VOLUMES-1], 4)
			//+ timeDelta*densityAir*cpAir*wind/ra*topBC);
		} else {
			Variables.errorEnergy = internalEnergyNew - (internalEnergy-Variables.internalEnergiesOld[NUM_CONTROL_VOLUMES-1]) - timeDelta*( -densityWater*cwWater*( 0.5*(Variables.darcyVelocities[NUM_CONTROL_VOLUMES -1])*(topBC-273.15 + Variables.temperatures[NUM_CONTROL_VOLUMES -2]-273.15) - 0.5*Math.abs(Variables.darcyVelocities[NUM_CONTROL_VOLUMES -1])*(topBC-273.15 - Variables.temperatures[NUM_CONTROL_VOLUMES -2]+273.15) )
					- densityWater*cwWater*( -0.5*Variables.darcyVelocities[0]*(Variables.temperatures[0]-273.15+bottomBC-273.15) + 0.5*Math.abs(Variables.darcyVelocities[0])*(Variables.temperatures[0]-273.15-bottomBC+273.15) ) );
		}
		System.out.println();
	} //// MAIN CYCLE END ////


}  /// CLOSE Richards1d ///



