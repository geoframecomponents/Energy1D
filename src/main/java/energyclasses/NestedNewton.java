package energyclasses;

import physicalquantities.Variables;
import soilparameters.SoilParameters;
import swrc.SoilWaterRetentionCurve;

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


/**
 * This class carries out the Nested-Newton algorithm
 * (A NESTED NEWTON-TYPE ALGORITHM FOR FINITE VOLUME METHODS SOLVING RICHARDS' EQUATION IN MIXED FORM, Casulli V., Zanolli P., Journal Scientific Computing, 2010)
 *  @author Niccolo' Tubini
 */

public class NestedNewton {
	private double outerResidual;
	private double innerResidual;
	
	int nestedNewton;
	int MAXITER_NEWT;
	int NUM_CONTROL_VOLUMES;
	
	double newtonTolerance;
	
	double[] psis;
	double[] mainDiagonal;
	double[] upperDiagonal;
	double[] lowerDiagonal;
	double[] rhss;
	double[] dx;
	
	//double[] thetaS;
	//double[] thetaR;
	//double[] par1SWRC;
	//double[] par2SWRC;
	//double[] par3SWRC;
	//double[] par4SWRC;
	//double[] par5SWRC;
	//double[] alphaSpecificStorage;
	//double[] betaSpecificStorage;
	
	double[] fs;
	double[] fks;
	double[] bb;
	double[] cc;
	double[] dis;
	double[] dpsis;
	double[] psis_outer;
	
	SoilWaterRetentionCurve soilWaterRetentionCurve;
	TotalDepth totalDepth;
	Thomas thomasAlg = new Thomas();

	
	
	/**
	 * @param nestedNewton control parameter to choose between simple Newton method (0), or the nested Newton one (1)
	 * @param newtonTolerance prefixed tolerance representing the maximum mass balance error allowed  
	 * @param MAXITER_NEWT prefixed maximum number of iteration
	 * @param NUM_CONTROL_VOLUMES number of control volumes
	 * @param soilPar is the class to compute the soil hydraulic properties
	 * @param totalDepth is the class to compute the total water depth
	 * @param par1SWRC vector containing the first parameter of the SWRC, it is a vector of length NUM_CONTROL_VOLUMES-1
	 * @param par2SWRC vector containing the second parameter of the SWRC, it is a vector of length NUM_CONTROL_VOLUMES-1
	 * @param thetaR vector containing the adimensional residual water contentfor each control volume, it is a vector of length NUM_CONTROL_VOLUMES-1
	 * @param thetaS vector containing the adimensional water content at saturation for each control volume, it is a vector of length NUM_CONTROL_VOLUMES-1
	 */
	public NestedNewton(int nestedNewton, double newtonTolerance, int MAXITER_NEWT, int NUM_CONTROL_VOLUMES, double[] dx, SoilWaterRetentionCurve soilWaterRetentionCurve, TotalDepth totalDepth){
		
		this.nestedNewton = nestedNewton;
		this.newtonTolerance = newtonTolerance;
		this.MAXITER_NEWT = MAXITER_NEWT;
		this.NUM_CONTROL_VOLUMES = NUM_CONTROL_VOLUMES;
		this.soilWaterRetentionCurve = soilWaterRetentionCurve;
		this.totalDepth = totalDepth;
		this.dx = dx;
		
		psis          = new double[this.NUM_CONTROL_VOLUMES];
		fs			  = new double[this.NUM_CONTROL_VOLUMES];
		fks			  = new double[this.NUM_CONTROL_VOLUMES];
		bb            = new double[this.NUM_CONTROL_VOLUMES]; 
		cc 			  = new double[this.NUM_CONTROL_VOLUMES];
		dis			  = new double[this.NUM_CONTROL_VOLUMES];
		dpsis		  = new double[this.NUM_CONTROL_VOLUMES];
		psis_outer	  = new double[this.NUM_CONTROL_VOLUMES];
	}
	
	
	
	/**
	 * @param psis vector contains the suction values, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param upperDiagonal upper diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param mainDiagonal main diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param lowerDiagonal lower diagonal of the coefficient matrix A of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 * @param rhss right hand side term of the linear system, it is a vector of length NUM_CONTROL_VOLUMES
	 */
	public void set(double[] mainDiagonal, double[] upperDiagonal, double[] lowerDiagonal, double[] rhss){
		
		//this.psis = psis;
		this.mainDiagonal = mainDiagonal;
		this.upperDiagonal = upperDiagonal;
		this.lowerDiagonal = lowerDiagonal;
		this.rhss = rhss;

	}

	
	
	public void solver(){

		// Initial guess of psis
		for(int i = 0; i < NUM_CONTROL_VOLUMES; i++) {
			if(i==NUM_CONTROL_VOLUMES-1) {
				Variables.waterSuctions[i] = Math.max(Variables.waterSuctions[i],0.1);
				//System.out.println(i +"   "+psis[i]);
			} else {
				Variables.waterSuctions[i] = Math.min(Variables.waterSuctions[i], SoilParameters.psiStar1[i] );
				//System.out.println(i +"   "+psis[i]);
			}
		}

		//// OUTER CYCLE ////
		for(int i = 0; i < MAXITER_NEWT; i++) {
			// I have to assign 0 to outerResidual otherwise I will take into account of the previous error
			outerResidual = 0.0;
			for(int j = 0; j < NUM_CONTROL_VOLUMES; j++) {
				if(j==0) {
					fs[j] = soilWaterRetentionCurve.waterContent(j)*dx[j] - rhss[j] + mainDiagonal[j]*Variables.waterSuctions[j] + upperDiagonal[j]*Variables.waterSuctions[j+1];
					dis[j] = soilWaterRetentionCurve.dWaterContent(j);
					//System.out.println(j+" "+fs[j]);
				} else if(j==NUM_CONTROL_VOLUMES-1) {
					fs[j] = totalDepth.totalDepth(Variables.waterSuctions[j]) - rhss[j] + lowerDiagonal[j]*Variables.waterSuctions[j-1] + mainDiagonal[j]*Variables.waterSuctions[j];
					dis[j] = totalDepth.dTotalDepth(Variables.waterSuctions[j]);
					//System.out.println(j+" "+fs[j]);
				} else {
					fs[j] = soilWaterRetentionCurve.waterContent(j)*dx[j] - rhss[j] + lowerDiagonal[j]*Variables.waterSuctions[j-1] + mainDiagonal[j]*Variables.waterSuctions[j] + upperDiagonal[j]*Variables.waterSuctions[j+1];
					dis[j] = soilWaterRetentionCurve.dWaterContent(j);
					//System.out.println(j+" "+soilPar.waterContent(psis[j],j));
					//System.out.println(j+" "+fs[j]);
				}
			
				outerResidual += fs[j]*fs[j];
			}
			outerResidual = Math.pow(outerResidual,0.5);  
			//System.out.println("   Outer iteration " + i + " with residual " +  outerResidual);
			if(outerResidual < newtonTolerance) {
				break;
			}
			if(nestedNewton == 0){
				bb = mainDiagonal.clone();
				cc = upperDiagonal.clone();
				for(int y = 0; y < NUM_CONTROL_VOLUMES; y++) {
					bb[y] += dis[y];
				}
				thomasAlg.set(cc,bb,lowerDiagonal,fs);
				dpsis = thomasAlg.solver();

				//// PSIS UPDATE ////
				for(int s = 0; s < NUM_CONTROL_VOLUMES; s++) {
					Variables.waterSuctions[s] = Variables.waterSuctions[s] - dpsis[s];
				}
			}else{
				psis_outer = Variables.waterSuctions.clone();

				// Initial guess of psis
				for(int j = 0; j < NUM_CONTROL_VOLUMES; j++) {
					if(j==NUM_CONTROL_VOLUMES-1) {
						Variables.waterSuctions[j] = Math.max(Variables.waterSuctions[j],0.11);
					} else {
						Variables.waterSuctions[j] = Math.max(Variables.waterSuctions[j], SoilParameters.psiStar1[j]);
					}
				}

				//// INNER CYCLE ////
				for(int j = 0; j < MAXITER_NEWT; j++) {
					// I have to assign 0 to innerResidual otherwise I will take into account of the previous error
					innerResidual = 0.0; 
					for(int l=0; l < NUM_CONTROL_VOLUMES; l++) {
						if(l==0) {
							fks[l] = soilWaterRetentionCurve.pIntegral(l)*dx[l] - ( soilWaterRetentionCurve.qIntegral(psis_outer[l],l) + soilWaterRetentionCurve.q(psis_outer[l],l)*(Variables.waterSuctions[l] - psis_outer[l]) )*dx[l] - this.rhss[l] + mainDiagonal[l]*Variables.waterSuctions[l] + upperDiagonal[l]*Variables.waterSuctions[l+1];
							dis[l] = ( soilWaterRetentionCurve.p(l) - soilWaterRetentionCurve.q(psis_outer[l],l) )*dx[l];
							//System.out.println(l+" "+fks[l]);
							//System.out.println(l+" "+dis[l]);
						} else if(l==NUM_CONTROL_VOLUMES-1) {
							//fks[l] = totalDepth.pIntegral(psis[l]) - ( totalDepth.qIntegral(psis_outer[l]) + totalDepth.q(psis_outer[l])*(psis[l] - psis_outer[l]) ) - this.rhss[l] + lowerDiagonal[l]*psis[l-1] + mainDiagonal[l]*psis[l];
							//dis[l] = totalDepth.p(psis[l]) - totalDepth.q(psis_outer[l]);
							fks[l] = totalDepth.pIntegral(Variables.waterSuctions[l])- this.rhss[l] + lowerDiagonal[l]*Variables.waterSuctions[l-1] + mainDiagonal[l]*Variables.waterSuctions[l];
							dis[l] = totalDepth.p(Variables.waterSuctions[l]);
							//System.out.println(l+" "+fks[l]);
							//System.out.println(l+" "+totalDepthJordanDecomposition.p(psis[l]));
						} else {
							fks[l] = soilWaterRetentionCurve.pIntegral(l)*dx[l] - ( soilWaterRetentionCurve.qIntegral(psis_outer[l],l) + soilWaterRetentionCurve.q(psis_outer[l],l)*(Variables.waterSuctions[l] - psis_outer[l]) )*dx[l] - this.rhss[l]  + lowerDiagonal[l]*Variables.waterSuctions[l-1] + mainDiagonal[l]*Variables.waterSuctions[l] + upperDiagonal[l]*Variables.waterSuctions[l+1];
							dis[l] = ( soilWaterRetentionCurve.p(l) - soilWaterRetentionCurve.q(psis_outer[l],l) )*dx[l];
							//System.out.println(l+" "+fks[l]);
							//System.out.println(l+" "+dis[l]);
						}

						innerResidual += fks[l]*fks[l];
					}
					innerResidual = Math.pow(innerResidual,0.5);

					//System.out.println("     -Inner iteration " + j + " with residual " +  innerResidual);    

					if(innerResidual < newtonTolerance) {
						break;
					}

					//// THOMAS ALGORITHM////
					// Attention: the main diagonal of the coefficient matrix must not change!! The same for the upper diagonal

					bb = mainDiagonal.clone();
					cc = upperDiagonal.clone();
					for(int y = 0; y < NUM_CONTROL_VOLUMES; y++) {
						bb[y] += dis[y];
					}
					thomasAlg.set(cc,bb,lowerDiagonal,fks);
					dpsis = thomasAlg.solver();

					//// PSIS UPDATE ////
					for(int s = 0; s < NUM_CONTROL_VOLUMES; s++) {
						Variables.waterSuctions[s] = Variables.waterSuctions[s] - dpsis[s];
					}
				}
			} //// INNER CYCLE END ////
		}
		//return psis;
	}
	
}
