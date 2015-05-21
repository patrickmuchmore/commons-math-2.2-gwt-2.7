package org.feno.estimate.client.model;

import org.feno.estimate.client.mcmc.TestData;

import jdistlib.Normal;
import jdistlib.math.approx.ApproximationFunction;

public class TwoCompartment extends Model {
  
  private static final double[] vDot = {30, 50, 100, 300};
  private static final double[] vDotSd = {3.1, 1.4, 0.8, 0.5};

  public enum Param {
    CANO,
    JAWNO,
    DAWNO,
  }

  //  public TwoCompartment() {
  //    this(Double.NaN, Double.NaN, Double.NaN);
  //  }

  public TwoCompartment(
      Prior canoPrior, double canoInitial, double canoPropSd,
      Prior jawnoPrior, double jawnoInitial, double jawnoPropSd,
      Prior dawnoPrior, double dawnoInitial, double dawnoPropSd) {
    super(
        new Parameter(Param.CANO.name(), canoPrior, canoInitial, canoPropSd), 
        new Parameter(Param.JAWNO.name(), jawnoPrior, jawnoInitial, jawnoPropSd), 
        new Parameter(Param.DAWNO.name(), dawnoPrior, dawnoInitial, dawnoPropSd)
        );
  }

  private Double getParamVal(Param parm) {
    return super.getParameter(parm.name()).getValue();
  }

  @Override
  public double fenoEstimate(double vdot) {
    //twoCompartmentFeNO = function(JawNO, DawNO, CANO, Vdot) {
    //  return(JawNO/DawNO + (CANO-JawNO/DawNO)*exp(-DawNO/Vdot))
    //}
    double jawOverDaw = getParamVal(Param.JAWNO)/getParamVal(Param.DAWNO);
    return jawOverDaw + (getParamVal(Param.CANO) - jawOverDaw) * Math.exp(-getParamVal(Param.DAWNO)/vdot);
  }

  @Override
  public double logLikelihood() {
    double logLikelihood = 0;
    for(int i=0; i<TestData.simulated.length; i++) {
      logLikelihood += Normal.density(
          TestData.simulated[i][1], 
          fenoEstimate(TestData.simulated[i][0]), 
          ApproximationFunction.linear(TestData.simulated[i][0], vDot, vDotSd, 0, Double.POSITIVE_INFINITY), 
          true);
    }
    return logLikelihood;
  }

}