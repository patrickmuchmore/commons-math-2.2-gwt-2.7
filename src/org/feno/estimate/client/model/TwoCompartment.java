package org.feno.estimate.client.model;

import java.util.Map;

public class TwoCompartment extends Model {

  public enum Parm {
    CANO,
    JAWNO,
    DAWNO,
  }
  
  public TwoCompartment() {
    super(Parm.CANO.name(), Parm.JAWNO.name(), Parm.DAWNO.name());
  }

  private Double getParmVal(Parm parm) {
    return super.getParameterValue(parm.name());
  }
  
  @Override
  public double fenoEstimate(double vdot) {
    //twoCompartmentFeNO = function(JawNO, DawNO, CANO, Vdot) {
    //  return(JawNO/DawNO + (CANO-JawNO/DawNO)*exp(-DawNO/Vdot))
    //}
    double jawOverDaw = getParmVal(Parm.JAWNO)/getParmVal(Parm.DAWNO);
    return jawOverDaw + (getParmVal(Parm.CANO) - jawOverDaw) * Math.exp(-getParmVal(Parm.DAWNO)/vdot);
  }

  @Override
  public double logLikelihood() {
    // TODO Auto-generated method stub
    return 0;
  }
  
}
