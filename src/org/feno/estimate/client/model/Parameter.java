package org.feno.estimate.client.model;

public class Parameter {

  private final String name;
  private final Prior priorDistribution;
  private double value;
  private final double propSd;
  
  //public Parameter(String parameterName, GenericDistribution priorDist) {
  //  this(parameterName, priorDist, Double.NaN);
  //}
  
  public Parameter(
      String parameterName, 
      Prior priorDist,
      double initialValue,
      double proposalSd) {
    
    name = parameterName;
    priorDistribution = priorDist;
    value = initialValue;
    propSd = proposalSd;
    
  }
  
  public double getValue() {
    return value;
  }
  
  public void setValue(double newValue) {
    value = newValue;
  }
  
  public String getName() {
    return name;
  }
  
  public double getPropSd() {
    return propSd;
  }
  
  public double priorDensity(double x, boolean log) {
    return priorDistribution.density(x, log);
  }
  
  public double reflectProposal(double proposedValue) {
    while(proposedValue < priorDistribution.getMin() || proposedValue > priorDistribution.getMax()) {
      if(proposedValue < priorDistribution.getMin()) {
        proposedValue = 2*priorDistribution.getMin() - proposedValue;
      } else {
        proposedValue = 2*priorDistribution.getMax() - proposedValue;
      }
    }
    return proposedValue;
  }
  
}
