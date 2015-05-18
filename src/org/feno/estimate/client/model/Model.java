package org.feno.estimate.client.model;

import java.util.Arrays;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

public abstract class Model {
  
  private final HashSet<String> parameterNames;
  private HashMap<String, Double> nameValueMap = new HashMap<String, Double>();

  protected Model(Set<String> names) {
    parameterNames = new HashSet<String>(names);
    Iterator<String> namesIter = names.iterator();
    while(namesIter.hasNext()) {
      nameValueMap.put(namesIter.next(), Double.NaN);
    }
  }
  
  protected Model(String... names) {
    this(new HashSet<String>(Arrays.asList(names)));
  }

  private void checkName(String name) {
    if(!parameterNames.contains(name)) {
      throw new IllegalArgumentException("Unkown parameter key " + name);
    }
  }
   
  public void setParameterValue(String name, Double value) {
    checkName(name);
    nameValueMap.put(name, value);
  }
  
  public Double getParameterValue(String name) {
    checkName(name);
    return nameValueMap.get(name); 
  }
  
  public void setParameterValues(String[] names, Double[] values) {
    if(names.length != values.length) {
      throw new IllegalArgumentException("Array of names and array of values have different lengths.");
    }
    for(int i=0; i<names.length; i++) {
      checkName(names[i]);
      setParameterValue(names[i], values[i]);
    }
  }
  
  public HashMap<String, Double> getNameValueMap() {
    return nameValueMap;
  }

  public abstract double fenoEstimate(double vdot);
  
  public abstract double logLikelihood();
  
}
