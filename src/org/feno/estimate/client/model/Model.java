package org.feno.estimate.client.model;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;

import org.feno.estimate.client.mcmc.PosteriorSample;

import jdistlib.Normal;
import jdistlib.rng.RandomEngine;

import com.google.gwt.typedarrays.client.Float64ArrayNative;

public abstract class Model {
  
  private final HashMap<String, Parameter> parameters = new HashMap<String, Parameter>();

  public Model(Iterable<Parameter> params) {
    for(Parameter param : params) {
      parameters.put(param.getName(), param);
    }
    if(parameters.size() == 0) {
      throw new IllegalArgumentException("A Model must have at least one Parameter.");
    }
  }
  
  public Model(Parameter... params) {
    this(Arrays.asList(params));
  }
  
  public Parameter getParameter(String name) {
    return parameters.get(name);
  }
  
  public Collection<Parameter> getParameters() {
    return parameters.values();
  }
  
  public double getParameterValue(String name) {
    return parameters.get(name).getValue();
  }
  
  public void setParameterValue(String name, double value) {
    parameters.get(name).setValue(value);
  }

  public abstract double fenoEstimate(double vdot);
  
  public abstract double logLikelihood();
  
  public HashMap<String, PosteriorSample> samplePosteriorToo(int steps, RandomEngine rng) {
    
    HashMap<String, PosteriorSample> parameterSamples = new HashMap<String, PosteriorSample>(parameters.size());
   
    for(Parameter param : parameters.values()) {
      parameterSamples.put(param.getName(), new PosteriorSample(steps));
    }
    
    HashMap<String, Double> currentParameterValues = new HashMap<String, Double>(parameters.size());
    double currentLogLikelihood = logLikelihood();
    double proposedLogLikelihood = Double.NEGATIVE_INFINITY;
    
    for(int i=0; i<steps; i++) {
      
      for(Parameter param : parameters.values()) {
        currentParameterValues.put(param.getName(), param.getValue());
        param.setValue(param.reflectProposal(Normal.random(param.getValue(), param.getPropSd(), rng)));
      }
      
      proposedLogLikelihood = logLikelihood();
      
      if(Math.log(rng.nextDouble()) < proposedLogLikelihood - currentLogLikelihood) {
        currentLogLikelihood = proposedLogLikelihood;
      } else {
        for(Parameter param : parameters.values()) {
          param.setValue(currentParameterValues.get(param.getName()));
        }
      }
      
      for(Parameter param : parameters.values()) {
        parameterSamples.get(param.getName()).setSampleValue(i, param.getValue());
      }
    }
    
    return parameterSamples;
  }
 
  public HashMap<String, Float64ArrayNative> samplePosterior(int steps, RandomEngine rng) {
    
      HashMap<String, Float64ArrayNative> parameterSamples = new HashMap<String, Float64ArrayNative>(parameters.size());
     
      for(Parameter param : parameters.values()) {
        parameterSamples.put(param.getName(), Float64ArrayNative.create(steps));
      }
      
      HashMap<String, Double> currentParameterValues = new HashMap<String, Double>(parameters.size());
      double currentLogLikelihood = logLikelihood();
      double proposedLogLikelihood = Double.NEGATIVE_INFINITY;
      
      for(int i=0; i<steps; i++) {
        
        for(Parameter param : parameters.values()) {
          currentParameterValues.put(param.getName(), param.getValue());
          param.setValue(param.reflectProposal(Normal.random(param.getValue(), param.getPropSd(), rng)));
        }
        
        proposedLogLikelihood = logLikelihood();
        
        if(Math.log(rng.nextDouble()) < proposedLogLikelihood - currentLogLikelihood) {
          currentLogLikelihood = proposedLogLikelihood;
        } else {
          for(Parameter param : parameters.values()) {
            param.setValue(currentParameterValues.get(param.getName()));
          }
        }
        
        for(Parameter param : parameters.values()) {
          parameterSamples.get(param.getName()).set(i, param.getValue());
        }
      }
      
      return parameterSamples;
    }
  
}
