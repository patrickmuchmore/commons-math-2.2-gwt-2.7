package org.feno.estimate.client.model;

public interface Prior {

  public double density(double x, boolean log);
  public double getMin();
  public double getMax();
  
  public class Uniform extends jdistlib.Uniform implements Prior {
    private final double min;
    private final double max;

    public Uniform(double min, double max) {
      super(min, max);
      this.min = min;
      this.max = max;
    }

    public double getMin() {
      return min;
    }
    
    public double getMax() {
      return max;
    }
  }
}
