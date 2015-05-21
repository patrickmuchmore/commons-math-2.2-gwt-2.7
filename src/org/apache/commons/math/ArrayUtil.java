package org.apache.commons.math;

public class ArrayUtil {
  
  public static final double[] cloneDouble(double[] orig) {
    double[] copy = new double[orig.length];
    System.arraycopy(orig, 0, copy, 0, orig.length);
    return copy;
  }

  public static final Object[] cloneObject(Object[] orig) {
    Object[] copy = new Object[orig.length];
    System.arraycopy(orig, 0, copy, 0, orig.length);
    return copy;
  }
  
}
