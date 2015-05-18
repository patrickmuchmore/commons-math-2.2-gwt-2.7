package ca.uol.aig.fftpack;

public class DoubleFFT_1D extends ComplexDoubleFFT {

  public DoubleFFT_1D(int n) {
    super(n);
  }

  public void realForwardFull(double[] new_kords) {
    if(new_kords.length % 2 != 0) {
      throw new IllegalStateException("Array length must be a multiple of 2.");
    }
    int halfLength = new_kords.length/2;
    for(int i=1; i<halfLength; i+=2) {
      new_kords[halfLength+i] = new_kords[i];
      new_kords[i] = 0;
    }
    super.ft(new_kords);
  }

  public void complexInverse(double[] kords, boolean b) {
    super.bt(kords);
//    ComplexDoubleFFT complexFFT = new ComplexDoubleFFT(kords.length);
//    kords = padWithZeros(kords);
//    complexFFT.bt(kords);
  }
  
  private static double[] padWithZeros(double[] values) {
    double[] padded = new double[values.length*2];
    System.arraycopy(values, 0, padded, 0, values.length);
    //System.arraycopy((double[]) new int[values.length], 0, padded, values.length+1, values.length);
    return padded;
  }

}
