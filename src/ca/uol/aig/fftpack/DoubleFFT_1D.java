package ca.uol.aig.fftpack;

public class DoubleFFT_1D extends RealDoubleFFT {

  public DoubleFFT_1D(int n) {
    super(n);
  }

  public void realForwardFull(double[] new_kords) {
    super.ft(new_kords);
  }

  public void complexInverse(double[] kords, boolean b) {
    if(kords.length % 2 != 0) {
      throw new IllegalStateException("Array argument does not have an even number of entries");
    }
    ComplexDoubleFFT complexFFT = new ComplexDoubleFFT(kords.length/2);
    complexFFT.bt(kords);
  }

}
