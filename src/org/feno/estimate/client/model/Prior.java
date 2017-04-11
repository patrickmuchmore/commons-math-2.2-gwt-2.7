package org.feno.estimate.client.model;

public interface Prior {

	public double density(double x, boolean log);

	public double getMin();

	public double getMax();

	public class Uniform implements Prior {
		private final double min, max, density, logDensity;

		public Uniform(double min, double max) {
			this.min = min;
			this.max = max;
			this.density = 1 / (this.max - this.min);
			this.logDensity = Math.log(density);
		}

		public double getMin() {
			return min;
		}

		public double getMax() {
			return max;
		}

		@Override
		public double density(double x, boolean log) {
			if (min <= x && x <= max) {
				return log ? logDensity : density;
			} else {
				return log ? Double.NEGATIVE_INFINITY : 0.0;
			}
		}
	}
}
