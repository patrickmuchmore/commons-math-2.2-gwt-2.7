package org.feno.estimate.client.mcmc;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import com.google.gwt.typedarrays.client.Float64ArrayNative;

public class PosteriorSample {

	private final int sampleSize;
	private final Float64ArrayNative sampleValues;
	private final SummaryStatistics sampleStats;

	public PosteriorSample(int size) {
		sampleSize = size;
		sampleValues = Float64ArrayNative.create(size);
		sampleStats = new SummaryStatistics();
	}

	public int getSampleSize() {
		return sampleSize;
	}

	public SummaryStatistics getSampleStats() {
		return sampleStats;
	}

	public Float64ArrayNative getSampleValues() {
		return sampleValues;
	}

	public void setSampleValue(int index, double value) {
		sampleValues.set(index, value);
		sampleStats.addValue(value);
	}

}
