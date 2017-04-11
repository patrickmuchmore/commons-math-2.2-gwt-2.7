package org.feno.estimate.client.plot;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import com.googlecode.gwt.charts.client.DataSource;
import com.googlecode.gwt.charts.client.corechart.ColumnChart;

public class Histogram extends ColumnChart {

	private final DataSource dataSource;
	private final SummaryStatistics sampleStats;
	private static final SummaryStatistics getSampleStats(double[] values) {
		SummaryStatistics sampleStats = new SummaryStatistics();
		for (int i = 0; i < values.length; i++) {
			sampleStats.addValue(values[i]);
		}
		return sampleStats;
	}

	// Prevent direct instantiation.
	@SuppressWarnings("unused")
	private Histogram() {
		this(new double[] { Double.NaN });
	}

	public Histogram(double[] values) {
		this(values, getSampleStats(values));
	}

	public Histogram(double[] values, SummaryStatistics summaryStatistics) {
		sampleStats = summaryStatistics;
		dataSource = setDataSource(values, sampleStats);

	}

	private DataSource setDataSource(double values[], SummaryStatistics sampleStats) {

		return null;
	}

	public void draw() {
		super.draw(dataSource);
	}

}
