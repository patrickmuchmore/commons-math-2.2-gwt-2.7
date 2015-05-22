package org.feno.estimate.client.plot;

import org.apache.commons.lang3.math.Fraction;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.feno.estimate.client.Util;
import org.feno.estimate.client.mcmc.PosteriorSample;

import com.google.gwt.user.client.ui.RequiresResize;
import com.googlecode.gwt.charts.client.DataSource;
import com.googlecode.gwt.charts.client.corechart.ColumnChart;

public class Histogram extends ColumnChart {

  private final DataSource dataSource;
  private final SummaryStatistics sampleStats;
  private double binWidth;
  private int nBins;
  private int[] binCounts;
  private double[] binUpperLimits;
 
  private static final SummaryStatistics getSampleStats(double[] values) {
    SummaryStatistics sampleStats = SummaryStatistics.newInstance();
    for(int i=0; i<values.length; i++) {
      sampleStats.addValue(values[i]);
    }
    return sampleStats;
  }
  //Prevent direct instantiation.
  @SuppressWarnings("unused")
  private Histogram() {
    this(new double[]{Double.NaN});
  } 
  
  public Histogram(double[] values) {
    this(values, getSampleStats(values));
  }
  
  public Histogram(double[] values, SummaryStatistics summaryStatistics) {
    sampleStats = summaryStatistics;
    binWidth = 3.5*sampleStats.getStandardDeviation() * 
        Math.pow(sampleStats.getN(), -Fraction.ONE_THIRD.doubleValue());
    nBins = (int) Math.ceil((sampleStats.getMax()-sampleStats.getMin())/binWidth);
    dataSource = setDataSource(values, sampleStats);
    
  }

//  public Histogram(PosteriorSample sample) {
//    this(sample.getSampleValues()., sample.getSampleStats());
//  }

  private DataSource setDataSource(double values[], SummaryStatistics sampleStats) {


    return null;
  }

  public void draw() {
    super.draw(dataSource);
  }

}
