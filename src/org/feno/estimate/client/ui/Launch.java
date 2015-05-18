package org.feno.estimate.client.ui;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Vector;

import org.feno.estimate.client.Util;
import org.feno.estimate.client.model.TwoCompartment;

import jama.Matrix;
import jdistlib.Normal;
import jdistlib.math.density.Density;
import jdistlib.rng.MersenneTwister;
import jdistlib.rng.RandomCMWC;

import com.google.gwt.cell.client.NumberCell;
import com.google.gwt.cell.client.TextCell;
import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.Scheduler.ScheduledCommand;
import com.google.gwt.dom.client.Style.Unit;
import com.google.gwt.user.cellview.client.CellTable;
import com.google.gwt.user.cellview.client.Column;
import com.google.gwt.user.cellview.client.DataGrid;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.*;
import com.googlecode.gwt.charts.client.ChartLoader;
import com.googlecode.gwt.charts.client.ChartPackage;
import com.googlecode.gwt.charts.client.ColumnType;
import com.googlecode.gwt.charts.client.DataTable;
import com.googlecode.gwt.charts.client.corechart.Histogram;
import com.googlecode.gwt.charts.client.corechart.HistogramOptions;
import com.googlecode.gwt.charts.client.corechart.ScatterChart;
import com.googlecode.gwt.charts.client.options.Legend;
import com.googlecode.gwt.charts.client.options.LegendPosition;
import com.googlecode.gwt.charts.client.util.ChartHelper;
import org.feno.estimate.client.model.TwoCompartment.Parm;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class Launch implements EntryPoint {
  private DockLayoutPanel appPanel;
  private Histogram histogram;

  /**
   * This is the entry point method.
   */
  public void onModuleLoad() {
    //Window.enableScrolling(false);
    //Window.setMargin("0px");
    
    TwoCompartment tc = new TwoCompartment();
    HashMap nvmap = tc.getNameValueMap();
    //Util.debugger();
    FlowPanel topPanel = new FlowPanel();
    topPanel.add(new HTML("<h1>Welcome to FeNOnline</h1>"));
    getAppPanel().addNorth(topPanel, 10);
    
//    MenuBar menu = new MenuBar();
//    menu.setHeight("auto");
//    menu.setWidth("auto");
//    MenuBar fileMenu = new MenuBar(true);
//    fileMenu.addItem("Nothing here yet", (ScheduledCommand) null);
//    menu.addItem("File", fileMenu);
//    getAppPanel().addNorth(menu, 2);

    //LayoutPanel rowContent = new LayoutPanel();
    //rowContent.setWidth("10em");
    
    class OptionPair {
      private String name;
      private Number value;
      
      OptionPair(String name, Number value) {
        this.name = name;
        this.value = value;
      }
      
      public String getName() {
        return name;
      }
      public Number getValue() {
        return value;
      }
    }
    
    DataGrid<OptionPair> optionsGrid = new DataGrid<OptionPair>();

    Column<OptionPair, String> nameCol = new Column<OptionPair, String>(new TextCell()) {
      @Override
      public String getValue(OptionPair pair) {
        return pair.getName();
      }
    };
    
    Column<OptionPair, Number> valueCol = new Column<OptionPair, Number>(new NumberCell()) {
      @Override
      public Number getValue(OptionPair pair) {
        return pair.getValue();
      }
    };
    
    Vector<OptionPair> optionPairs = new Vector<OptionPair>();
    OptionPair chainSteps = new OptionPair("Steps", 100);
    optionPairs.add(chainSteps);
    optionsGrid.addColumn(nameCol, "Option Name");
    optionsGrid.addColumn(valueCol, "Current Value");
    optionsGrid.setRowData(optionPairs);
    
    
//    rowContent.add(steps);
//    rowContent.setWidgetLeftWidth(steps, 0, Unit.PCT, 50, Unit.PCT);
//    rowContent.setWidgetTopHeight(steps, 0, Unit.EM, 2, Unit.EM);

//    Label stepsLabel = new Label("Steps");
//    optionsGrid.setWidget(0, 1, stepsLabel);
//    stepsLabel.setSize("75%", "75%");
//    rowContent.add(stepsLabel);
//    rowContent.setWidgetRightWidth(stepsLabel, 0, Unit.PCT, 50, Unit.PCT);
//    rowContent.setWidgetTopHeight(stepsLabel, 0, Unit.EM, 2, Unit.EM);
    
    LayoutPanel leftPanel = new LayoutPanel();
    leftPanel.add(optionsGrid);
    //leftPanel.setWidgetTopHeight(rowContent, 0, Unit.EM, 2, Unit.EM);
    getAppPanel().addWest(optionsGrid, 20);
    RootLayoutPanel.get().add(getAppPanel());
    
    // Create the API Loader
    ChartLoader chartLoader = new ChartLoader(ChartPackage.CORECHART);
    chartLoader.loadApi(new Runnable() {

      @Override
      public void run() {
        getAppPanel().add(getHistogram());
        drawHistogram();
      }
    });
  }
  
  private DockLayoutPanel getAppPanel() {
    if (appPanel == null) {
      appPanel = new DockLayoutPanel(Unit.EM);
    }
    return appPanel;
  }
  
  private Histogram getHistogram() {
    if (histogram == null) {
      histogram = new Histogram();
    }
    return histogram;
  }
  
  private void drawHistogram() {
    //RandomCMWC rng = new RandomCMWC();
    MersenneTwister rng = new MersenneTwister(31415);
    rng.setSeed(31415);
    
    int nDoubles = 1000000;
    
    DataTable randomSample = DataTable.create();
    int valueCol = randomSample.addColumn(ColumnType.NUMBER, "Sample values");
    randomSample.addRows(nDoubles);
    
    for(int i=1; i<nDoubles; i++) {
      randomSample.setCell(i, valueCol, Normal.random_standard(rng));
    }

    // Set options
    HistogramOptions options = HistogramOptions.create();
    options.setTitle("A histogram of " + String.valueOf(nDoubles) + " standard normal sample points");
    options.setLegend(Legend.create(LegendPosition.IN));

    // Draw the chart
    getHistogram().draw(randomSample, options);
  }
  
}
