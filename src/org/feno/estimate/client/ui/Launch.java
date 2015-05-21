package org.feno.estimate.client.ui;

import java.util.HashSet;
import java.util.Vector;
import java.util.HashMap;

import jdistlib.rng.MersenneTwister;

import org.feno.estimate.client.Util;
import org.feno.estimate.client.mcmc.PosteriorSample;
import org.feno.estimate.client.model.Parameter;
import org.feno.estimate.client.model.Prior.Uniform;
import org.feno.estimate.client.model.TwoCompartment;

import com.google.gwt.cell.client.NumberCell;
import com.google.gwt.cell.client.TextCell;
import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.Scheduler;
import com.google.gwt.core.client.Scheduler.ScheduledCommand;
import com.google.gwt.dom.client.Style.Cursor;
import com.google.gwt.dom.client.Style.Unit;
import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.event.logical.shared.ResizeEvent;
import com.google.gwt.typedarrays.client.Float64ArrayNative;
import com.google.gwt.user.cellview.client.Column;
import com.google.gwt.user.cellview.client.DataGrid;
import com.google.gwt.user.client.ui.*;
import com.googlecode.gwt.charts.client.*;
import com.googlecode.gwt.charts.client.corechart.Histogram;
import com.googlecode.gwt.charts.client.corechart.HistogramOptions;
import com.googlecode.gwt.charts.client.event.ReadyEvent;
import com.googlecode.gwt.charts.client.event.ReadyHandler;
import com.googlecode.gwt.charts.client.options.Legend;
import com.googlecode.gwt.charts.client.options.LegendPosition;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class Launch implements EntryPoint {
  private DockLayoutPanel appPanel;
  private LayoutPanel plotPanel = new LayoutPanel();
  private final int steps = 100000;
  private final int plotSampleRatio = 4;
  private HashMap<String, Float64ArrayNative> parameterSamples;
  private TwoCompartment twoCompModel = new TwoCompartment(
      new Uniform(0,100), 2, 0.75,
      new Uniform(-5000,5000), 800, 75,
      new Uniform(-500,500), 5, 7.5);
  private final PopupPanel popup = new PopupPanel(false, true); // Create a modal dialog box that will not auto-hide
  private final HTML popupLabel = new HTML();
  private final SimplePanel coverup = new SimplePanel();
  private ChartLoader chartLoader = new ChartLoader(ChartPackage.CORECHART);
  private final MersenneTwister mtRng = new MersenneTwister();
  private HashMap<String, Histogram> currentHistograms = new HashMap<String, Histogram>();
  private HashMap<String, HistogramOptions> currentHistogramOptions = new HashMap<String, HistogramOptions>();

  /**
   * This is the entry point method.
   */
  public void onModuleLoad() {

    popup.add(popupLabel);
    popup.setGlassEnabled(true); // Enable the glass panel
    //popup.getElement().getStyle().setCursor(Cursor.WAIT);
    coverup.getElement().getStyle().setBackgroundColor("#FFFFFF");

    FlowPanel topPanel = new FlowPanel();
    topPanel.add(new HTML("<h1>Welcome to FeNOnline</h1>"));
    getAppPanel().addNorth(topPanel, 6);

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
    OptionPair chainSteps = new OptionPair("Steps", steps);
    optionPairs.add(chainSteps);
    optionsGrid.addColumn(nameCol, "Option Name");
    optionsGrid.addColumn(valueCol, "Current Value");
    optionsGrid.setRowData(optionPairs); 

    LayoutPanel leftPanel = new LayoutPanel();
    Button runDemo = new Button("Run Demo");
    leftPanel.add(runDemo);
    leftPanel.setWidgetTopHeight(runDemo, 0, Unit.EM, 2, Unit.EM);
    leftPanel.setWidgetLeftWidth(runDemo, 30, Unit.PCT, 40, Unit.PCT);
    leftPanel.add(optionsGrid);
    leftPanel.setWidgetTopHeight(optionsGrid, 4, Unit.EM, 100, Unit.PCT);

    getAppPanel().addWest(leftPanel, 20);
    getAppPanel().add(plotPanel);
    RootLayoutPanel.get().add(getAppPanel());   

    runDemo.addClickHandler(new ClickHandler() {
      @Override
      public void onClick(ClickEvent event) {
        Util.showWaitCursor();
        popupLabel.setHTML("<h2>Generating MCMC samples...</h2>");
        popup.center(); // Center the popup and make it visible
        for(Histogram hist : currentHistograms.values()) {
          hist.removeFromParent();
        }
        Scheduler.get().scheduleDeferred(new ScheduledCommand() {    
          @Override
          public void execute() {
            demoCalculate();
//            Scheduler.get().scheduleDeferred(new ScheduledCommand() {    
//              @Override
//              public void execute() {
//                popupLabel.setHTML("<h2>Plotting MCMC samples...</h2>");
//                plotPanel.clear();
//                demoPlot();
//              }
//            });  
          }
        });
      }
    });
  }

  private void demoCalculate(){
    HashMap<String, PosteriorSample> paramSamples = twoCompModel.samplePosteriorToo(steps, mtRng);
    String text = "";
    for(Parameter param : twoCompModel.getParameters()) {
      text += "Param " + param.getName() + "Mean " + paramSamples.get(param.getName()).getSampleStats().getMean();
    }
    popupLabel.setText(text);
    //parameterSamples  = new HashMap<String, Float64ArrayNative>(twoCompModel.getParameters().size());
    //parameterSamples = twoCompModel.samplePosterior(steps, mtRng);
    //    Scheduler.get().scheduleDeferred(new ScheduledCommand() {    
    //      @Override
    //      public void execute() {
    //        Scheduler.get().scheduleDeferred(new ScheduledCommand() {    
    //          @Override
    //          public void execute() {
    //            parameterSamples = twoCompModel.samplePosterior(steps, mtRng);
    //          }
    //        });
    //      }
    //    });
  }


  private void demoPlot() {
    // Create the API Loader
    //ChartLoader chartLoader = new ChartLoader(ChartPackage.CORECHART);
    chartLoader.loadApi(new Runnable() {
      double[] leftPCTOffset = {0, 50, 25};
      double[] topPCTOffset = {0, 0, 50};

      @Override
      public void run() {
        currentHistograms.clear();
        currentHistogramOptions.clear();
        int panelIndex = 0;
        for(Parameter param : twoCompModel.getParameters()) {
          Histogram hist = new Histogram();
          HistogramOptions options = HistogramOptions.create();
          options.setTitle("Posterior sample for " + param.getName());
          options.setLegend(Legend.create(LegendPosition.NONE));
          currentHistograms.put(param.getName(), hist);
          currentHistogramOptions.put(param.getName(), options);
          //hist.setVisible(false);
          plotPanel.add(hist);
          plotPanel.setWidgetLeftWidth(hist, leftPCTOffset[panelIndex], Unit.PCT, 50, Unit.PCT);
          plotPanel.setWidgetTopHeight(hist, topPCTOffset[panelIndex], Unit.PCT, 50, Unit.PCT);
          //hist.draw(createDataTable(param), options);
          panelIndex += 1;
        }
        plotPanel.add(coverup);
        Scheduler.get().scheduleDeferred(new ScheduledCommand() {    
          @Override
          public void execute() {
            for(Parameter param : twoCompModel.getParameters()) {
              currentHistograms.get(param.getName())
              .draw(createDataTable(param.getName()), currentHistogramOptions.get(param.getName()));
            }
            Scheduler.get().scheduleDeferred(new ScheduledCommand() {    
              @Override
              public void execute() {
                plotPanel.remove(coverup);
                popup.hide();
                Util.showDefaultCursor();
              }
            });
          }
        });
      }
    });
  }

  private DockLayoutPanel getAppPanel() {
    if (appPanel == null) {
      appPanel = new DockLayoutPanel(Unit.EM);
    }
    return appPanel;
  }

  private DataTable createDataTable(String paramName) {
    Float64ArrayNative sampleArray = parameterSamples.get(paramName);
    DataTable plotDataTable = DataTable.create();
    int valueCol = plotDataTable.addColumn(ColumnType.NUMBER, paramName);
    int nRows = (int) Math.floor(steps/plotSampleRatio);
    plotDataTable.addRows(nRows);
    for(int i=0; i<nRows; i++) {
      plotDataTable.setCell(i, valueCol, sampleArray.get(i));
    }
    return plotDataTable;
  }

}
