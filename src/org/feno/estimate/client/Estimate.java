package org.feno.estimate.client;

import jdistlib.Normal;
import jdistlib.math.density.Density;
import jdistlib.rng.MersenneTwister;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.RootLayoutPanel;
import com.google.gwt.user.client.ui.SimpleLayoutPanel;
import com.googlecode.gwt.charts.client.ChartLoader;
import com.googlecode.gwt.charts.client.ChartPackage;
import com.googlecode.gwt.charts.client.ColumnType;
import com.googlecode.gwt.charts.client.DataTable;
import com.googlecode.gwt.charts.client.corechart.Histogram;
import com.googlecode.gwt.charts.client.corechart.HistogramOptions;
import com.googlecode.gwt.charts.client.corechart.ScatterChart;
import com.googlecode.gwt.charts.client.options.Legend;
import com.googlecode.gwt.charts.client.options.LegendPosition;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class Estimate implements EntryPoint {
  private SimpleLayoutPanel layoutPanel;
  private Histogram histogram;
  private ScatterChart scatterChart;
  
  public native void debugger() /*-{
    debugger;
  }-*/;

  /**
   * This is the entry point method.
   */
  public void onModuleLoad() {
    Window.enableScrolling(false);
    Window.setMargin("0px");

    RootLayoutPanel.get().add(getSimpleLayoutPanel());

    // Create the API Loader
    ChartLoader chartLoader = new ChartLoader(ChartPackage.CORECHART);
    chartLoader.loadApi(new Runnable() {

      @Override
      public void run() {
        getSimpleLayoutPanel().setWidget(getScatterChart());
        drawHistogram();
      }
    });
  }
  
  private SimpleLayoutPanel getSimpleLayoutPanel() {
    if (layoutPanel == null) {
      layoutPanel = new SimpleLayoutPanel();
    }
    return layoutPanel;
  }
  
  private Histogram getHistogram() {
    if (histogram == null) {
      histogram = new Histogram();
    }
    return histogram;
  }
  
  private ScatterChart getScatterChart() {
    if (scatterChart == null) {
      scatterChart = new ScatterChart();
    }
    return scatterChart;       
  }
  
  private void drawHistogram() {
    int nDoubles = 100;
    double[] someDoubles = new double[nDoubles];
    MersenneTwister mt = new MersenneTwister(31415);
    for(int i=1; i<nDoubles; i++) {
      someDoubles[i] = Normal.random_standard(mt);
      //RootPanel.get().add(new HTML("<div>"+String.valueOf(someDoubles[i])+"</div>"));
    }
    //debugger();
    Density sampleDensity = Density.density(someDoubles);
    DataTable densityData = DataTable.create();
    int xCol = densityData.addColumn(ColumnType.NUMBER, "xValues");
    int densityCol = densityData.addColumn(ColumnType.NUMBER, "densityValues");
    for(int i=0; i<sampleDensity.x.length; i++) {
      int row = densityData.addRow();
      densityData.setCell(row, xCol, sampleDensity.x[i]);
      densityData.setCell(row, densityCol, sampleDensity.y[i]);
    }
    
    scatterChart.draw(densityData);
    //DataTable dataTable = ChartHelper.arrayToDataTable(someDoubles, true);

    // Set options
    HistogramOptions options = HistogramOptions.create();
    options.setTitle("A histogram of " + String.valueOf(nDoubles) + " standard normal sample points");
    options.setLegend(Legend.create(LegendPosition.IN));

    // Draw the chart
    //histogram.draw(dataTable, options);
  }
  

//    int nDoubles = 10;
//    double[] someDoubles = new double[nDoubles];
//    MersenneTwister mt = new MersenneTwister(314);
//    for(int i=1; i<nDoubles; i++) {
//      someDoubles[i] = Normal.random(1, 4, mt);
//      RootPanel.get().add(new HTML("<div>"+String.valueOf(someDoubles[i])+"</div>"));
//    }
//    final Button sendButton = new Button("Send");
//    final TextBox nameField = new TextBox();
//    nameField.setText("GWT User");
//    final Label errorLabel = new Label();
//
//    // We can add style names to widgets
//    sendButton.addStyleName("sendButton");
//
//    // Add the nameField and sendButton to the RootPanel
//    // Use RootPanel.get() to get the entire body element
//    RootPanel.get("nameFieldContainer").add(nameField);
//    RootPanel.get("sendButtonContainer").add(sendButton);
//    RootPanel.get("errorLabelContainer").add(errorLabel);
//
//    // Focus the cursor on the name field when the app loads
//    nameField.setFocus(true);
//    nameField.selectAll();
//
//    // Create the popup dialog box
//    final DialogBox dialogBox = new DialogBox();
//    dialogBox.setText("Remote Procedure Call");
//    dialogBox.setAnimationEnabled(true);
//    final Button closeButton = new Button("Close");
//    // We can set the id of a widget by accessing its Element
//    closeButton.getElement().setId("closeButton");
//    final Label textToServerLabel = new Label();
//    final HTML serverResponseLabel = new HTML();
//    VerticalPanel dialogVPanel = new VerticalPanel();
//    dialogVPanel.addStyleName("dialogVPanel");
//    dialogVPanel.add(new HTML("<b>Sending name to the server:</b>"));
//    dialogVPanel.add(textToServerLabel);
//    dialogVPanel.add(new HTML("<br><b>Server replies:</b>"));
//    dialogVPanel.add(serverResponseLabel);
//    dialogVPanel.setHorizontalAlignment(VerticalPanel.ALIGN_RIGHT);
//    dialogVPanel.add(closeButton);
//    dialogBox.setWidget(dialogVPanel);
//
//    // Add a handler to close the DialogBox
//    closeButton.addClickHandler(new ClickHandler() {
//      public void onClick(ClickEvent event) {
//        dialogBox.hide();
//        sendButton.setEnabled(true);
//        sendButton.setFocus(true);
//      }
//    });
//
//    // Create a handler for the sendButton and nameField
//    class MyHandler implements ClickHandler, KeyUpHandler {
//      /**
//       * Fired when the user clicks on the sendButton.
//       */
//      public void onClick(ClickEvent event) {
//        sendNameToServer();
//      }
//
//      /**
//       * Fired when the user types in the nameField.
//       */
//      public void onKeyUp(KeyUpEvent event) {
//        if (event.getNativeKeyCode() == KeyCodes.KEY_ENTER) {
//          sendNameToServer();
//        }
//      }
//
//      /**
//       * Send the name from the nameField to the server and wait for a response.
//       */
//      private void sendNameToServer() {
//        // First, we validate the input.
//        errorLabel.setText("");
//        String textToServer = nameField.getText();
//        if (!FieldVerifier.isValidName(textToServer)) {
//          errorLabel.setText("Please enter at least four characters");
//          return;
//        }
//        
//        // Then, we send the input to the server.
//        sendButton.setEnabled(false);
//        textToServerLabel.setText(textToServer);
//        serverResponseLabel.setText("");
//        greetingService.greetServer(textToServer, new AsyncCallback<String>() {
//          public void onFailure(Throwable caught) {
//            // Show the RPC error message to the user
//            dialogBox.setText("Remote Procedure Call - Failure");
//            serverResponseLabel.addStyleName("serverResponseLabelError");
//            serverResponseLabel.setHTML(SERVER_ERROR);
//            dialogBox.center();
//            closeButton.setFocus(true);
//          }
//
//          public void onSuccess(String result) {
//            dialogBox.setText("Remote Procedure Call");
//            serverResponseLabel.removeStyleName("serverResponseLabelError");
//            serverResponseLabel.setHTML(result);
//            dialogBox.center();
//            closeButton.setFocus(true);
//          }
//        });
//      }
//    }
//
//    // Add a handler to send the name to the server
//    MyHandler handler = new MyHandler();
//    sendButton.addClickHandler(handler);
//    nameField.addKeyUpHandler(handler);
//  }
}
