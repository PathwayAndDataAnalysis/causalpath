package org.panda.causalpath.analyzer;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Ellipse2D;

// A class I made to compare p-values produced from comparing to a null distribution
// versus p-values produced from using the library function
public class pValueAnalysis extends JFrame {

    public pValueAnalysis(String title, double[] pVal1, double[] pVal2){
        super(title);

        XYDataset dataset = createDataSet(pVal1, pVal2);

        JFreeChart chart = ChartFactory.createScatterPlot("P-Values calculated in two different ways",
                "X-axis", "Y-axis", dataset);

        //Changes background color
        XYPlot plot = (XYPlot)chart.getPlot();

        XYItemRenderer renderer = plot.getRenderer();

// Set the size and shape of the points for the first series (Boys)
        renderer.setSeriesShape(0, new Ellipse2D.Double(-2, -2, 4, 4));
        renderer.setSeriesPaint(0, Color.BLUE);

// Set the size and shape of the points for the second series (Girls)
        renderer.setSeriesShape(1, new Ellipse2D.Double(-2, -2, 6, 6));
        renderer.setSeriesPaint(1, Color.RED);

        plot.setBackgroundPaint(new Color(255,228,196));


        // Create Panel
        ChartPanel panel = new ChartPanel(chart);
        setContentPane(panel);

    }

    public XYDataset createDataSet(double[] pVal1, double[] pVal2){
        XYSeriesCollection dataSet = new XYSeriesCollection();

        XYSeries series1 = new XYSeries("First Series");

        for(int i = 0; i < pVal1.length; i++){
            series1.add(pVal1[i], pVal2[i]);
        }

        dataSet.addSeries(series1);

        XYSeries series2 = new XYSeries("Second Series");

        for(double i = 0; i < 1; i = i + 0.001){
            series2.add(i, i);
        }

        dataSet.addSeries(series2);


        /*
        XYSeries s2 = new XYSeries("Second series");

        for(int i = 0; i < pVal2.length; i++){
            s2.add(i, pVal2[i]);
        }

        dataSet.addSeries(s2);

         */

        return dataSet;

    }
}
