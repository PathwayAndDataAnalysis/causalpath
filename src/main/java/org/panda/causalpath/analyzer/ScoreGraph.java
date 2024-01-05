package org.panda.causalpath.analyzer;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;
import org.panda.causalpath.resource.ProteomicsLoader;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Used for generating histogram. Probably remove before push
 */
public class ScoreGraph extends PValueCalculator {

    /**
     * Main goal of constructor is to initialize matrix, a 2-d array
     * where the first column is the change value for each peptide, and
     * every other column is each kinase's score for every peptide
     *
     * @param pL
     */

    public ScoreGraph.HistogramGenerator hE;
    public ScoreGraph(ProteomicsLoader pL, String directory) {
        super(pL);
        hE = this.new HistogramGenerator(originalMatrix, 10);
        hE.generateAndSaveHistogram(directory);


    }

    public void chartScores(){

    }

    @Override
    public HashMap<String, Double> getPValues() {
        return null;
    }

    @Override
    public HashMap<String, Double> getSignedPValues() {
        return null;
    }

    public double getMaxScore(){
        double max = 0;
        int row;
        int column;
        int saveLargestRow = 0;
        int saveLargestColumn = 0;
        for(column = 1; column < originalMatrix[0].length; column++){
            for(row = 0; row < originalMatrix.length; row++){
                if(originalMatrix[row][column] > max){
                    max = originalMatrix[row][column];
                    saveLargestRow = row;
                    saveLargestColumn = column;
                }
            }
        }

        String kinase = columnKinase.get(saveLargestColumn);

        System.out.println(kinase);


        return max;
    }

    public double getMinScore(){
        double min = Double.MAX_VALUE;
        for(int column = 1; column < originalMatrix[0].length; column++){
            for(int row = 0; row < originalMatrix.length; row++){
                if(originalMatrix[row][column] < min){
                    min = originalMatrix[row][column];
                }
            }
        }

        return min;
    }










        public class HistogramGenerator {
            private double[][] data;
            private int numberOfBins;
            private double[] allData;

            public HistogramGenerator(double[][] data, int numberOfBins) {
                this.data = data;
                this.numberOfBins = numberOfBins;
            }

            public void generateParticularKinaseHistogram(double[][] data,
                                                          java.util.List<Integer> listOfColumns, String outputDirectory){

                HistogramDataset dataset = new HistogramDataset();

                // Need enough space for all the columns data
                double[] d = new double[data.length * listOfColumns.size()];
                int cap = 0;

                for(int i = 0; i < listOfColumns.size(); i++){
                    int column = listOfColumns.get(i);
                    for(int row = 0; row < data.length; row++){
                        d[cap] = data[row][column];
                        cap++;
                    }


                }



                dataset.addSeries("Kinase Scores", d, numberOfBins);

                // Create the chart
                JFreeChart chart = ChartFactory.createHistogram(
                        "Histogram Example",
                        "Value",
                        "Frequency",
                        dataset,
                        PlotOrientation.VERTICAL,
                        true,
                        true,
                        false);


                NumberAxis xAxis = (NumberAxis) chart.getXYPlot().getDomainAxis();
                xAxis.setRange(0, 1);

                NumberAxis yAxis = (NumberAxis) chart.getXYPlot().getRangeAxis();
                yAxis.setRange(0, 10000);


                // Customize the chart
                chart.setBackgroundPaint(Color.WHITE);


                // Save the chart as a PNG image
                String outputPath = outputDirectory + File.separator + "partialKinaseScoreHistogram.png";
                File outputFile = new File(outputPath);
                try {
                    ChartUtils.saveChartAsPNG(outputFile, chart, 1000, 1000);
                    System.out.println("Chart saved successfully as " + outputPath);
                } catch (IOException e) {
                    System.err.println("Error saving chart: " + e.getMessage());
                }
            }
            public void writeToFile(String directory) {
                String outputPath = directory + File.separator + "histogram_data.csv";
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath))) {
                    // Write the histogram data
                    for (int i = 0; i < allData.length; i++) {
                        writer.write(String.valueOf(allData[i]));
                        writer.newLine();
                    }
                    System.out.println("Histogram data saved successfully.");
                } catch (IOException e) {
                    System.err.println("Error saving histogram data: " + e.getMessage());
                }
            }

            public void generateAndSaveHistogram(String outputDirectory) {
                // Concatenate values from the columns (excluding the first column) into a single array
                int numRows = data.length;
                int numColumns = data[0].length;
                int totalValues = (numRows) * (numColumns-1);
                System.out.println(totalValues);
                allData = new double[totalValues];
                int index = 0;
                for (int j = 1; j < data[0].length; j++) {
                    for (int i = 0; i < data.length; i++) {
                        allData[index++] = data[i][j];
                    }
                }

                // Create a dataset with all the data values


                HistogramDataset dataset = new HistogramDataset();
                dataset.addSeries("Data", allData, numberOfBins);

                // Create the chart
                JFreeChart chart = ChartFactory.createHistogram(
                        "Histogram Example",
                        "Value",
                        "Frequency",
                        dataset,
                        PlotOrientation.VERTICAL,
                        true,
                        true,
                        false);



                NumberAxis xAxis = (NumberAxis) chart.getXYPlot().getDomainAxis();
                xAxis.setRange(0, 5);

                /*
                NumberAxis yAxis = (NumberAxis) chart.getXYPlot().getRangeAxis();
                yAxis.setRange(0, 8000);

                 */


                // Customize the chart
                chart.setBackgroundPaint(Color.WHITE);


                // Save the chart as a PNG image
                String outputPath = outputDirectory + File.separator + "histogram.png";
                File outputFile = new File(outputPath);
                try {
                    ChartUtils.saveChartAsPNG(outputFile, chart, 1000, 1000);
                    System.out.println("Chart saved successfully as " + outputPath);
                } catch (IOException e) {
                    System.err.println("Error saving chart: " + e.getMessage());
                }
            }
        }


}
