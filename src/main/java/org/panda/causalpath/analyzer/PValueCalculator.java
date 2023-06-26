package org.panda.causalpath.analyzer;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.panda.causalpath.resource.ProteomicsLoader;
import org.panda.resource.KinaseLibrary;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

public abstract class PValueCalculator {
    ProteomicsLoader pL;
    KinaseLibrary kN;
    // A mapping to keep track of which kinase corresponds to which column
    HashMap<String, Integer> kinaseColumn;
    // A mapping to keep track of which column corresponds to which kinase
    HashMap<Integer, String> columnKinase;

    public double[][] matrix;

    // A real matrix containing the p-values generated
    RealMatrix r;

    /**
     * Method returns a map where each key is a kinase
     * and each value is its corresponding p-value
     * @return Map containing each kinase's unsigned p-value
     */
    public abstract HashMap<String, Double> getPValues();

    /**
     * Method returns a map where each key is a kinase
     * and each value is its corresponding signed p-value.
     * No sign is activated; whereas a negative sign means inactive
     * @return
     */
    public abstract HashMap<String, Double> getSignedPValues();

    HashMap<String, Double> pValueMap;

    HashMap<String, Double> signedPValueMap;

    /**
     * Main goal of constructor is to initialize matrix, a 2-d array
     * where the first column is the change value for each peptide, and
     * every other column is each kinase's score for every peptide
     * @param pL
     */
    public PValueCalculator(ProteomicsLoader pL){

        this.pL = pL;
        this.kN = new KinaseLibrary();

        int numRows = pL.getSeqChangeValMap().size();
        // One column reserved for the peptide change values
        int numColumns = kN.kinaseSet().size() + 1;

        matrix = new double[numRows][numColumns];

        // Initialize maps that keep track of which column refers to which kinase
        kinaseColumn = new HashMap<>();
        columnKinase = new HashMap<>();

        PearsonsCorrelation pC = new PearsonsCorrelation();

        setKinaseColumnMappings();

        initializeMatrix();





        convertToRanks();

        pValueMap = new HashMap<>();
        signedPValueMap = new HashMap<>();


    }

    private void setKinaseColumnMappings(){
        // 1 is the first valid column, as the first column is reserved for peptide change values
        int column = 1;
        for(String kinase: kN.kinaseSet()){
            kinaseColumn.put(kinase, column);
            columnKinase.put(column, kinase);
            column++;
        }
    }

    private void initializeMatrix(){

        HashMap<String, Double> seqChangeVal = pL.getSeqChangeValMap();

        int rowNum = 0;

        for(String sequence: seqChangeVal.keySet()){
            // Get the change value for this row
            Double peptideChangeValue = seqChangeVal.get(sequence);
            // Put it in the first column
            matrix[rowNum][0] = peptideChangeValue;
            // For the current sequence, obtain every kinase's score for this sequence (one row)
            HashMap<String, Double> kScores = kN.peptideScore(sequence);

            for(String kinase: kScores.keySet()){
                int column = kinaseColumn.get(kinase);
                matrix[rowNum][column] = kScores.get(kinase);
            }

            rowNum++;

        }
    }

    public void convertToRanks(){
        // An object used for converting a column to ranks
        RankingAlgorithm rankingObject = new NaturalRanking();


        for(int column = 0; column < matrix[0].length; column++){

            // Create an array to hold the column data
            double[] ranked = new double[matrix.length];

            // Extract one column of the matrix and then convert it to ranks
            for(int row = 0; row < matrix.length; row++){
                ranked[row] = matrix[row][column];
            }

            int row = 0;

            // Take the ranks and reinput them back into matrix
            for(double d: rankingObject.rank(ranked)){
                matrix[row][column] = d;
                row++;
            }
        }



    }

    /**
     * Not for use in program. Each class that extends PValueCalculator uses a unique approach
     * to generate the p-Values. This method can be used to compare the p-Values calculated
     * through this PValueCalculator's methodology and another PValueCalculator's methodology.
     * For each kinase's p-Value produced via the two methods(a,b) it will plot this point onto
     * the scatterplot with the expectation that these fall on the line y = x
     * @param otherPCalc The other PValueCalculator object
     * @param directory The directory to which you want the scatter plot to be saved
     */
    public void plotPValues(PValueCalculator otherPCalc, String directory){
        PValueCalculator.pValuePlotter pP = this.new pValuePlotter(otherPCalc);
        pP.plotAndSaveChart("P-Values Plotted", directory);
    }


    private class pValuePlotter{
        PValueCalculator otherPCalc;

        private final double[] xData;
        private final double[] yData;
        public pValuePlotter(PValueCalculator otherPCalc){
            this.otherPCalc = otherPCalc;

            Set<String> kinases = pValueMap.keySet();
            xData = new double[kinases.size()];
            yData = new double[kinases.size()];

            int cap = 0;

            for(String kinase: pValueMap.keySet()){
                xData[cap] = pValueMap.get(kinase);
                yData[cap] = otherPCalc.pValueMap.get(kinase);
                cap++;
            }
        }




        private void plotAndSaveChart(String fileName, String directory) {
            XYSeries series = createSeries();
            XYSeriesCollection dataset = new XYSeriesCollection(series);
            JFreeChart chart = createChart(dataset);

            saveChartAsPNG(chart, fileName, directory);
        }

        private XYSeries createSeries() {
            XYSeries series = new XYSeries("Data Series");
            for (int i = 0; i < xData.length; i++) {
                series.add(xData[i], yData[i]);
            }
            return series;
        }

        private JFreeChart createChart(XYSeriesCollection dataset) {
            JFreeChart chart = ChartFactory.createScatterPlot(
                    "Chart Title", // Chart title
                    "X Axis", // X-axis label
                    "Y Axis", // Y-axis label
                    dataset, // Dataset
                    PlotOrientation.VERTICAL, // Plot orientation
                    true, // Include legend
                    true, // Include tooltips
                    false // Include URLs
            );
            return chart;
        }

        private void saveChartAsPNG(JFreeChart chart, String fileName, String directory) {
            try {
                ChartUtils.saveChartAsPNG(new File(directory, fileName), chart, 800, 600);
                System.out.println("Chart saved as PNG successfully.");
            } catch (IOException e) {
                System.err.println("Error saving chart as PNG: " + e.getMessage());
            }
        }


    }







}
