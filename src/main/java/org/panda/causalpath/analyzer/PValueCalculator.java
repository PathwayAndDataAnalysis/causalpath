package org.panda.causalpath.analyzer;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.RankingAlgorithm;
import org.panda.causalpath.resource.ProteomicsLoader;
import org.panda.resource.KinaseLibrary;

import java.util.HashMap;
import java.util.List;

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

    public abstract HashMap<String, Double> getPValues();

    public abstract HashMap<String, Double> getSignedPValues();

    HashMap<String, Double> pValueMap;

    HashMap<String, Double> signedPValueMap;


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

    public void setKinaseColumnMappings(){
        // 1 is the first valid column, as the first column is reserved for peptide change values
        int column = 1;
        for(String kinase: kN.kinaseSet()){
            kinaseColumn.put(kinase, column);
            columnKinase.put(column, kinase);
            column++;
        }
    }

    public void initializeMatrix(){

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







}
