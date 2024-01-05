package org.panda.causalpath.analyzer;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.panda.causalpath.resource.ProteomicsLoader;
import org.panda.utility.statistics.FishersExactTest;

import java.util.*;



public class FischerTestPValueCalculator extends PValueCalculator {



    /**
     * Main goal of constructor is to initialize matrix, a 2-d array
     * where the first column is the change value for each peptide, and
     * every other column is each kinase's score for every peptide
     *
     * @param pL
     */
    public FischerTestPValueCalculator(ProteomicsLoader pL) {
        super(pL);
        setPValues();
    }

    // Method takes two columns, sorts objects by the thresholds provided, and calculates the fischer test value
    public double pValueByThreshold(int col1, int col2, int col1Threshold, int col2Threshold){
        return 0;
    }

    public void setPValues(){
        for(int column = 1; column < originalMatrix[0].length; column++){
            System.out.println("CURRENT COLUMN IS: " + column);
            double pVal = pValueByPercentile(0, column, 0.7, 0.7);
            pValueMap.put(columnKinase.get(column), pVal);
        }
    }

    public double CDFInflectionPoint(double delta){
        double[] thresholds = getThresholds(delta, originalMatrix);
        int numDataPoints = (int) thresholds[0];

        double[] differences = new double[thresholds.length - 1];

        double[] secondDifferentials = new double[differences.length - 1];

        for(int i = 0; i < thresholds.length; i++){
            thresholds[i] = numDataPoints - ((int) thresholds[i]);

        }

        for(int i = 0; i < thresholds.length - 1; i++){
            differences[i] = thresholds[i] - thresholds[i + 1];
        }

        double inflectionPoint = Double.POSITIVE_INFINITY;
        int inflectionIndex = -1;

        for(int i = 0; i < differences.length - 1; i++){
            thresholds[i] = differences[i] - differences[i + 1];
            if(Math.abs(thresholds[i]) < Math.abs(inflectionPoint)){
                System.out.print(thresholds[i] + " ");
                inflectionPoint = thresholds[i];
                inflectionIndex = i;
            }
        }

        return inflectionIndex * delta;


    }

    public double findThresholdSlope(double delta){
        double[] thresholds = getThresholds(delta, originalMatrix);
        int numDataPoints = (int) thresholds[0];

        double[] differences = new double[thresholds.length - 1];

        for(int i = 0; i < thresholds.length; i++){
            thresholds[i] = numDataPoints - ((int) thresholds[i]);

        }

        for(int i = 0; i < thresholds.length - 1; i++){
            differences[i] = Math.abs(thresholds[i] - thresholds[i + 1]);
        }

        int maxIndex = 0;
        double maxValue = Double.NEGATIVE_INFINITY;
        for(int i = 0; i < differences.length; i++){
            if(differences[i] > maxValue){
                maxValue = differences[i];
                maxIndex = i;
            }
        }

        for(int i = 0; i < differences.length; i++){
            System.out.println("From " + i * delta + " to " + ((i * delta) + delta) + " difference value " + differences[i]);
        }

        return maxIndex * delta;

    }

    public double[] getThresholds(double delta, double[][] originalMatrix){
        int numberOfColumns = originalMatrix[0].length;
        int numberOfRows = originalMatrix.length;
        // First column includes the change-values, so we must subtract one
        int numScoreData = (numberOfColumns - 1) * numberOfRows;
        double[] data = new double[numScoreData];
        int cap = 0;
        // Get all kinase scores and sort them
        for(int column = 1; column < originalMatrix[0].length; column++){
            for(int row = 0; row < originalMatrix.length; row++){
                data[cap] = originalMatrix[row][column];
                cap++;
            }
        }

        Arrays.sort(data);
        System.out.println(data.length);
        // Start with the largest score, for every 0.5 decrease in score,
        // calculate the number of datapoints to the right of the threshold

        // The maximum score, the starting threshold
        double maxScore = data[numScoreData - 1];


        // Guarantee that firstThreshold is a multiple of delta
        double firstThreshold = (maxScore % delta == 0) ? maxScore : (Math.floor(maxScore/delta) * delta);

        System.out.println("First threshold " + firstThreshold);

        // To store all thresholds
        int numThresholds = ((int)(firstThreshold/delta)) + 1;
        double[] thresholds = new double[numThresholds];

        // The last index of thresholds will correspond to the greatest threshold
        int thresholdsIndex = thresholds.length - 1;
        // This variable represents the first value in data that is to the left of the current threshold
        int outsideThresh = data.length - 1;
        // A counter maintaining how many data points are greater than the current threshold
        int greaterThanCurrThresh = 0;

        System.out.println(thresholds.length + " " + firstThreshold/delta);

        for(double currentThreshold = firstThreshold; currentThreshold >= 0; currentThreshold = currentThreshold - delta){
            System.out.println(currentThreshold + " " + thresholdsIndex);

            /*
            Given the current threshold, examine the largest datapoint that was less than the previous threshold
            If it is greater than the current threshold, count it, and then repeat this process by checking the
            next smallest kinase score to see if it is greater than the currentThreshold
             */
            while (data[outsideThresh] >= currentThreshold) {
                greaterThanCurrThresh++;
                outsideThresh--;

                // If this is true, we have found the rightmost threshold that encompasses
                // 100% of the data points, so the remaining thresholds will all have greaterThanCurrThresh =
                // all data points
                if(outsideThresh == -1){
                    // Iterate over the remaining thresholds and put in all data points as greaterThanCurrThresh
                    for(int curr = thresholdsIndex; curr >= 0; curr--){
                        thresholds[curr] = greaterThanCurrThresh;
                    }
                    // Completed procedure

                    // Delete for later
                    for(int i = 0; i < thresholds.length; i++){
                        System.out.println("For a threshold of: " + i * delta + " the number greater is " + thresholds[i]);
                    }
                    return thresholds;
                }

            }
            // Assign the number of values counted to the right of the current threshold
            thresholds[thresholdsIndex] = greaterThanCurrThresh;
            thresholdsIndex--;



        }
        // Delete for later
        for(int i = 0; i < thresholds.length; i++){
            System.out.println("For a threshold of: " + i * delta + " the number greater is " + thresholds[i]);
        }
        return thresholds;
    }



    public double findThreshold(double delta){

        double[] thresholds = getThresholds(delta, originalMatrix);





        // For each pair of consecutive measurements take their difference.
        // This difference represents the number of added data points for a 0.5 decrease in score
        // If this measurement is maximum then this is probably a good threshold

        double maxDifference = Double.NEGATIVE_INFINITY;
        int bestThreshold = -1;
        for(int maxThreshIndex = thresholds.length - 1; maxThreshIndex >= 1; maxThreshIndex--){
            double currDifference = Math.abs(thresholds[maxThreshIndex] - thresholds[maxThreshIndex - 1]);
            if(currDifference > maxDifference){
                maxDifference = currDifference;
                bestThreshold = maxThreshIndex;
            }
        }

        return bestThreshold * delta;
    }

    /**
     *  Method takes two columns, objects will be sorted into table based on the percentile of attribute
     * @param peptideChangeColumn
     * @param kinaseScoreColumn
     * @param highScorePercentThreshold You must be above this percent in order to be considered a high score
     * @param increasedPercentThreshold You must be above this percent in order to be increased
     * @return
     */
    public double pValueByPercentile(int peptideChangeColumn, int kinaseScoreColumn,
                                     double highScorePercentThreshold, double increasedPercentThreshold){

        // Variables to keep track of table values
        int increasedPeptideHighScore = 0;
        int increasedPeptideLowScore = 0;
        int notIncPeptideHighScore = 0;
        int notIncPeptideLowScore = 0;

        // A list to accumulate every pair of attributes (treat them as a single object)
        List<Pair> objects = new ArrayList<Pair>();
        // A comparator to sort based on first attribute
        byCol1Value c1 = new byCol1Value();
        byCol2Value c2 = new byCol2Value();

        for(int row = 0; row < originalMatrix.length; row++){
            objects.add(new Pair(originalMatrix[row][peptideChangeColumn], originalMatrix[row][kinaseScoreColumn]));
        }

        // Sort to easily take desired percentile
        Collections.sort(objects, c1);

        double quantityHighScorers = objects.size() * (1 - highScorePercentThreshold);
        double quantityIncreasedPeptides = objects.size() * (1 - increasedPercentThreshold);

        quantityHighScorers = Math.round(quantityHighScorers);
        quantityIncreasedPeptides = Math.round(quantityIncreasedPeptides);

        double quantityLowScorers = objects.size() - quantityHighScorers;
        double quantityNotIncreasedPeptides = objects.size() - quantityIncreasedPeptides;

        HashSet<Pair> notIncreasedPeptides = new HashSet<>();

        // Get the objects with notIncreasedPeptides into a map
        // If you are not a member of this map, that object has an increased peptide
        for(int i = 0; i < quantityNotIncreasedPeptides; i++){
            notIncreasedPeptides.add(objects.get(i));
        }

        // Sort by the kinase scores of the objects
        Collections.sort(objects, c2);

        for(int i = 0; i < objects.size(); i++){
            // Quantity low scorers is the number of objects that will be considered low scorers
            // based on the threshold provided by the user. For example, if there are 10 objects,
            // and the user specifies that 70th percentile and above is high-scores, then quantity
            // of low scorers would be 7. Because we have sorted by score, we can simply obtain the
            // quantityLowScorers lowest scores, and these are the low-scorers based on the user provided
            // percentile
            if(i < quantityLowScorers){
                // Now that we know that this object is a low-scorer, is the corresponding peptide
                // increased or not increased (simply check the hashmap created earlier)
                if(notIncreasedPeptides.contains(objects.get(i))){
                    notIncPeptideLowScore++;
                }
                else{
                    increasedPeptideLowScore++;
                }
            }
            else{
                // The first quantityLowScorer objects are objects with low scorers,
                // everything past that is a high scorer
                // We will simply check whether the object corresponds to a not
                // increased peptide or an increased peptide
                if(notIncreasedPeptides.contains(objects.get(i))){
                    notIncPeptideHighScore++;
                }
                else{
                    increasedPeptideHighScore++;
                }
            }
        }

        System.out.println(increasedPeptideHighScore + " " + increasedPeptideLowScore + " " + notIncPeptideHighScore + " " + notIncPeptideLowScore);

        return FishersExactTest.calcPositiveDepPval(increasedPeptideHighScore, increasedPeptideLowScore, notIncPeptideHighScore, notIncPeptideLowScore);






    }

    @Override
    public HashMap<String, Double> getPValues() {
        return pValueMap;
    }

    @Override
    public HashMap<String, Double> getSignedPValues() {
        return null;
    }


    public static class Pair{
        double col1Value;
        double col2Value;

        public Pair(double col1Value, double col2Value){
            this.col1Value = col1Value;
            this.col2Value = col2Value;
        }
    }

    public static class byCol1Value implements Comparator<Pair> {

        @Override
        public int compare(Pair o1, Pair o2) {
            return Double.compare(o1.col1Value, o2.col1Value);
        }
    }

    public static class byCol2Value implements Comparator<Pair>{
        @Override
        public int compare(Pair o1, Pair o2) {
            return Double.compare(o1.col2Value, o2.col2Value);
        }


    }
}
