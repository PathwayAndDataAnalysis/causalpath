package org.panda.causalpath.analyzer;

import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.panda.causalpath.resource.ProteomicsLoader;

import java.util.*;

public class EmpiricalPValueCalculator extends PValueCalculator{

    SpearmansCorrelation sP;


    public EmpiricalPValueCalculator(ProteomicsLoader pL) {
        super(pL);
        this.sP = new SpearmansCorrelation(new BlockRealMatrix(matrix));
        initializePValueMap();


    }

    @Override
    public HashMap<String, Double> getPValues() {
        return pValueMap;
    }

    @Override
    public HashMap<String, Double> getSignedPValues() {
        return signedPValueMap;
    }

    private List<Double> getNullDistribution(){
        List<Double> list1 = new ArrayList<>();
        List<Double> list2 = new ArrayList<>();

        for(int i = 0; i < matrix.length; i++){
            list1.add((double)i);
            list2.add((double)i);
        }

        return generateNullDistribution(list1, list2);
    }



    private void initializePValueMap(){

        RealMatrix r = sP.getCorrelationMatrix();

        List<Double> nullDistribution = getNullDistribution();




        for(int i = 1; i < matrix[0].length; i++){
            String currKinase = columnKinase.get(i);
            double kinaseCorrelationCoeff = r.getEntry(0, i);
            double pValue = getPValue(nullDistribution, kinaseCorrelationCoeff);
            pValueMap.put(currKinase, pValue);
            double signedPValue = pValue >= 0 ? pValue: pValue * -1;
            signedPValueMap.put(currKinase, signedPValue);

        }
    }

    private double getPValue(List<Double> nullDistribution, double observation){

        int count = 0;

        for (Double nullVal : nullDistribution) {
            if (Math.abs(nullVal) >= Math.abs(observation)) {
                count++;
            }
        }

        return count/ ((double) nullDistribution.size());

    }

    private List<Double> generateNullDistribution(List<Double> changeValues, List<Double> peptideScoresList){

        ArrayList<Double> correlationCoeffecients = new ArrayList<>();

        for(int i = 0; i < 5000; i++){

            // Produce and shuffle a deep copy; avoid shuffling actual list!
            List<Double> shuffledPeptideScores = copyList(peptideScoresList);
            Collections.shuffle(shuffledPeptideScores);

            // Calculate correlation coeffecient

            double coeff = listCorrelationCoeffecient(changeValues, shuffledPeptideScores);

            correlationCoeffecients.add(coeff);

        }

        Collections.sort(correlationCoeffecients);

        return correlationCoeffecients;
    }

    private List<Double> copyList(List<Double> originalList){

        return new ArrayList<>(originalList);
    }

    private double listCorrelationCoeffecient(List<Double> ls1, List<Double> ls2){
        double[] ls1Arr = listToArray(ls1);
        double[] ls2Arr = listToArray(ls2);

        return sP.correlation(ls1Arr, ls2Arr);
    }

    private double[] listToArray(List<Double> list){
        double[] arr = new double[list.size()];
        int cap = 0;
        for(double d: list){
            arr[cap] = d;
            cap++;
        }

        return arr;
    }
}
