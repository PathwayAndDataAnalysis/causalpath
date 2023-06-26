package org.panda.causalpath.analyzer;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.panda.causalpath.resource.ProteomicsLoader;

import java.util.HashMap;

public class LibraryPValueCalculator extends PValueCalculator{

    PearsonsCorrelation pC;

    public LibraryPValueCalculator(ProteomicsLoader pL) {
        super(pL);
        pC = new PearsonsCorrelation(matrix);
        setPValues();
    }

    private void setPValues(){
        RealMatrix r = pC.getCorrelationMatrix();
        RealMatrix r1 = pC.getCorrelationPValues();

        int numColumns = matrix[0].length;

        // Iterate over all columns corresponding to kinase scores
        // and obtain the correlation coeffecient (first column is peptide change values)
        for(int i = 1; i < numColumns; i++){
            // Find out which kinase this column corresponds to
            String currKinase = columnKinase.get(i);
            double correlationCoeff = r.getEntry(0, i);
            double pValue = r1.getEntry(0, i);
            // Negative correlation coeffecient => inactive in test vs control
            double signedPValue = correlationCoeff >= 0 ? pValue: pValue * -1;
            pValueMap.put(currKinase, pValue);
            signedPValueMap.put(currKinase, signedPValue);


        }
    }

    @Override
    public HashMap<String, Double> getPValues() {
        return pValueMap;
    }

    @Override
    public HashMap<String, Double> getSignedPValues() {
        return signedPValueMap;
    }
}
