package org.panda.causalpath.analyzer;

import org.panda.causalpath.resource.ProteomicsLoader;
import org.panda.utility.statistics.FDR;

import java.io.File;
import java.io.FileWriter;
import java.io.Writer;
import java.util.Map;

public class KinaseActivityAnalyzer {
    private PValueCalculator pCalc;
    private Map<String, Double> signedPValues;
    public KinaseActivityAnalyzer(ProteomicsLoader pL, String directory){
        pCalc = new LibraryPValueCalculator(pL);
        signedPValues = pCalc.getSignedPValues();
        resultsToFile(directory);
    }
    public void resultsToFile(String directory){

        // qVals is a map for each kinase, and it will calculate the FDR
        // assuming that the kinase's pvalue is used as the threshold
        Map<String, Double> qVals= FDR.getQVals(pCalc.pValueMap, null);
        File file = new File(directory, "kinaseAnalysis");
        try{
            Writer w = new FileWriter(file);
            w.write("Kinase" + "\t" + "pValues" +"\t" + "QValue" + "\n");
            for(String kinase: signedPValues.keySet()){
                w.write(kinase + "\t" + signedPValues.get(kinase) +"\t" + qVals.get(kinase));
                w.write("\n");
            }

            w.close();
        }
        catch(Exception e){
            System.out.println("Error occured when writing to file: " + e.getMessage());
        }

    }
}
