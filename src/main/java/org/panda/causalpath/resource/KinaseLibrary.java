package org.panda.causalpath.resource;

import org.panda.resource.FileServer;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class KinaseLibrary extends FileServer {

    /*
    A map which can be used to determine the probability that a kinase has some
    amino acid at some position
     */
    private Map<String, Map<Integer, Map<String, Double>>> probabilityMap;

    private int maximumSiteVal;

    private int minimumSiteVal;

    @Override
    public String[] getLocalFilenames() {
        return new String[]{"KinaseLibraryMatrices.txt"};
    }

    // Returns the minimum location value
    public int minimumLocationValue(){
        return minimumSiteVal;
    }

    // Returns the maximum location value
    public int maximumLocationValue(){
        return maximumSiteVal;
    }


    // Returns all kinases within the map
    public Set<String> kinaseSet(){
        return probabilityMap.keySet();
    }

    /* Public method for user to access the probability that if kinase phosphorlyated
     a given peptide, what are the odds the peptide has given aminoAcid at specified location
     */
    public double aminoAcidProbability(String kinase, int location, String aminoAcid){

        // Check for user requesting information that does not exist in the map
        if(!probabilityMap.containsKey(kinase)){
            throw new IllegalArgumentException("No such kinase available.");
        }
        else if(!probabilityMap.get(kinase).containsKey(location)){
            throw new IllegalArgumentException("No such location available.");
        }
        else if(!probabilityMap.get(kinase).get(location).containsKey(aminoAcid)){
            throw new IllegalArgumentException("No such amino acid available");
        }

        return probabilityMap.get(kinase).get(location).get(aminoAcid);

    }

    @Override
    public boolean load() throws IOException {
        // Get the path to the file
        Path p = (Paths.get(this.locateInBase(this.getLocalFilenames()[0])));

        // Read all the lines to begin processing
        List<String> lines = Files.readAllLines(p);

        // Process the first line in order to get column headers (do not hardcode these)
        String[] firstLineSplit = lines.get(0).split("\t");

        // This next part of the code will find the locations/amino acids from the column headers
        // We will store these in two arrays

        int[] locations = new int[lines.size()];
        String[] aminoAcids = new String[lines.size()];

        splitColumnHeaders(firstLineSplit, locations, aminoAcids);

        // Assign maximum location value and minimum location value
        maximumSiteVal = findMax(locations);
        minimumSiteVal = findMin(locations);


        // Process row-by-row to fill out maps

         /* A 3 layered map, which can be used to access raw data regarding a kinase's
       preference for a specific amino acid at some location
        */
        Map<String, Map<Integer, Map<String, Double>>> kinaseLocationMap = new HashMap<>();

        // Assign probabilityMap
        probabilityMap = new HashMap<>();


        for(int row = 1; row < lines.size(); row++){
            // Split the current row or line by the delimiter tab
            String[] splitLine = lines.get(row).split("\t");

            // Determine the kinase for curr row, and put it into the kinase-Location map
            String currKinase = splitLine[0];

            kinaseLocationMap.put(currKinase, new HashMap<>());
            probabilityMap.put(currKinase, new HashMap<>());

            // Map to keep track of sums of measurements for each location(in order to compute probabilities)
            Map<Integer, Double> sumByLocation = new HashMap<>();

            for(int cell = 1; cell < splitLine.length; cell++){

                // Use the column headers to determine the location/amino acid for current cell
                int currLocation = locations[cell];
                String currAminoAcid = aminoAcids[cell];
                double measurement = Double.parseDouble(splitLine[cell]);

                // Update the sums of measurements for each location
                if(sumByLocation.containsKey(currLocation)){
                    Double partialSum = sumByLocation.get(currLocation);
                    Double newSum = partialSum + measurement;
                    sumByLocation.put(currLocation, newSum);
                }
                else{
                    sumByLocation.put(currLocation, measurement);
                }

                // Put the measurement into the map appropriately
                if(!kinaseLocationMap.get(currKinase).containsKey(currLocation)){
                    kinaseLocationMap.get(currKinase).put(currLocation, new HashMap<>());
                    probabilityMap.get(currKinase).put(currLocation, new HashMap<>());
                }

                kinaseLocationMap.get(currKinase).get(currLocation).put(currAminoAcid, measurement);
                probabilityMap.get(currKinase).get(currLocation).put(currAminoAcid, measurement);

            }

            /* Once we have processed a single row, we may go ahead and calculate the probability that
               an amino acid has occured at a given location
             */

            // Iterate through all locations for this kinase
            for(int location: kinaseLocationMap.get(currKinase).keySet()){
                // For each location for the given kinase, we will iterate through aa-raw measurement map
                Map<String, Double> aminoAcidValue = kinaseLocationMap.get(currKinase).get(location);

                // For each amino acid, compute the probability by dividing corresponding value
                // by the sum of the measurements for that site
                for(String aa: aminoAcidValue.keySet()){
                    double aaMeasurement = probabilityMap.get(currKinase).get(location).get(aa);
                    double aaProbability = aaMeasurement/sumByLocation.get(location);
                    probabilityMap.get(currKinase).get(location).put(aa, aaProbability);
                }


            }
        }

        return true;
    }


    public static void main(String[] args){

        // Load

        KinaseLibrary kN = new KinaseLibrary();
        try{
            kN.load();
        }
        catch(Exception e){

        }

        System.out.println("Is it true that the tests were passed?: " + kN.testHandler());
    }

    /*
    Testing method, takes some segment of the data from excel, it manually calculates the probability, then
    compares it with the value in the map.
     */
    private boolean testProbabilityMap(String Kinase, String[] dataSegment, String[] correspondingHeader){

        // Find the corresponding locations and amino acids for this data segment
        int[] locations = new int[correspondingHeader.length];
        String[] aminoAcids = new String[correspondingHeader.length];
        splitColumnHeaders(correspondingHeader, locations, aminoAcids);

        // A hashmap to maintain the manually computed sum for each location
        HashMap<Integer, Double> sumByLocation = new HashMap<>();


        double sum = 0;

        for(int i = 0; i < dataSegment.length; i++){
            double currMeasurement = Double.parseDouble(dataSegment[i]);

            if(sumByLocation.containsKey(locations[i])){
                Double partialSum = sumByLocation.get(locations[i]);
                Double newSum = partialSum + currMeasurement;
                sumByLocation.put(locations[i], newSum);
            }
            else{
                sumByLocation.put(locations[i], currMeasurement);
            }
        }

        for(int i = 0; i < locations.length; i++){

            double manuallyComputedProbability = Double.parseDouble(dataSegment[i])/sumByLocation.get(locations[i]);

            boolean validError = manuallyComputedProbability - probabilityMap.get(Kinase).get(locations[i]).get(aminoAcids[i]) < 0.0000001 &&
                    manuallyComputedProbability - probabilityMap.get(Kinase).get(locations[i]).get(aminoAcids[i]) > -0.0000001;

            if(!validError){
                System.out.println("Test failure at column " + i + "where cell value is: " + Double.parseDouble(dataSegment[i]) +
                        " and the calculated sum of measurements is " + sumByLocation.get(locations[i]));
                System.out.println(manuallyComputedProbability + " " + probabilityMap.get(Kinase).get(locations[i]).get(aminoAcids[i]));
                System.out.println(manuallyComputedProbability - probabilityMap.get(Kinase).get(locations[i]).get(aminoAcids[i]));
                return false;
            }
        }
        return true;
    }

    public boolean testHandler(){

        // Test #1
        String[] ds1 = {"5683377.86", "9080772.95", "8340840.22", "10883246.95", "8186154.23", "7555351.79", "5841037.6", "6726258.4", "5794366.62", "5560729.33", "6387081.93", "7418541.16", "6976863.18", "8319964.94", "10588197.24", "12875913.33", "6178303.19", "9242426.42", "7731860.37", "8130753.93", "8006580.23", "8006580.23", "10039583.51", "7958867.31", "10857162.09", "9473304.23", "11466629.34", "9472229.52", "7808415.78", "6145822.45", "6211134.39", "6735254.52", "7197377.31", "6378000.66", "6882198.26", "8495568.55", "8370676.86", "11899998.53", "15019256.73", "7990390.94", "9149076.87", "12280572.35", "8466437.59", "8550929.04", "8550929.04", "8461969.98", "7893630.45", "12341329.54", "5944089.31", "9439142.51", "9163620.4", "6872721.84", "3803083.24", "4114322.07", "3281708.07", "4837111.52", "5128332.52", "5133141.96", "6967519.69", "8615566.13", "14118447.67", "43189607.09", "5960485.23", "7663342.47", "5417073.9", "7601312.84", "7860535.17", "7860535.17", "8954818.68", "3594447.85", "3978731.61", "2608708.96", "3950832.04", "9450264.91", "4586496.58", "2023682.92", "2016602.06", "2266725.34", "2188762.8", "2530417.73", "2308419.4", "2883789.22", "4248220.34", "16373437.14", "146172558.6", "7567379.18", "4296821.91", "3407441.02", "3195616.4", "4660708.62", "4660708.62", "4432927.54", "6104437.77", "5618741.13", "5568670", "8329791.91", "8376059.15", "6398860.01", "4904681.74", "5198540.82", "10755585.68", "11285236.69", "7207505.51", "10853139.28", "4791017.71", "8897219.55", "5768673.36", "9428310.25", "8854185.69", "9383652.54", "11303201.82", "12748320.66", "4798936.23", "4798936.23", "18243735.87", "1570050.67", "4817906.3", "3530110.16", "13593499.82", "10018001.58", "6904104.97", "8300042.84", "12391701.34", "15242214.98", "13093979.88", "13596272.37", "8452920.14", "9417994.61", "4952508.91", "1794713.07", "3165114.16", "5469980.93", "6346712.55", "8523094.45", "5924121.94", "11965552.32", "11965552.32", "6221615.33", "6187683", "13923765.79", "7964465.97", "13763518.63", "9211406.25", "5923056.04", "6073157.6", "8137327.3", "6222667.55", "6597541.73", "7595807.16", "21996928.96", "6986969.28", "10348689.08", "2570453.2", "4039190.54", "4992070.17", "5148360.13", "7879384.3", "8182293.36", "10007316.82", "10007316.82", "12171087.94", "12164972.2", "13800089.39", "8459006.06", "14577163.74", "22444340.11", "7800360.93", "7234269.2", "8073261.69", "7576863.68", "7023019.46", "11680945.08", "12666320.57", "11212145.91", "12320805.5", "7874559.56", "7666434.04", "7864492.23", "7452652.89", "7673844.45", "7943308.09", "9963238.5", "9963238.5", "16243733.68", "8179334", "8900275.8", "7283974.75", "11543035.5", "9089075.19", "7825897.21", "5761667.26", "7859246.98", "5475055.46", "8466023.95", "6375030.33", "8001396.18", "10273142.04", "9946645.41", "10051110.92", "10434776.86", "8957963.67", "10614474.87", "7650278.9", "9278532.9", "12892094.54", "12892094.54", "12066131.31"};
        String[] ch1 = {"-5P", "-5G", "-5A", "-5C", "-5S", "-5T", "-5V", "-5I", "-5L", "-5M", "-5F", "-5Y", "-5W", "-5H", "-5K", "-5R", "-5Q", "-5N", "-5D", "-5E", "-5s", "-5t", "-5y", "-4P", "-4G", "-4A", "-4C", "-4S", "-4T", "-4V", "-4I", "-4L", "-4M", "-4F", "-4Y", "-4W", "-4H", "-4K", "-4R", "-4Q", "-4N", "-4D", "-4E", "-4s", "-4t", "-4y", "-3P", "-3G", "-3A", "-3C", "-3S", "-3T", "-3V", "-3I", "-3L", "-3M", "-3F", "-3Y", "-3W", "-3H", "-3K", "-3R", "-3Q", "-3N", "-3D", "-3E", "-3s", "-3t", "-3y", "-2P", "-2G", "-2A", "-2C", "-2S", "-2T", "-2V", "-2I", "-2L", "-2M", "-2F", "-2Y", "-2W", "-2H", "-2K", "-2R", "-2Q", "-2N", "-2D", "-2E", "-2s", "-2t", "-2y", "-1P", "-1G", "-1A", "-1C", "-1S", "-1T", "-1V", "-1I", "-1L", "-1M", "-1F", "-1Y", "-1W", "-1H", "-1K", "-1R", "-1Q", "-1N", "-1D", "-1E", "-1s", "-1t", "-1y", "1P", "1G", "1A", "1C", "1S", "1T", "1V", "1I", "1L", "1M", "1F", "1Y", "1W", "1H", "1K", "1R", "1Q", "1N", "1D", "1E", "1s", "1t", "1y", "2P", "2G", "2A", "2C", "2S", "2T", "2V", "2I", "2L", "2M", "2F", "2Y", "2W", "2H", "2K", "2R", "2Q", "2N", "2D", "2E", "2s", "2t", "2y", "3P", "3G", "3A", "3C", "3S", "3T", "3V", "3I", "3L", "3M", "3F", "3Y", "3W", "3H", "3K", "3R", "3Q", "3N", "3D", "3E", "3s", "3t", "3y", "4P", "4G", "4A", "4C", "4S", "4T", "4V", "4I", "4L", "4M", "4F", "4Y", "4W", "4H", "4K", "4R", "4Q", "4N", "4D", "4E", "4s", "4t", "4y"};
        String k1 = "AURA";
        boolean t1Result = testProbabilityMap(k1, ds1, ch1);


        // Test #2
        String[] ds2 = {"68698710.41", "59218262.54", "62565581.73", "47032202.15", "61616336.29", "56098619.96", "36743399.58", "30306423.09", "36355517.58", "48712490.14", "35163265.44", "37418799.37", "34950672.85", "66081810.06", "116520252.7", "121205816.5", "45898493.49", "42068899.07", "22026856.89", "28658828.25", "33792337.14", "33792337.14", "24941029.1"};
        String[] ch2 = {"-4P", "-4G", "-4A", "-4C", "-4S", "-4T", "-4V", "-4I", "-4L", "-4M", "-4F", "-4Y", "-4W", "-4H", "-4K", "-4R", "-4Q", "-4N", "-4D", "-4E", "-4s", "-4t", "-4y"};
        String k2 = "AMPKA2";
        boolean t2Result = testProbabilityMap(k2, ds2, ch2);

        String[] ds3 = {"22272478.91", "22500845.89", "23558498.69", "39058983.9", "31924743.39", "22374007.03", "18700629.54", "25303509.11", "29550644.39", "30948918.73", "28096580.81", "28949953.72", "26184388.76", "19668318.36", "31610796.05", "39965803.53", "19090456.7", "13983687.19", "8335811.64", "8296149.21", "14020833.18", "14020833.18", "13866540.32", "42088765.77", "30904902.53", "30350282.03", "33106508.08", "39842313.28", "26209596.25", "18002039.71", "13854103.2", "15388927.07", "23394112.75", "18224990.42", "22420527.45", "25903390.43", "34608218.6", "54081968", "56703233.63", "24234077.15", "19975765.79", "9733746.02", "11564934.49", "16580845.36", "16580845.36", "14164545.25", "21467121.04", "21963389.94", "25893977.65", "29837460.38", "49474787.16", "30280431.01", "12647389.21", "12164409.65", "13345903.01", "18408992.19", "22740430.1", "22213527.31", "27687057.22", "40771789.95", "84756425.69", "142560226.8", "26905817.91", "21287128.03", "8244859.96", "8028130.26", "12518275.78", "12518275.78", "11743434.71", "20903925.45", "33782290.03", "48038309.79", "52522536.6", "190144751.5", "63795112.22", "24542842.77", "14828650.07", "14170442.11", "16663972.35", "29807930.71", "28450238.44", "33573715.81", "26247524.8", "26190910.56", "33803634.62", "17590499.15", "21173629.12", "20476731.83", "16539395.64", "11323779.77", "11323779.77", "9190740.12", "40452739.46", "33889020.98", "28066678.56", "32507426.9", "39627954.92", "25294252.49", "20427810.5", "17625407.37", "25014455.16", "29818188.69", "22675182.54", "26320403.39", "16791624.22", "30788306.85", "39957210.24", "31500003.96", "29219893.72", "26737505.61", "14105913.62", "11335393.23", "8685699.96", "8685699.96", "23515867.32", "4900676.66", "12818533.65", "15811649.38", "65106961.1", "40667017.35", "29362700.29", "22116658.58", "24824815.76", "28331657.98", "72963006.22", "48876596.54", "28956593.73", "22374295.88", "19416966.4", "15089883.36", "25887298.58", "28869249.4", "22383344.76", "8248332.11", "9581577.43", "21994885.36", "21994885.36", "10152480.27", "13772369.72", "13068641.21", "16016238.22", "31575032.32", "73429518.06", "47680735.32", "15972741.82", "19125396.39", "20806524.2", "26772403.26", "47685810.26", "101817902.2", "55677663.14", "42038411.39", "9475811.07", "11709944.53", "28346003.1", "29022475.41", "14220304.77", "13252781.94", "15199950.46", "15199950.46", "14209262.35", "12478112.15", "26354964.77", "22507607.39", "59541676.45", "29187272.28", "18848017.6", "15459914.59", "16052048.57", "18877699.82", "20925849.28", "60331289.32", "33393901.78", "33627453.15", "38755319.3", "24026074.13", "28375060.26", "33133668.24", "83061252.63", "41677708.8", "17790029.15", "11708948.67", "11708948.67", "17683807.38", "14686701.54", "12615845.22", "13416051.92", "21063139.5", "12304131.4", "12447586.63", "19666072.84", "72417183.13", "68706937.82", "33033860", "29104620.49", "12551069.36", "13615920.47", "16122204.93", "16319167.3", "13442694.21", "14839748.84", "15052284.81", "9207766.28", "6714313.85", "10146132.4", "10146132.4", "6439578.58"};
        String[] ch3 = {"-5P", "-5G", "-5A", "-5C", "-5S", "-5T", "-5V", "-5I", "-5L", "-5M", "-5F", "-5Y", "-5W", "-5H", "-5K", "-5R", "-5Q", "-5N", "-5D", "-5E", "-5s", "-5t", "-5y", "-4P", "-4G", "-4A", "-4C", "-4S", "-4T", "-4V", "-4I", "-4L", "-4M", "-4F", "-4Y", "-4W", "-4H", "-4K", "-4R", "-4Q", "-4N", "-4D", "-4E", "-4s", "-4t", "-4y", "-3P", "-3G", "-3A", "-3C", "-3S", "-3T", "-3V", "-3I", "-3L", "-3M", "-3F", "-3Y", "-3W", "-3H", "-3K", "-3R", "-3Q", "-3N", "-3D", "-3E", "-3s", "-3t", "-3y", "-2P", "-2G", "-2A", "-2C", "-2S", "-2T", "-2V", "-2I", "-2L", "-2M", "-2F", "-2Y", "-2W", "-2H", "-2K", "-2R", "-2Q", "-2N", "-2D", "-2E", "-2s", "-2t", "-2y", "-1P", "-1G", "-1A", "-1C", "-1S", "-1T", "-1V", "-1I", "-1L", "-1M", "-1F", "-1Y", "-1W", "-1H", "-1K", "-1R", "-1Q", "-1N", "-1D", "-1E", "-1s", "-1t", "-1y", "1P", "1G", "1A", "1C", "1S", "1T", "1V", "1I", "1L", "1M", "1F", "1Y", "1W", "1H", "1K", "1R", "1Q", "1N", "1D", "1E", "1s", "1t", "1y", "2P", "2G", "2A", "2C", "2S", "2T", "2V", "2I", "2L", "2M", "2F", "2Y", "2W", "2H", "2K", "2R", "2Q", "2N", "2D", "2E", "2s", "2t", "2y", "3P", "3G", "3A", "3C", "3S", "3T", "3V", "3I", "3L", "3M", "3F", "3Y", "3W", "3H", "3K", "3R", "3Q", "3N", "3D", "3E", "3s", "3t", "3y", "4P", "4G", "4A", "4C", "4S", "4T", "4V", "4I", "4L", "4M", "4F", "4Y", "4W", "4H", "4K", "4R", "4Q", "4N", "4D", "4E", "4s", "4t", "4y"};
        String k3 = "BRSK2";
        boolean t3Result = testProbabilityMap(k3, ds3, ch3);

        return t1Result && t2Result;

    }

    private int findMax(int[] arr){
        int max = arr[0];
        for(int i = 0; i < arr.length; i++){
            if(arr[i] > max){
                max = arr[i];
            }
        }
        return max;
    }

    private int findMin(int[] arr){
        int min = arr[0];
        for(int i = 0; i < arr.length; i++){
            if(arr[i] < min){
                min = arr[i];
            }
        }
        return min;
    }



    private void splitColumnHeaders(String[] firstLine, int[] sites, String[] aminoAcids){
        for(int i = 0; i < firstLine.length; i++){
            int startIndex = 0;
            boolean foundDivider = false;
            while(startIndex < firstLine[i].length() && !foundDivider)
                if(firstLine[i].charAt(startIndex) == '-' || Character.isDigit(firstLine[i].charAt(startIndex))){
                    startIndex++;
                }
                else{
                    int siteValue = Integer.parseInt(firstLine[i].substring(0, startIndex));
                    String aminoAcidValue = firstLine[i].substring(startIndex);
                    sites[i] = siteValue;
                    aminoAcids[i] = aminoAcidValue;
                    foundDivider = true;
                }

        }
    }
}
