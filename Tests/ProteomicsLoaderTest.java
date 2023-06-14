import org.junit.Before;
import org.junit.Test;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.ProteinSite;
import org.panda.causalpath.data.SiteModProteinData;
import org.panda.causalpath.resource.ProteomicsLoader;
import org.panda.resource.UniProtSequence;
import org.panda.resource.siteeffect.Feature;

import java.util.*;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * The purpose of the test class is to test the getSequence method.
 * It acomplishes this by passing a dummy data map to a ProteomicLoader
 * object's constructor, and then checking to see if getSequence() method
 * will recover the appropriate sequences and return them in a hashmap.
 */
public class ProteomicsLoaderTest {


    Map<String, Set<ExperimentData>> dataMapDummy;
    ProteomicsLoader pL;

    public static final String humanID = "9606";

    /**
     * SetUp method will set up a dummyDataMap object, containing
     * several different SiteModProteinData Objects. There are four
     * main cases that were to be tested; SiteModProteinData objects containing sitemaps that have:
     * single prot single site(valid),
     * single prot multiple site(invalid), multiple prot each with one site(valid),
     * multiple prot multiple sites(invalid). SiteModProtein objects that meet these
     * cases were created, then added into the datamap.
     *
     * For each SiteModProteinData object, we set the changeDetector to a DummyChangeDetector
     * so we could assign it a known changeValue (the DummyChangeDetector allows the setting of the changeValue,
     * whereas other ChangeDetector implementations do not)
     * The purpose of this was so that we could easily verify that for SiteModProteinData objects that
     * meet the valid conditions, these object's associated changeValues could be found in the hashMap returned
     * by the getSequence method.
     */
    @Before
    public void setUp(){

        dataMapDummy = new HashMap<>();


        // Creating SiteModProteinData objects

        // SiteModProteinData with SiteMap containing protein STAT3 (Single Protein, Single Site)
        ExperimentData item1 = createSiteModProteinData(-2.445136194872271, proteinNamesSet("STAT3"), true, toSet(705));
        // SiteModProteinData with SiteMap containing protein MAP2K1 (Single Protein, Multiple Site)
        ExperimentData item2 = createSiteModProteinData(-0.27898862494155197, proteinNamesSet("MAP2K1"), true, toSet(218, 217, 221, 222) );
        // SiteModProteinData (Multiple Protein, Multiple Site)
        ExperimentData item3 = createSiteModProteinData(-0.619769844905208, proteinNamesSet("MAPK1", "MAPK3") , true, toSet(185, 187), toSet(204, 202));
        // SiteModProteinData with Multiple Protein Single Site
        ExperimentData item4 = createSiteModProteinData(7.619689476352495, proteinNamesSet("AKT2", "AKT1", "AKT3"), true, toSet(474), toSet(472), toSet(473));
        // SiteModProteinData containing protein YAP1 (Single Protein, Single Site)
        ExperimentData item5 = createSiteModProteinData(8.024169754065225, proteinNamesSet("YAP1"), true, toSet(127));
        // SiteModProteinData object that is not of type Phosphorylation (check that program only picks out SMPD objects with mod Phosphorylation)
        ExperimentData item6 = createSiteModProteinData(0, proteinNamesSet("dummyProtein"), false, toSet(50));

        // Input into dataMap
        dataMapDummy.put("geneA",  toSet(item1));
        dataMapDummy.put("geneB", toSet(item2));
        dataMapDummy.put("geneC", toSet(item3));
        dataMapDummy.put("geneD", toSet(item4));
        dataMapDummy.put("geneE", toSet(item5, item6));



        // ProteomicLoader object with dummy constructor passed dataMap
        pL = new ProteomicsLoader(dataMapDummy);
    }


    /**
     * Method takes a variable number of arguments of type t and returns
     * a set. Useful for easy creation of sets without repeatedly calling add
     * @param items
     * @return
     * @param <T>
     */
    public <T> Set<T> toSet(T... items){
        HashSet<T> s = new HashSet<>(Arrays.asList(items));
        Collections.addAll(s, items);
        return s;
    }

    /**
     * Method takes a variable number of proteinNames as strings and returns
     * a list containing all of them. Useful for quickly creating a datastructure
     * containing several proteinNames.
     * @param proteinNames
     * @return
     */
    public List<String> proteinNamesSet(String... proteinNames){
        ArrayList<String> s = new ArrayList<String>();
        Collections.addAll(s, proteinNames);

        return s;
    }

    /**
     * Constructs a siteMap object with the arguments passed.
     * @param proteinNames The names of the proteins(keys) of this sitemap. Note that the order of the list matters.
     * @param arrayOfSiteSets An array where each element is a set containing integers phosphorylated locations
     * @return
     */
    public Map<String, Set<ProteinSite>> createSiteMap(List<String> proteinNames, Set<Integer>[] arrayOfSiteSets){
        // Create a HashMap to store the siteMap
        Map<String, Set<ProteinSite>> siteMap = new HashMap<>();

        // Iterate over the proteins
        for(int i = 0; i < proteinNames.size(); i++){
            String currProtein = proteinNames.get(i);

            HashSet<ProteinSite> pSite = new HashSet<>();

            // Protein names and listofSiteSets correspond to eachother sequentially
            Set<Integer> sitesForProtein = arrayOfSiteSets[i];

            for(int j: sitesForProtein){
                pSite.add(new ProteinSite(j, "n/a", 1));
            }

            siteMap.put(currProtein, pSite);
        }

        return siteMap;
    }

    public String getSequence(String name, int location){
        UniProtSequence uI = UniProtSequence.get();
        name = uI.getNameOfSymbol(name, humanID);
        return uI.getSeqAround(name, 5, 4, location);
    }

    /**
     * Method responsible for creating SiteModProteinData objects appropriately
     * @param changeValue changeValue for the object, so that we can check to verify if this exists in the HashMap returned by getSequence
     * @param proteinNames names for the proteins, where order matters(they correspond via index with site sets in arrayOfSiteSets)
     * @param isPhosphorylation Whether this object's feature is Phosphorylation or not
     * @param arrayOfSiteSets The siteSet for each given protein in proteinNames
     * @return The SiteModProteinData object instantiated based off the arguments passed
     */
    public ExperimentData createSiteModProteinData(double changeValue, List<String> proteinNames, boolean isPhosphorylation, Set<Integer> ... arrayOfSiteSets){
        SiteModProteinData s;
        if(isPhosphorylation){
            s = new SiteModProteinData("empty", new HashSet<String>(), Feature.PHOSPHORYLATION);
            s.setChDet(new DummyChangeDetector(changeValue));
        }
        else{
            // This will help test to make sure that the method won't wrongly pick up the wrong type of ExperimentData
            // As this will break the code
            s = new SiteModProteinData("differentID", new HashSet<String>(), Feature.ACETYLATION);
        }

        // Sets up the SiteModProteinData object's sitemap appropriately
        s.setSiteMap(createSiteMap(proteinNames, arrayOfSiteSets));

        return s;
    }


    @Test
    public void testGetSequences(){

        // Test that single protein multiple locations is not put into map
        assertFalse(pL.getSeqChangeValMap().containsValue(-0.27898862494155197));
        // Test that multiple protein multiple site is not put into map
        assertFalse(pL.getSeqChangeValMap().containsValue(-0.619769844905208));
        // Test that multiple protein single site where aa sequence DOES NOT exactly match IS NOT put into map
        assertFalse(pL.getSeqChangeValMap().containsValue(7.619689476352495));

        // Test that program inputs sequence from single protein single site case into map
        assertTrue(pL.getSeqChangeValMap().containsKey((getSequence("STAT3", 705))));

        // Check that the corresponding changeValue is correct
        assertTrue( (pL.getSeqChangeValMap().get(getSequence("STAT3", 705))) == -2.445136194872271);






        // Upon debugging, it appears that UniProt object returning null for getNameFromSymbol. Working on this,
        // believe issue has to do with the FTP file
        assertTrue(pL.getSeqChangeValMap().containsKey((getSequence("YAP1", 127))));

        assertTrue((pL.getSeqChangeValMap().get(getSequence("YAP1", 127))) == 8.024169754065225);

        assertTrue(pL.getSeqChangeValMap().size() == 2);
    }

}
