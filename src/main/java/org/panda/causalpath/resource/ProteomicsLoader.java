package org.panda.causalpath.resource;

import org.panda.causalpath.analyzer.CausalityHelper;
import org.panda.causalpath.analyzer.OneDataChangeDetector;
import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;
import org.panda.resource.UniProtSequence;
import org.panda.resource.siteeffect.Feature;
import org.panda.resource.tcga.ProteomicsFileRow;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Coverts the RPPAData in resource project to appropriate objects for this project and serves them.
 */
public class ProteomicsLoader {
    /**
     * Map from genes to related data.
     */
    Map<String, Set<ExperimentData>> dataMap;

    /**
     * Set of all data. This collection supposed to hold everything in the dataMap's values.
     */
    Set<ExperimentData> datas;

    /**
     * Uniprot organism ID for humans
     */
    public static final String humanID = "9606";

    public ProteomicsLoader(Collection<ProteomicsFileRow> rows, Map<DataType, Double> stdevThresholds) {

        dataMap = new HashMap<>();
        datas = new HashSet<>();
        rows.stream().distinct().forEach(r ->
        {
            ExperimentData ed = r.isActivity() ? new ActivityData(r) :
                    r.isSiteSpecific() ? new SiteModProteinData(r) :
                            r.isRNA() ? new RNAData(r) :
                                    r.isMetabolite() ? new MetaboliteData(r) :
                                            new ProteinData(r);

            if (stdevThresholds != null && ed instanceof NumericData) {
                DataType type = ed.getType();
                Double thr = stdevThresholds.get(type);
                if (thr != null) {
                    double sd = Summary.stdev(((NumericData) ed).vals);
                    if (Double.isNaN(sd) || sd < thr) return;
                }
            }

            for (String sym : ed.getGeneSymbols()) {
                if (!dataMap.containsKey(sym)) dataMap.put(sym, new HashSet<>());

                // check if there is already some data with the same ID
                for (ExperimentData data : dataMap.get(sym)) {
                    if (data.getId().equals(ed.getId())) {
                        throw new RuntimeException("Proteomic data has non-unique IDs: " + ed.getId());
                    }
                }

                dataMap.get(sym).add(ed);
                datas.add(ed);
            }
        });
    }

    /**
     * Constructor to be used for testing purposes only.
     *
     * @param dataMap dummy dataMap object for testing of getSequences method
     */
    public ProteomicsLoader(Map<String, Set<ExperimentData>> dataMap) {
        this.dataMap = dataMap;
    }


    /**
     * Looks through dataMap and obtains all sequences corresponding to proteinSiteModData
     * objects (with feature Phosphorylation) whose sitemaps have A) a single protein with
     * a single site from which it is possible to derive a sequence or B) multiple proteins
     * with a single site for each protein, for which the sequences extracted for each protein-site
     * pair in the siteMap is the same
     *
     * @return Returns a map of all sequences and their corresponding changeValues
     */
    public HashMap<String, Double> getSeqChangeValMap() {

        // Set to store all obtained sequences
        HashMap<String, Double> seqChangeVal = new HashMap<>();

        for (String gene : dataMap.keySet()) {
            // Get the corresponding value in the map for this gene
            Set<ExperimentData> geneData = dataMap.get(gene);

            for (ExperimentData data : geneData) {
                // This line will identify if an experimentData object is a SiteModProteinData object
                // and it will also verify that the feature of this object is phosphorylation
                if (isPhosphorylationData(data)) {
                    SiteModProteinData siteModData = (SiteModProteinData) data;
                    Map<String, Set<ProteinSite>> siteMap = siteModData.getSiteMap();

					/* Two cases will be handled separately (via if and else block)
					   if block will handle multiple proteins within sitemap case;
					   else handles single prot in sitemap
					 */
                    if (siteMap.size() > 1) {

                        // We will only consider this siteMap if each protein has exactly one site
                        if ((validNumberProtSites(siteMap))) {


						/* Call method that will check if every protein in siteMap has the same sequence around its corresponding
						   site. Method call will return null should the sequence extracted be different
						   for the different proteins
						 */

                            String sequence = getSequenceMultiProtein(siteMap);

                            if (sequence != null) {
                                seqChangeVal.put(sequence, siteModData.getChangeValue());
                            }
                        }


                    } else {
                        // Handles case where it is a single protein
                        for (String protein : siteMap.keySet()) {
                            Set<ProteinSite> proteinSiteSet = siteMap.get(protein);
                            // If we identify a protein that has more than one phosphorylated site,
                            // we will not include it in the set of sequences
                            if (proteinSiteSet.size() == 1) {
                                for (ProteinSite pS : proteinSiteSet) {
                                    // Identify the location
                                    int location = pS.getSite();
                                    // Get the uniprot name of the protein
                                    String uniprotName = UniProtSequence.get().getNameOfSymbol(protein, humanID);
                                    String seqAroundSite = UniProtSequence.get().getSeqAround(uniprotName, 5, 4, location);
                                    if (seqAroundSite != null) {
                                        seqChangeVal.put(seqAroundSite, siteModData.getChangeValue());
                                    }
                                }
                            }

                        }
                    }
                }

            }

        }
        return seqChangeVal;
    }


    /**
     * Method sorts the sequences in sequenceChangeValueMap by their
     * changeValue, and then returns the sorted list
     *
     * @return list of sequences in ascending order based on their change value
     */
    public List<String> getSortedSequences() {
        List<String> sortedKeys = getSeqChangeValMap().entrySet()
                .stream()
                .sorted(Map.Entry.comparingByValue())
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());

        return sortedKeys;
    }


    /**
     * This method will iterate over a sitemap which contains multiple proteins
     * but single sites for each protein. It will check to make sure that every protein - site pair
     * in the map shares the same sequence. If there ar any differences, it will return null,
     * else if the sequences are the same it will return true.
     *
     * @param siteMap
     * @return sequence as a string; or null if the sequence is not common between proteins
     */
    private String getSequenceMultiProtein(Map<String, Set<ProteinSite>> siteMap) {

        String sequence = null;

        // UniProtSequence object to obtain sequences

        for (String protein : siteMap.keySet()) {
            Set<ProteinSite> phosphorylatedSites = siteMap.get(protein);
            for (ProteinSite proteinSite : phosphorylatedSites) {

                int location = proteinSite.getSite();
                // Obtain the UNIPROT name for this protein
                String nameOrID = UniProtSequence.get().getNameOfSymbol(protein, humanID);

                // If statement runs prior to initial assignment to sequence variable
                // I.e. First iteration of the for loop
                if (sequence == null) {
                    sequence = UniProtSequence.get().getSeqAround(nameOrID, 5, 4, location);
					/* Method call may return null if there is it is not possible to recover
					 sequence of such length. If the if statement runs,
					 it cannot be true that every protein shares the same sequence surrounding location
					 as for current protein we cannot obtain a sequence. Thus return null.
					 */
                    if (sequence == null) {
                        return null;
                    }
                } else {
                    // Else runs for all iterations after first iteration

                    // currSeq refers to the sequence surrounding the location for the current protein
                    String currSeq = UniProtSequence.get().getSeqAround(nameOrID, 5, 4, location);

                    // currSeq may be null in case where it is not possible to recover sequence of desired length
                    // In this case, it is not true that all the proteins share the same sequence, thus ret null
                    if (currSeq == null) {
                        return null;
                    }
                    // If currSeq does not match prior sequence, then it is not true that all proteins share
                    // same seq, thus ret null
                    else if (!currSeq.equals(sequence)) {
                        return null;
                    }

                }
            }

        }

        return sequence;
    }

    /**
     * Method that verfies the experimental data is of type SiteModProteinData
     * with modification phosphorylation
     *
     * @param data ExperimentalData object to verify
     * @return true if the ExperimentalData object is of type SiteModProteinData with modification Phosphorylation, else false
     */
    private boolean isPhosphorylationData(ExperimentData data) {
        if (data instanceof SiteModProteinData) {
            SiteModProteinData siteModData = (SiteModProteinData) data;
            return siteModData.getModification().equals(Feature.PHOSPHORYLATION);
        }

        return false;
    }

    /**
     * Method that iterates over a siteMap and verifies that each protein in the
     * siteMap has a corresponding set of protein sites that is exactly 1 in size
     *
     * @param siteMap
     * @return true if each protein has ProteinSite list of length 1, false otherwise
     */
    private boolean validNumberProtSites(Map<String, Set<ProteinSite>> siteMap) {
        for (String protein : siteMap.keySet()) {
            // If a protein exists with multiple sites on it that have been phosphorylated
            // it is invalid
            if (siteMap.get(protein).size() > 1) {
                return false;
            }
        }

        return true;
    }


    public void addRepeatData(Collection<ProteomicsFileRow> rows, Map<DataType, Double> stdevThresholds) {
        rows.stream().distinct().forEach(r ->
        {
            ExperimentData ed = r.isActivity() ? new ActivityData(r) :
                    r.isSiteSpecific() ? new SiteModProteinData(r) : new ProteinData(r);

            if (stdevThresholds != null && ed instanceof NumericData) {
                DataType type = ed.getType();
                Double thr = stdevThresholds.get(type);
                if (thr != null) {
                    double sd = Summary.stdev(((NumericData) ed).vals);
                    if (Double.isNaN(sd) || sd < thr) return;
                }
            }

            for (String sym : ed.getGeneSymbols()) {
                if (!dataMap.containsKey(sym)) dataMap.put(sym, new HashSet<>());

                ExperimentData orig = null;
                for (ExperimentData data : dataMap.get(sym)) {
                    if (data.getId().equals(ed.getId())) {
                        orig = data;
                        break;
                    }
                }

                if (orig != null) {
                    orig.addRepeatData(ed);
                } else {
                    dataMap.get(sym).add(ed);
                    datas.add(ed);
                }
            }
        });
    }

    public void printStDevHistograms() {
        printStDevHistograms(datas);
    }

    public void printStDevHistograms(Set<ExperimentData> datas) {
        System.out.println("\nSt-dev histograms:");
        Map<DataType, Histogram> hMap = new HashMap<>();
        datas.stream().filter(d -> d instanceof NumericData).map(d -> (NumericData) d).forEach(d ->
        {
            DataType type = d.getType();
            if (!hMap.containsKey(type)) {
                Histogram h = new Histogram(0.05);
                h.setBorderAtZero(true);
                hMap.put(type, h);
            }
            hMap.get(type).count(Summary.stdev(d.vals));
        });
        hMap.keySet().forEach(k ->
        {
            System.out.println("type = " + k);
            hMap.get(k).print();
        });
    }

    /**
     * Adds the related data to the given relations.
     */
    public void decorateRelations(Set<Relation> relations) {
        Map<String, GeneWithData> map = collectExistingData(relations);
        CausalityHelper ch = new CausalityHelper();
        for (Relation rel : relations) {
            if (rel.sourceData == null) {
                if (map.containsKey(rel.source)) {
                    rel.sourceData = map.get(rel.source);
                } else {
                    GeneWithData gwd = new GeneWithData(rel.source);
                    gwd.addAll(dataMap.get(rel.source));
                    map.put(gwd.getId(), gwd);
                }
            } else rel.sourceData.addAll(dataMap.get(rel.source));

            if (rel.targetData == null) {
                if (map.containsKey(rel.target)) {
                    rel.targetData = map.get(rel.target);
                } else {
                    GeneWithData gwd = new GeneWithData(rel.target);
                    gwd.addAll(dataMap.get(rel.target));
                    map.put(gwd.getId(), gwd);
                }
            } else rel.targetData.addAll(dataMap.get(rel.target));

            rel.chDet = ch;
        }
    }

    private Map<String, GeneWithData> collectExistingData(Set<Relation> relations) {
        Map<String, GeneWithData> map = new HashMap<>();

        for (Relation rel : relations) {
            if (rel.sourceData != null) map.put(rel.sourceData.getId(), rel.sourceData);
            if (rel.targetData != null) map.put(rel.targetData.getId(), rel.targetData);
        }

        return map;
    }

    /**
     * Puts the given change detector to the data that is filtered by the given selector.
     */
    public void associateChangeDetector(OneDataChangeDetector chDet, DataSelector selector) {
        datas.stream().filter(selector::select).forEach(d -> d.setChDet(chDet));
    }

    public void initMissingDataForProteins() {
        Optional<ProteinData> opt = datas.stream().filter(d -> d instanceof ProteinData)
                .map(d -> (ProteinData) d).findAny();

        if (!opt.isPresent()) return;

        int size = opt.get().vals.length;

//		int[] totalProtCnt = new int[size];
//		int[] phospProtCnt = new int[size];

        boolean[] hasTotalProt = new boolean[size];
        boolean[] hasPhospProt = new boolean[size];

        Arrays.fill(hasTotalProt, false);
        Arrays.fill(hasPhospProt, false);

        datas.stream().filter(d -> d instanceof ProteinData).map(d -> (ProteinData) d).forEach(d ->
        {
            if (d.getType() == DataType.PROTEIN) {
                for (int i = 0; i < size; i++) {
                    if (!Double.isNaN(d.vals[i])) {
                        hasTotalProt[i] = true;
//						totalProtCnt[i]++;
                    }
                }
            } else if (d.getType() == DataType.PHOSPHOPROTEIN) {
                for (int i = 0; i < size; i++) {
                    if (!Double.isNaN(d.vals[i])) {
                        hasPhospProt[i] = true;
//						phospProtCnt[i]++;
                    }
                }
            }
        });

        datas.stream().filter(d -> d instanceof ProteinData).map(d -> (ProteinData) d)
                .forEach(d -> d.initPresenceData(d.getType() == DataType.PHOSPHOPROTEIN ? hasPhospProt : hasTotalProt));

//		System.out.println("Total protein counts");
//		Histogram h = new Histogram(100, ArrayUtil.toDouble(totalProtCnt));
//		h.print();
//
//		System.out.println("\nPhosphoprotein counts");
//		h = new Histogram(100, ArrayUtil.toDouble(phospProtCnt));
//		h.print();
    }

    /**
     * Function to filter experiment data.
     */
    public interface DataSelector {
        boolean select(ExperimentData data);
    }
}
