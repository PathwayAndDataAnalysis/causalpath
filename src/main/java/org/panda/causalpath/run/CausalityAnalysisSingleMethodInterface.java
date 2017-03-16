package org.panda.causalpath.run;

import org.panda.causalpath.analyzer.CausalitySearcher;
import org.panda.causalpath.analyzer.RelationTargetCompatibilityChecker;
import org.panda.causalpath.analyzer.ThresholdDetector;
import org.panda.causalpath.data.ActivityData;
import org.panda.causalpath.data.ProteinData;
import org.panda.causalpath.network.GraphWriter;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationAndSelectedData;
import org.panda.causalpath.resource.ProteomicsFileReader;
import org.panda.causalpath.resource.ProteomicsLoader;
import org.panda.causalpath.resource.NetworkLoader;
import org.panda.resource.PhosphoSitePlus;
import org.panda.resource.ResourceDirectory;
import org.panda.resource.tcga.ProteomicsFileRow;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * This class provides a single basic method interface to the causality analysis, which is specifically designed to be
 * used by tools in other programming languages, such as R.
 *
 * @author Ozgun Babur
 */
public class CausalityAnalysisSingleMethodInterface
{
	/**
	 * Reads the proteomics platform and data files, and generates a ChiBE SIF graph.
	 *
	 * @param platformFile Name of the antibody reference file
	 * @param idColumn Column name of IDs
	 * @param symbolsColumn Column name of gene symbols
	 * @param sitesColumn Column name for phosphorylation sites
	 * @param effectColumn Column name for effect of the site on activity
	 * @param valuesFile Name of the measurements file
	 * @param valueColumn Name of the values column in the measurements file
	 * @param valueThreshold The value threshold to be considered as significant
	 * @param graphType Either "compatible" or "conflicting"
	 * @param siteMatchStrict option to enforce matching a phosphorylation site in the network with
	 *     the annotation of antibody
	 * @param siteMatchProximityThreshold when site matching is on, this parameter sets the proxomity threshold for a
	 *     site number in the relation to match the site whose change is observed in the data
	 * @param siteEffectProximityThreshold when the site effect is not known, we can approximate it with the known
	 *    effect of proximate sites. This parameter sets the proximity threshold for using the proximate sites for that
	 *    prediction.
	 * @param geneCentric Option to produce a gene-centric or an antibody-centric graph
	 * @param outputFilePrefix If the user provides xxx, then xxx.sif and xxx.format are generated
	 * @param customNetworkDirectory The directory that the network will be downloaded and SignedPC
	 *                               directory will be created in. Pass null to use default
	 * @throws IOException
	 */
	public static void generateCausalityGraph(String platformFile, String idColumn,
		String symbolsColumn, String sitesColumn, String effectColumn, String valuesFile,
		String valueColumn, double valueThreshold, String graphType, boolean siteMatchStrict,
		int siteMatchProximityThreshold, int siteEffectProximityThreshold, boolean geneCentric,
		String outputFilePrefix, String customNetworkDirectory) throws IOException
	{
		if (customNetworkDirectory != null) ResourceDirectory.set(customNetworkDirectory);

		// Read platform file
		List<ProteomicsFileRow> rows = ProteomicsFileReader.readAnnotation(platformFile, idColumn, symbolsColumn,
			sitesColumn, effectColumn);

		// Read values
		List<String> vals = Collections.singletonList(valueColumn);
		ProteomicsFileReader.addValues(rows, valuesFile, idColumn, vals, 0D, false);

		// Fill-in missing effect from PhosphoSitePlus
		PhosphoSitePlus.get().fillInMissingEffect(rows, siteEffectProximityThreshold);

		generateCausalityGraph(rows, valueThreshold, graphType, siteMatchStrict, siteMatchProximityThreshold,
			geneCentric, outputFilePrefix);
	}

	/**
	 * For the given proteomics data, generates a ChiBE SIF graph.
	 *
	 * @param rows The proteomics data rows that are read from an external source
	 * @param valueThreshold The value threshold to be considered as significant
	 * @param graphType Either "compatible" or "conflicting"
	 * @param siteMatchStrict option to enforce matching a phosphorylation site in the network with
	 *                       the annotation of antibody
	 * @param geneCentric Option to produce a gene-centric or an antibody-centric graph
	 * @param outputFilePrefix If the user provides xxx, then xxx.sif and xxx.format are generated
	 * @throws IOException
	 */
	public static void generateCausalityGraph(Collection<ProteomicsFileRow> rows, double valueThreshold,
		String graphType, boolean siteMatchStrict, int siteMatchProximityThreshold, boolean geneCentric,
		String outputFilePrefix) throws IOException
	{
		ProteomicsLoader loader = new ProteomicsLoader(rows);
		// Associate change detectors
		loader.associateChangeDetector(new ThresholdDetector(valueThreshold, ThresholdDetector.AveragingMethod.ARITHMETIC_MEAN), data -> data instanceof ProteinData);
		loader.associateChangeDetector(new ThresholdDetector(0.1, ThresholdDetector.AveragingMethod.ARITHMETIC_MEAN), data -> data instanceof ActivityData);

		// Prepare relation-target compatibility checker
		RelationTargetCompatibilityChecker rtcc = new RelationTargetCompatibilityChecker();
		rtcc.setForceSiteMatching(siteMatchStrict);
		rtcc.setSiteProximityThreshold(siteMatchProximityThreshold);

		// Load signed relations
		Set<Relation> relations = NetworkLoader.load();
		loader.decorateRelations(relations, rtcc);

		// Prepare causality searcher
		CausalitySearcher cs = new CausalitySearcher(rtcc);
		if (graphType.toLowerCase().startsWith("conflict")) cs.setCausal(false);

		// Search causal or conflicting relations
		Set<RelationAndSelectedData> relDat =  cs.run(relations);

		GraphWriter writer = new GraphWriter(relDat);
		writer.setUseGeneBGForTotalProtein(true);

		// Generate output
		if (geneCentric) writer.writeSIFGeneCentric(outputFilePrefix);
		else writer.writeSIFDataCentric(outputFilePrefix);
	}

	// Test in class. Bad practice. Tsk tsk tsk
	public static void main(String[] args) throws IOException
	{
		generateCausalityGraph("/home/babur/Documents/Temp/temp/abdata-chibe.txt", "ID1", "Symbols", "Sites",
			"Effect", "/home/babur/Documents/Temp/temp/ovcar4_dif_drug_sig.txt", "change", 0.001,
			"compatible", true, 0, 0, false, "/home/babur/Documents/Temp/temp/out", null);
	}
}
