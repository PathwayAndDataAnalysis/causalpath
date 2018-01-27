package org.panda.causalpath.run;

import org.panda.causalpath.analyzer.CausalitySearcher;
import org.panda.causalpath.analyzer.ThresholdDetector;
import org.panda.causalpath.data.ActivityData;
import org.panda.causalpath.data.ProteinData;
import org.panda.causalpath.network.GraphWriter;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.resource.NetworkLoader;
import org.panda.causalpath.resource.ProteomicsFileReader;
import org.panda.causalpath.resource.ProteomicsLoader;
import org.panda.resource.ResourceDirectory;
import org.panda.resource.siteeffect.SiteEffectCollective;
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
	 * @param doSiteMatch option to enforce matching a phosphorylation site in the network with
	 *     the annotation of antibody
	 * @param siteMatchProximityThreshold when site matching is on, this parameter sets the proxomity threshold for a
	 *     site number in the relation to match the site whose change is observed in the data
	 * @param siteEffectProximityThreshold when the site effect is not known, we can approximate it with the known
	 *    effect of proximate sites. This parameter sets the proximity threshold for using the proximate sites for that
	 *    prediction.
	 * @param geneCentric Option to produce a gene-centric or an antibody-centric graph
	 * @param colorSaturationValue The value that maps to the most saturated color
	 * @param outputFilePrefix If the user provides xxx, then xxx.sif and xxx.format are generated
	 * @param customNetworkDirectory The directory that the network will be downloaded and SignedPC
	 *                               directory will be created in. Pass null to use default
	 * @throws IOException
	 */
	public static void generateCausalityGraph(String platformFile, String idColumn,
		String symbolsColumn, String sitesColumn, String effectColumn, String valuesFile,
		String valueColumn, double valueThreshold, String graphType, boolean doSiteMatch,
		int siteMatchProximityThreshold, int siteEffectProximityThreshold, boolean geneCentric,
		double colorSaturationValue, String outputFilePrefix, String customNetworkDirectory) throws IOException
	{
		if (customNetworkDirectory != null) ResourceDirectory.set(customNetworkDirectory);

		// Read platform file
		List<ProteomicsFileRow> rows = ProteomicsFileReader.readAnnotation(platformFile, idColumn, symbolsColumn,
			sitesColumn, effectColumn);

		// Read values
		List<String> vals = Collections.singletonList(valueColumn);
		ProteomicsFileReader.addValues(rows, valuesFile, idColumn, vals, 0D, false);

		// Fill-in missing effect
		SiteEffectCollective sec = new SiteEffectCollective();
		sec.fillInMissingEffect(rows, siteEffectProximityThreshold);

		generateCausalityGraph(rows, valueThreshold, graphType, doSiteMatch, siteMatchProximityThreshold,
			geneCentric, colorSaturationValue, outputFilePrefix);
	}

	/**
	 * For the given proteomics data, generates a ChiBE SIF graph.
	 *
	 * @param rows The proteomics data rows that are read from an external source
	 * @param valueThreshold The value threshold to be considered as significant
	 * @param graphType Either "compatible" or "conflicting"
	 * @param doSiteMatch option to enforce matching a phosphorylation site in the network with
	 *                       the annotation of antibody
	 * @param geneCentric Option to produce a gene-centric or an antibody-centric graph
	 * @param colorSaturationValue The value that maps to the most saturated color
	 * @param outputFilePrefix If the user provides xxx, then xxx.sif and xxx.format are generated
	 * @throws IOException
	 */
	public static void generateCausalityGraph(Collection<ProteomicsFileRow> rows, double valueThreshold,
		String graphType, boolean doSiteMatch, int siteMatchProximityThreshold, boolean geneCentric,
		double colorSaturationValue, String outputFilePrefix) throws IOException
	{
		ProteomicsLoader loader = new ProteomicsLoader(rows, null);
		// Associate change detectors
		loader.associateChangeDetector(new ThresholdDetector(valueThreshold, ThresholdDetector.AveragingMethod.ARITHMETIC_MEAN), data -> data instanceof ProteinData);
		loader.associateChangeDetector(new ThresholdDetector(0.1, ThresholdDetector.AveragingMethod.ARITHMETIC_MEAN), data -> data instanceof ActivityData);

		// Load signed relations
		Set<Relation> relations = NetworkLoader.load();
		loader.decorateRelations(relations);

		// Prepare causality searcher
		CausalitySearcher cs = new CausalitySearcher(!graphType.toLowerCase().startsWith("conflict"));
		cs.setForceSiteMatching(doSiteMatch);
		cs.setSiteProximityThreshold(siteMatchProximityThreshold);

		// Search causal or conflicting relations
		Set<Relation> relDat =  cs.run(relations);

		GraphWriter writer = new GraphWriter(relDat);
		writer.setUseGeneBGForTotalProtein(true);
		writer.setColorSaturationValue(colorSaturationValue);

		// Generate output
		if (geneCentric) writer.writeSIFGeneCentric(outputFilePrefix);
		else writer.writeSIFDataCentric(outputFilePrefix, cs.getInferenceUnits());
	}

	// Test in class. Bad practice. Tsk tsk tsk
	public static void main(String[] args) throws IOException
	{
		generateCausalityGraph("/home/babur/Documents/Temp/temp/platform.txt", "ID1", "Symbols", "Sites",
			"Effect", "/home/babur/Documents/Temp/temp/values.txt", "change", 0.001,
			"conflicting", false, 0, 0, false, 10.0, "/home/babur/Documents/Temp/temp/out", null);
	}
}
