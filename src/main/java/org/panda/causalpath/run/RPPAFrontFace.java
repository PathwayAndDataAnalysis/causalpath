package org.panda.causalpath.run;

import org.panda.causalpath.analyzer.CausalitySearcher;
import org.panda.causalpath.analyzer.ThresholdDetector;
import org.panda.causalpath.data.Activity;
import org.panda.causalpath.data.ActivityData;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.ProteinData;
import org.panda.causalpath.network.GraphWriter;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationAndSelectedData;
import org.panda.causalpath.resource.RPPAFileReader;
import org.panda.causalpath.resource.RPPALoader;
import org.panda.causalpath.resource.SignedPCUser;
import org.panda.resource.PhosphoSitePlus;
import org.panda.resource.ResourceDirectory;
import org.panda.resource.tcga.RPPAData;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * This class reads the RPPA platform and data files, and generates a ChiBE SIF graph.
 *
 * @author Ozgun Babur
 */
public class RPPAFrontFace
{
	/**
	 * Reads the RPPA platform and data files, and generates a ChiBE SIF graph.
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
	 *                       the annotation of antibody
	 * @param geneCentric Option to produce a gene-centric or an antibody-centric graph
	 * @param outputFilePrefix If the user provides xxx, then xxx.sif and xxx.format are generated
	 * @param customNetworkDirectory The directory that the network will be downloaded and SignedPC
	 *                               directory will be created in. Pass null to use default
	 * @throws IOException
	 */
	public static void generateRPPAGraphs(String platformFile, String idColumn,
		String symbolsColumn, String sitesColumn, String effectColumn, String valuesFile,
		String valueColumn, double valueThreshold, String graphType, boolean siteMatchStrict,
		boolean geneCentric, String outputFilePrefix, String customNetworkDirectory)
		throws IOException
	{
		if (customNetworkDirectory != null) ResourceDirectory.set(customNetworkDirectory);

		// Read platform file
		List<RPPAData> rppas = RPPAFileReader.readAnnotation(platformFile, idColumn, symbolsColumn,
			sitesColumn, effectColumn);

		// Read values
		List<String> vals = Collections.singletonList(valueColumn);
		RPPAFileReader.addValues(rppas, valuesFile, idColumn, vals, 0D);

		// Fill-in missing effect from PhosphoSitePlus
		PhosphoSitePlus.get().fillInMissingEffect(rppas, 0);

		generateRPPAGraphs(rppas, valueThreshold, graphType, siteMatchStrict, geneCentric, outputFilePrefix);
	}

	/**
	 * For the given RPPA data, generates a ChiBE SIF graph.
	 *
	 * @param rppas The RPPA data that is read from en external source
	 * @param valueThreshold The value threshold to be considered as significant
	 * @param graphType Either "compatible" or "conflicting"
	 * @param siteMatchStrict option to enforce matching a phosphorylation site in the network with
	 *                       the annotation of antibody
	 * @param geneCentric Option to produce a gene-centric or an antibody-centric graph
	 * @param outputFilePrefix If the user provides xxx, then xxx.sif and xxx.format are generated
	 * @throws IOException
	 */
	public static void generateRPPAGraphs(Collection<RPPAData> rppas, double valueThreshold, String graphType,
		boolean siteMatchStrict, boolean geneCentric, String outputFilePrefix)
		throws IOException
	{
		RPPALoader loader = new RPPALoader(rppas);
		// Associate change detectors
		loader.associateChangeDetector(new ThresholdDetector(valueThreshold), data -> data instanceof ProteinData);
		loader.associateChangeDetector(new ThresholdDetector(0.1), data -> data instanceof ActivityData);

		// Load signed relations
		Set<Relation> relations = SignedPCUser.getSignedPCRelations();
		loader.decorateRelations(relations);

		// Prepare causality searcher
		CausalitySearcher cs = new CausalitySearcher();
		cs.setForceSiteMatching(siteMatchStrict);

		// Search causal or conflicting relations
		Set<RelationAndSelectedData> relDat = graphType.toLowerCase().startsWith("conflict") ?
			cs.selectConflictingRelations(relations) : cs.selectCausalRelations(relations);

		GraphWriter writer = new GraphWriter(relDat);
		writer.setUseGeneBGForTotalProtein(true);

		// Generate output
		if (geneCentric) writer.writeGeneCentric(outputFilePrefix);
		else writer.writeDataCentric(outputFilePrefix);
	}

	// Test in class. Bad practice. Tsk tsk tsk
	public static void main(String[] args) throws IOException
	{
		generateRPPAGraphs("/home/ozgun/Documents/JQ1/abdata-chibe.txt", "ID1", "Symbols", "Sites",
			"Effect", "/home/ozgun/Documents/JQ1/ovcar4_dif_drug_sig.txt", "change", 0.001,
			"compatible", true, false, "/home/ozgun/Temp/temp", null);
	}
}
