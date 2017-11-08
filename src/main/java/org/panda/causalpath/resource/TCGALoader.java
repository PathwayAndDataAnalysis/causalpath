package org.panda.causalpath.resource;

import org.panda.causalpath.analyzer.OneDataChangeDetector;
import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;
import org.panda.resource.tcga.*;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Loads experiment data from TCGA files that are downloadable from Broad Firehose.
 */
public class TCGALoader
{
	/**
	 * Reader for copy number data.
	 */
	private CNAReader cnaReader;

	/**
	 * Reader for RNA expression data.
	 */
	private ExpressionReader expReader;

	/**
	 * Reader for mutation data.
	 */
	private MutationReader mutReader;

	/**
	 * Reader for RPPA data.
	 */
	private RPPAReader rppaReader;

	// todo add methylation data

	/**
	 * RPPA reader is different from other readers. A row may correspond to multiple gene symbols. This cache is used
	 * to prevent duplicate creation of experiment data for the same row.
	 */
	private Map<String, ProteinData> rppaCache;

	/**
	 * When a data is already been read for a gene, this cache stores it so that creation of redundant data objects is
	 * avoided on multiple access.
	 */
	private Map<String, Set<ExperimentData>> dataCache;

	/**
	 * This set of genes are used to ignore RNA expression, CNA or methylation data when there is a total protein
	 * measurement exists.
	 *
	 */
	private Set<String> genesWithTotalProteinData;

	/**
	 * Names of samples loaded.
	 */
	private String[] samples;

	private Map<String, Integer> mutationEffectMap;

	// Names of data files under the TCGA data directory

	private static final String COPY_NUMBER_FILE = "/copynumber.txt";
	private static final String MUTATION_FILE = "/mutation.maf";
	private static final String EXPRESSION_FILE = "/expression.txt";
	private static final String RPPA_FILE = "/rppa.txt";

	// todo: Make loading parameterized by passing the types of desired data.
	public TCGALoader(String dir)
	{
		try{cnaReader = new CNAReader(dir + COPY_NUMBER_FILE, false, 0);} catch (FileNotFoundException e){}
		try{expReader = new ExpressionReader(dir + EXPRESSION_FILE);} catch (FileNotFoundException e){}
		try{mutReader = new MutationReader(dir + MUTATION_FILE);} catch (IOException e){}
		try{rppaReader = new RPPAReader(dir + RPPA_FILE);} catch (FileNotFoundException e){}
		this.samples = getUnionSamples();
		if (rppaReader != null) rppaCache = new HashMap<>();
		dataCache = new HashMap<>();
		genesWithTotalProteinData = findGenesWithTotalProteinData();
	}

	public void setSamples(String[] samples)
	{
		this.samples = samples;
	}

	public void loadMutationEffectMap(String filename) throws IOException
	{
		this.mutationEffectMap = Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> Integer.valueOf(t[1])));
	}

	public void setMutationEffectMap(Map<String, Integer> mutationEffectMap)
	{
		this.mutationEffectMap = mutationEffectMap;
	}

	/**
	 * Gets the set of experiment data for the given symbol.
	 */
	public Set<ExperimentData> getData(String symbol)
	{
		// if asked before, return from cache
		if (dataCache.containsKey(symbol)) return dataCache.get(symbol);

		Set<ExperimentData> set = new HashSet<>();

		if (rppaReader != null)
		{
			Set<ProteomicsFileRow> rowSet = rppaReader.getAssociatedData(symbol, samples);

			rowSet.stream().filter(row -> !rppaCache.containsKey(row.id)).forEach(row -> {
				ProteinData d = row.isPhospho() ? new PhosphoProteinData(row) : new ProteinData(row);
				rppaCache.put(d.id, d);
			});

			rowSet.stream().map(rppa -> rppaCache.get(rppa.id)).forEach(set::add);
		}

		if (cnaReader != null)
		{
			int[] val = cnaReader.getGeneAlterationArray(symbol, samples);
			if (val != null)
			{
				CNAData d = new CNAData(symbol + "-cna", symbol);
				d.data = new SingleCategoricalData[val.length];
				for (int i = 0; i < val.length; i++)
				{
					d.data[i] = new CNA(val[i]);
				}
				set.add(d);
			}
		}

		if (expReader != null)
		{
			double[] val = expReader.getGeneAlterationArray(symbol, samples);
			if (val != null)
			{
				RNAData d = new RNAData(symbol + "-rna", symbol);
				d.vals = val;
				set.add(d);
			}
		}

		if (mutReader != null && mutationEffectMap != null && mutationEffectMap.containsKey(symbol))
		{
			MutationData d = getMutations(symbol, mutationEffectMap.get(symbol));
			if (d != null) set.add(d);
		}

		dataCache.put(symbol, set);
		return set;
	}

	/**
	 * Associated the related data to the given set of relations.
	 */
	public void decorateRelations(Set<Relation> relations)
	{
		relations.stream().map(r -> new GeneWithData[]{r.sourceData, r.targetData}).flatMap(Arrays::stream).distinct()
			.forEach(d -> d.addAll(getData(d.getId())));
	}

	/**
	 * Fetches the mutation data from the reader and prepares for analysis.
	 * @param gene symbol
	 * @param mutEffectSign is it activating (1) or inhibiting (-1) mutation
	 */
	private MutationData getMutations(String gene, int mutEffectSign)
	{
		List<MutTuple>[] mutsList = mutReader.getMutations(gene, samples);
		if (mutsList == null) return null;

		MutationData mutData = new MutationData(gene + "-mut", gene, mutEffectSign);
		mutData.data = new SingleCategoricalData[samples.length];

		for (int i = 0; i < mutsList.length; i++)
		{
			int categ = 0;
			if (mutsList[i] == null)
			{
				categ = SingleCategoricalData.ABSENT;
			}
			else if (!mutsList[i].isEmpty())
			{
				categ = 1;
			}
			mutData.data[i] = new Mutation(categ, mutsList[i]);
		}
		return mutData;
	}

	/**
	 * Gets all the samples that exists in at least one type of dataset.
	 */
	public String[] getUnionSamples()
	{
		Set<String> set = new HashSet<>();
		if (cnaReader != null) set.addAll(cnaReader.getSamples());
		if (expReader != null) set.addAll(expReader.getSamples());
		if (mutReader != null) set.addAll(mutReader.getSamples());
		if (rppaReader != null) set.addAll(rppaReader.getSamples());

		List<String> list = new ArrayList<>(set);
		Collections.sort(list);
		return list.toArray(new String[list.size()]);
	}

	public Set<String> getGenesWithTotalProteinData()
	{
		return genesWithTotalProteinData;
	}

	/**
	 * Gets the set of genes with total protein RPPA measurement.
	 */
	private Set<String> findGenesWithTotalProteinData()
	{
		if (rppaReader == null) return Collections.emptySet();

		Set<String> genes = new HashSet<>();
		Set<String> sampleSet = rppaReader.getSamples();
		String[] samples = new String[]{sampleSet.iterator().next()};

		rppaReader.getGenes().stream().forEach(gene ->
			rppaReader.getAssociatedData(gene, samples).stream().filter(data -> !data.isPhospho())
				.forEach(data -> genes.addAll(data.genes)));

		return genes;
	}

	/**
	 * Puts the given change detector to the data that is filtered by the given selector.
	 */
	public void associateChangeDetector(OneDataChangeDetector chDet, DataSelector selector)
	{
		dataCache.values().stream().flatMap(Collection::stream).filter(selector::select).forEach(d -> d.setChDet(chDet));
	}

	/**
	 * Function to filter experiment data.
	 */
	public interface DataSelector
	{
		boolean select(ExperimentData data);
	}
}
