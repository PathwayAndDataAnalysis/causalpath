package org.panda.causalpath.resource;

import org.panda.causalpath.analyzer.OneDataChangeDetector;
import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;
import org.panda.resource.MatrixOfValuesDatasetReader;
import org.panda.resource.tcga.*;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Loads custom RNA data.
 */
public class RNALoader
{
	/**
	 * Reader for RPPA data.
	 */
	private MatrixOfValuesDatasetReader reader;

	/**
	 * When a data is already been read for a gene, this cache stores it so that creation of redundant data objects is
	 * avoided on multiple access.
	 */
	private Map<String, Set<ExperimentData>> dataCache;

	/**
	 * Names of samples loaded.
	 */
	private String[] samples;

	public RNALoader(String filename) throws IOException
	{
		reader = new MatrixOfValuesDatasetReader(filename);
		reader.load();
		Set<String> samples = reader.getSamples();
		this.samples = samples.toArray(new String[samples.size()]);
		dataCache = new HashMap<>();
	}

	public void setSamples(String[] samples)
	{
		this.samples = samples;
	}

	/**
	 * Gets the set of experiment data for the given symbol.
	 */
	public Set<ExperimentData> getData(String symbol)
	{
		// if asked before, return from cache
		if (dataCache.containsKey(symbol)) return dataCache.get(symbol);

		Set<ExperimentData> set = new HashSet<>();

		double[] val = reader.getGeneAlterationArray(symbol, samples);
		if (val != null)
		{
			RNAData d = new RNAData(symbol + "-rna", symbol);
			d.vals = val;
			set.add(d);
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
