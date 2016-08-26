package org.panda.causalpath.resource;

import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;
import org.panda.resource.tcga.*;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * Created by babur on 4/5/16.
 */
public class TCGALoader
{
	private CNAReader cnaReader;
	private ExpressionReader expReader;
	private MutationReader mutReader;
	private RPPAReader rppaReader;
	// todo add methylation data

	private Map<String, ProteinData> rppaCache;
	private Map<String, Set<ExperimentData>> dataCache;
	private Set<String> genesWithTotalProteinData;

	private String[] samples;

	private static final String COPY_NUMBER_FILE = "/copynumber.txt";
	private static final String MUTATION_FILE = "/mutation.maf";
	private static final String EXPRESSION_FILE = "/expression.txt";
	private static final String RPPA_FILE = "/rppa.txt";

	public TCGALoader(String dir)
	{
//		try{cnaReader = new CNAReader(dir + COPY_NUMBER_FILE, false, 0);} catch (FileNotFoundException e){}
//		try{expReader = new ExpressionReader(dir + EXPRESSION_FILE);} catch (FileNotFoundException e){}
		try{mutReader = new MutationReader(dir + MUTATION_FILE);} catch (IOException e){}
		try{rppaReader = new RPPAReader(dir + RPPA_FILE);} catch (FileNotFoundException e){}
		this.samples = getUnionSamples();
		if (rppaReader != null) rppaCache = new HashMap<>();
		dataCache = new HashMap<>();
		genesWithTotalProteinData = findGenesWithTotalProteinData();
	}

	public Set<ExperimentData> getData(String symbol)
	{
		if (dataCache.containsKey(symbol)) return dataCache.get(symbol);

		Set<ExperimentData> set = new HashSet<>();
		if (rppaReader != null)
		{
			Set<RPPAData> rppaSet = rppaReader.getAssociatedData(symbol, samples);

			rppaSet.stream().filter(rppa -> !rppaCache.containsKey(rppa.id)).forEach(rppa -> {
				ProteinData d = rppa.isPhospho() ? new PhosphoProteinData(rppa) : new ProteinData(rppa);
				rppaCache.put(d.id, d);
			});

			rppaSet.stream().map(rppa -> rppaCache.get(rppa.id)).forEach(set::add);
		}
		if (cnaReader != null)
		{
			int[] val = cnaReader.getGeneAlterationArray(symbol, samples);
			if (val != null)
			{
				CNAData d = new CNAData(symbol + "-cna", symbol);
				d.data = new SingleQData[val.length];
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
				ExpressionData d = new ExpressionData(symbol + "-rna", symbol);
				d.vals = val;
				set.add(d);
			}
		}
		if (mutReader != null)
		{
			MutationData d = getMutations(symbol, MutationClassifier.isActivatedByMutation(symbol));
			if (d != null) set.add(d);
		}
		dataCache.put(symbol, set);
		return set;
	}

	public void decorateRelations(Set<Relation> relations)
	{
		for (Relation relation : relations)
		{
			relation.sourceData.addAll(getData(relation.source));
			relation.targetData.addAll(getData(relation.target));
		}
	}

	private MutationData getMutations(String gene, boolean activated)
	{
		List<MutTuple>[] mutsList = mutReader.getMutations(gene, samples);
		if (mutsList == null) return null;

		MutationData mutData = new MutationData(gene + "-mut", gene);
		mutData.data = new SingleQData[samples.length];

		for (int i = 0; i < mutsList.length; i++)
		{
			int categ = 0;
			if (mutsList[i] == null)
			{
				categ = SingleQData.ABSENT;
			}
			else if (!mutsList[i].isEmpty())
			{
				categ = activated ? 1 : -1;
			}
			mutData.data[i] = new Mutation(categ, mutsList[i]);
		}
		return mutData;
	}

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

}
