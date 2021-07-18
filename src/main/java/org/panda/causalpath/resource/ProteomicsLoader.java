package org.panda.causalpath.resource;

import org.panda.causalpath.analyzer.CausalityHelper;
import org.panda.causalpath.analyzer.OneDataChangeDetector;
import org.panda.causalpath.data.*;
import org.panda.causalpath.network.Relation;
import org.panda.resource.tcga.ProteomicsFileRow;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;

import java.util.*;

/**
 * Coverts the RPPAData in resource project to appropriate objects for this project and serves them.
 */
public class ProteomicsLoader
{
	/**
	 * Map from genes to related data.
	 */
	Map<String, Set<ExperimentData>> dataMap;

	/**
	 * Set of all data. This collection supposed to hold everything in the dataMap's values.
	 */
	Set<ExperimentData> datas;

	/**
	 * Initializes self using a set of RPPAData.
	 */
	public ProteomicsLoader(Collection<ProteomicsFileRow> rows, Map<DataType, Double> stdevThresholds)
	{
		dataMap = new HashMap<>();
		datas = new HashSet<>();
		rows.stream().distinct().forEach(r ->
		{
			ExperimentData ed = r.isActivity() ? new ActivityData(r) :
				r.isSiteSpecific() ? new SiteModProteinData(r) :
				r.isRNA() ? new RNAData(r) :
				r.isMetabolite() ? new MetaboliteData(r) :
					new ProteinData(r);

			if (stdevThresholds != null && ed instanceof NumericData)
			{
				DataType type = ed.getType();
				Double thr = stdevThresholds.get(type);
				if (thr != null)
				{
					double sd = Summary.stdev(((NumericData) ed).vals);
					if (Double.isNaN(sd) || sd < thr) return;
				}
			}

			for (String sym : ed.getGeneSymbols())
			{
				if (!dataMap.containsKey(sym)) dataMap.put(sym, new HashSet<>());

				// check if there is already some data with the same ID
				for (ExperimentData data : dataMap.get(sym))
				{
					if (data.getId().equals(ed.getId()))
					{
						throw new RuntimeException("Proteomic data has non-unique IDs: " + ed.getId());
					}
				}

				dataMap.get(sym).add(ed);
				datas.add(ed);
			}
		});
	}

	public void addRepeatData(Collection<ProteomicsFileRow> rows, Map<DataType, Double> stdevThresholds)
	{
		rows.stream().distinct().forEach(r ->
		{
			ExperimentData ed = r.isActivity() ? new ActivityData(r) :
				r.isSiteSpecific() ? new SiteModProteinData(r) : new ProteinData(r);

			if (stdevThresholds != null && ed instanceof NumericData)
			{
				DataType type = ed.getType();
				Double thr = stdevThresholds.get(type);
				if (thr != null)
				{
					double sd = Summary.stdev(((NumericData) ed).vals);
					if (Double.isNaN(sd) || sd < thr) return;
				}
			}

			for (String sym : ed.getGeneSymbols())
			{
				if (!dataMap.containsKey(sym)) dataMap.put(sym, new HashSet<>());

				ExperimentData orig = null;
				for (ExperimentData data : dataMap.get(sym))
				{
					if (data.getId().equals(ed.getId()))
					{
						orig = data;
						break;
					}
				}

				if (orig != null)
				{
					orig.addRepeatData(ed);
				}
				else
				{
					dataMap.get(sym).add(ed);
					datas.add(ed);
				}
			}
		});
	}

	public void printStDevHistograms()
	{
		printStDevHistograms(datas);
	}
	public void printStDevHistograms(Set<ExperimentData> datas)
	{
		System.out.println("\nSt-dev histograms:");
		Map<DataType, Histogram> hMap = new HashMap<>();
		datas.stream().filter(d -> d instanceof NumericData).map(d -> (NumericData) d).forEach(d ->
		{
			DataType type = d.getType();
			if (!hMap.containsKey(type))
			{
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
	public void decorateRelations(Set<Relation> relations)
	{
		Map<String, GeneWithData> map = collectExistingData(relations);
		CausalityHelper ch = new CausalityHelper();
		for (Relation rel : relations)
		{
			if (rel.sourceData == null)
			{
				if (map.containsKey(rel.source))
				{
					rel.sourceData = map.get(rel.source);
				}
				else
				{
					GeneWithData gwd = new GeneWithData(rel.source);
					gwd.addAll(dataMap.get(rel.source));
					map.put(gwd.getId(), gwd);
				}
			}
			else rel.sourceData.addAll(dataMap.get(rel.source));

			if (rel.targetData == null)
			{
				if (map.containsKey(rel.target))
				{
					rel.targetData = map.get(rel.target);
				}
				else
				{
					GeneWithData gwd = new GeneWithData(rel.target);
					gwd.addAll(dataMap.get(rel.target));
					map.put(gwd.getId(), gwd);
				}
			}
			else rel.targetData.addAll(dataMap.get(rel.target));

			rel.chDet = ch;
		}
	}

	private Map<String, GeneWithData> collectExistingData(Set<Relation> relations)
	{
		Map<String, GeneWithData> map = new HashMap<>();

		for (Relation rel : relations)
		{
			if (rel.sourceData != null) map.put(rel.sourceData.getId(), rel.sourceData);
			if (rel.targetData != null) map.put(rel.targetData.getId(), rel.targetData);
		}

		return map;
	}

	/**
	 * Puts the given change detector to the data that is filtered by the given selector.
	 */
	public void associateChangeDetector(OneDataChangeDetector chDet, DataSelector selector)
	{
		datas.stream().filter(selector::select).forEach(d -> d.setChDet(chDet));
	}

	public void initMissingDataForProteins()
	{
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
			if (d.getType() == DataType.PROTEIN)
			{
				for (int i = 0; i < size; i++)
				{
					if (!Double.isNaN(d.vals[i]))
					{
						hasTotalProt[i] = true;
//						totalProtCnt[i]++;
					}
				}
			}
			else if (d.getType() == DataType.PHOSPHOPROTEIN)
			{
				for (int i = 0; i < size; i++)
				{
					if (!Double.isNaN(d.vals[i]))
					{
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
	public interface DataSelector
	{
		boolean select(ExperimentData data);
	}
}
