package org.panda.causalpath.analyzer;

import org.panda.causalpath.data.*;
import org.panda.causalpath.network.GraphFilter;
import org.panda.causalpath.network.Relation;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Tuple;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This class matches the experiment data with the pathway relations, and detect potential causality.
 */
public class CausalitySearcher implements Cloneable
{
	/**
	 * If that is false, then we are interested in conflicting relations.
	 */
	private int causal;

	/**
	 * If this is false, then we don't care if the target site of a relation and the site on the data matches.
	 */
	private boolean forceSiteMatching;

	/**
	 * When site matching is on, this parameter determines the largest value of the mismatch between relation and data
	 * to consider it as a match. Set this to 0 for the most strict match.
	 */
	private int siteProximityThreshold;

	/**
	 * When this is true, the data whose effect is 0, but potentially can be part of the causal network if known, are
	 * collected.
	 */
	boolean collectDataWithMissingEffect;

	/**
	 * The interesting subset of phosphorylation data with unknown effect.
	 */
	private Set<SiteModProteinData> dataNeedsAnnotation;

	/**
	 * When this is true, the data that are used for inference of causality or conflict are saved in a set.
	 */
	boolean collectDataUsedForInference;

	/**
	 * The set of data that is used during inference of causality.
	 */
	private Map<Relation, Set<ExperimentData>> dataUsedForInference;

	/**
	 * Set of pairs of data that contributed to inference.
	 */
	private Map<Relation, Set<List<ExperimentData>>> pairsUsedForInference;

	/**
	 * If true, then an expression relation has to have an Activity data at its source.
	 */
	protected boolean mandateActivityDataUpstreamOfExpression;

	/**
	 * The data types that can be used for evidence of expression change. This is typically protein or rna or both.
	 */
	protected DataType[] expressionEvidence;

	/**
	 * When true, this class uses only the strongest changing proteomic data with known effect as general activity
	 * evidence.
	 */
	protected boolean useStrongestProteomicsDataForActivity;

	/**
	 * When true, if a node has activity data, other data types are ignored for providing activity evidence.
	 */
	protected boolean prioritizeActivityData;

	/**
	 * Data types that indicate activity change.
	 */
	Set<DataType> generalActivityChangeIndicators;

	/**
	 * A graph filter if needed to use at the end of network inference.
	 */
	GraphFilter graphFilter;

	Map<Relation, Set<ExperimentData>> affectingSourceDataMap;
	Map<Relation, Set<ExperimentData>> explainableTargetDataMap;


	/**
	 * Constructor with the reasoning type.
	 * @param causal true:causal, false:conflicting
	 */
	public CausalitySearcher(boolean causal)
	{
		setCausal(causal);
		this.forceSiteMatching = true;
		this.siteProximityThreshold = 0;
		this.collectDataWithMissingEffect = true;
		this.collectDataUsedForInference = true;
		this.mandateActivityDataUpstreamOfExpression = false;
		this.useStrongestProteomicsDataForActivity = false;

		this.generalActivityChangeIndicators = new HashSet<>(Arrays.asList(DataType.PROTEIN, DataType.PHOSPHOPROTEIN,
			DataType.ACETYLPROTEIN, DataType.METHYLPROTEIN, DataType.METABOLITE, DataType.ACTIVITY));
	}

	public void initRelationDataMappingMemory()
	{
		this.affectingSourceDataMap = new HashMap<>();
		this.explainableTargetDataMap = new HashMap<>();
	}

	public CausalitySearcher copy()
	{
		try
		{
			CausalitySearcher cs = (CausalitySearcher) this.clone();
			cs.pairsUsedForInference = null;
			cs.dataUsedForInference = null;
			cs.dataNeedsAnnotation = null;
			cs.generalActivityChangeIndicators = new HashSet<>(generalActivityChangeIndicators);
			cs.setCollectDataUsedForInference(false);
			cs.setCollectDataWithMissingEffect(false);
			return cs;
		}
		catch (CloneNotSupportedException e)
		{
			throw new RuntimeException(e);
		}
	}

	/**
	 * Finds compatible or conflicting relations. The relations have to be associated with experiment data. Both the
	 * experiment data and the relations have to be associated with related change detectors.
	 */
	public Set<Relation> run(Set<Relation> relations)
	{
//		printSizeOfRelationsBetweenSignificantData(relations);

		// Initialize collections

		if (collectDataWithMissingEffect)
		{
			if (dataNeedsAnnotation == null) dataNeedsAnnotation = new HashSet<>();
			else dataNeedsAnnotation.clear();
		}
		if (collectDataUsedForInference)
		{
			if (dataUsedForInference == null)
			{
				dataUsedForInference = new HashMap<>();
				pairsUsedForInference = new HashMap<>();
			}
			else
			{
				dataUsedForInference.clear();
				pairsUsedForInference.clear();
			}
		}

		// This is where magic happens
		Set<Relation> results = relations.stream().filter(this::satisfiesCriteria).collect(Collectors.toSet());

		// If a subset of the results is desired, trim it
		if (graphFilter != null)
		{
			results = graphFilter.postAnalysisFilter(results);

			// remove unnecessary entries in the collected data
			if (collectDataUsedForInference)
			{
				Set<Relation> removedRels = new HashSet<>(dataUsedForInference.keySet());
				removedRels.removeAll(results);
				removedRels.forEach(r ->
				{
					dataUsedForInference.remove(r);
					pairsUsedForInference.remove(r);
				});
			}
		}

		return results;
	}

	/**
	 * Checks if the relation explains/conflicts the associated data.
	 * @param relation relation to check
	 * @return true if explanatory
	 */
	public boolean satisfiesCriteria(Relation relation)
	{
		// Get data of target gene this relation can explain the change
		if (explainableTargetDataMap != null && !explainableTargetDataMap.containsKey(relation))
			explainableTargetDataMap.put(relation, getExplainableTargetDataWithSiteMatch(relation));

		Set<ExperimentData> td = explainableTargetDataMap == null ?
			getExplainableTargetDataWithSiteMatch(relation) : explainableTargetDataMap.get(relation);

		if (!td.isEmpty())
		{
			// Get data of the source gene that can be cause of this relation
			if (affectingSourceDataMap != null && !affectingSourceDataMap.containsKey(relation))
				affectingSourceDataMap.put(relation, getAffectingSourceData(relation));

			Set<ExperimentData> sd = affectingSourceDataMap == null ?
				getAffectingSourceData(relation) : affectingSourceDataMap.get(relation);

			if (!sd.isEmpty())
			{
				return satisfiesCriteria(sd, relation, td);
			}
		}
		return false;
	}

	/**
	 * Checks if the relation has potential to explain/conflict with the associated data to source and targets, but
	 * without evaluating data values.
	 *
	 * @param relation relation to check
	 * @return true if the relation has considerable data at both sides
	 */
	public boolean hasConsiderableDownstreamData(Relation relation)
	{
		return !getExplainableTargetDataWithSiteMatch(relation).isEmpty();
	}

	/**
	 * Checks if the relation has potential to explain/conflict with the associated data to source and targets, but
	 * without evaluating data values.
	 *
	 * @param relation relation to check
	 * @return true if the relation has considerable data at both sides
	 */
	public boolean hasConsiderableData(Relation relation)
	{
		if (hasConsiderableDownstreamData(relation))
		{
			// Get data of the source gene that can be cause of this relation
			Set<ExperimentData> sd = getAffectingSourceData(relation);
			return !sd.isEmpty();
		}
		return false;
	}

	/**
	 * Checks if any pair from the given source and target data can be explained by the given relation.
	 * @param sd source data
	 * @param rel the relation
	 * @param td target data
	 * @return true if any sd td pair is explained/conflicted by the given relation
	 */
	private boolean satisfiesCriteria(Set<ExperimentData> sd, Relation rel, Set<ExperimentData> td)
	{
		boolean satisfies = false;

		for (ExperimentData sourceData : sd)
		{
			for (ExperimentData targetData : td)
			{
				if (satisfiesCriteria(rel, sourceData, targetData))
				{
					satisfies = true;
				}
			}
		}

		return satisfies;
	}

	/**
	 * Checks if the relation can explain/conflict the given source target data pair.
	 * @param rel the relation
	 * @param sourceData the source data
	 * @param targetData the target data
	 * @return true if the relation can explain/conflict the given data pair
	 */
	private boolean satisfiesCriteria(Relation rel, ExperimentData sourceData, ExperimentData targetData)
	{
		int e = rel.chDet.getChangeSign(sourceData, targetData) * rel.getSign();

		if (e != 0 && collectDataWithMissingEffect && sourceData.getEffect() == 0)
		{
			dataNeedsAnnotation.add((SiteModProteinData) sourceData);
		}
		else if (sourceData.getEffect() * e == causal)
		{
			if (collectDataUsedForInference)
			{
				if (!dataUsedForInference.containsKey(rel))
				{
					dataUsedForInference.put(rel, new HashSet<>());
					pairsUsedForInference.put(rel, new HashSet<>());
				}

				dataUsedForInference.get(rel).add(sourceData);
				dataUsedForInference.get(rel).add(targetData);
				pairsUsedForInference.get(rel).add(Arrays.asList(sourceData, targetData));
			}
			return true;
		}
		return false;
	}

	/**
	 * Gets the source data that match with the given relation and target data.
	 * @param rel the relation
	 * @param target target data
	 * @return the set of matching source data
	 */
	public Set<ExperimentData> getSatisfyingSourceData(Relation rel, ExperimentData target)
	{
		return getAffectingSourceData(rel).stream()
			.filter(source -> satisfiesCriteria(rel, source, target))
			.collect(Collectors.toSet());
	}

	/**
	 * Gets the set of target data that the given relation can potentially explain its change. but that is without
	 * checking the value changes, only by data types.
	 * @param rel the relation
	 * @return the set of target data explainable by the relation
	 */
	public Set<ExperimentData> getExplainableTargetData(Relation rel)
	{
		if (rel.type.affectsPhosphoSite)
		{
			return getEvidenceForPhosphoChange(rel.targetData);
		}
		else if (rel.type.affectsTotalProt)
		{
			return getEvidenceForExpressionChange(rel.targetData);
		}
		else if (rel.type.affectsGTPase)
		{
			return getEvidenceForGTPaseChange(rel.targetData);
		}
		else if (rel.type.affectsAcetylSite)
		{
			return getEvidenceForAcetylChange(rel.targetData);
		}
		else if (rel.type.affectsMethlSite)
		{
			return getEvidenceForMethylChange(rel.targetData);
		}
		else if (rel.type.affectsMetabolite)
		{
			return getEvidenceForMetaboliteChange(rel.targetData);
		}

		throw new RuntimeException("Code should not reach here. Is there a new relation type to handle?");
	}

	/**
	 * Gets the set of target data that the given relation can potentially explain its change. That is without
	 * checking the value changes, but site matching (if relevant) is evaluated.
	 * @param rel the relation
	 * @return the set of target data explainable by the relation
	 */
	public Set<ExperimentData> getExplainableTargetDataWithSiteMatch(Relation rel)
	{
		Set<ExperimentData> datas = getExplainableTargetData(rel);

		return datas.stream().filter(d -> !(d instanceof SiteModProteinData) ||
			isTargetSiteCompatible(rel, (SiteModProteinData) d)).collect(Collectors.toSet());
	}


	/**
	 * Checks if target sites match for the relation and the data.
	 * @param rel the relation
	 * @param target target data
	 * @return true if there is a match or no match needed
	 */
	public boolean isTargetSiteCompatible(Relation rel, SiteModProteinData target)
	{
		return !forceSiteMatching || !CollectionUtil.intersectionEmpty(
				rel.getTargetWithSites(siteProximityThreshold), target.getGenesWithSites());
	}

	/**
	 * Gets the source data that can be cause of this relation. This is without checking any change in values, only by
	 * data types.
	 * @param rel the relation
	 * @return the set of source data that can be affecting
	 */
	public Set<ExperimentData> getAffectingSourceData(Relation rel)
	{
		if (rel.type.affectsPhosphoSite || rel.type.affectsAcetylSite || rel.type.affectsMethlSite)
		{
			return getUpstreamEvidenceForSiteSpecificChange(rel.sourceData);
		}
		else if (rel.type.affectsTotalProt)
		{
			return getUpstreamEvidenceForExpressionChange(rel.sourceData);
		}
		else if (rel.type.affectsGTPase)
		{
			return getUpstreamEvidenceForGTPaseChange(rel.sourceData);
		}
		else if (rel.type.affectsMetabolite)
		{
			return getUpstreamEvidenceForMetaboliteChange(rel.sourceData);
		}

		throw new RuntimeException("Code should not reach here. There must be a new relation type to handle.");
	}

	/**
	 * Gets the evidence for activity change of a gene in terms of the associated data. This method does not evaluate a
	 * change in values but only assesses the valid data types.
	 * @param gene the gene
	 * @return the data with potential to indicate activity change
	 */
	public Set<ExperimentData> getGeneralActivationEvidence(GeneWithData gene)
	{
		Set<ExperimentData> set = new HashSet<>();

		for (DataType type : generalActivityChangeIndicators)
		{
			set.addAll(gene.getData(type));
		}

		set = set.stream().filter(d -> d.getEffect() != 0).collect(Collectors.toSet());

		if (useStrongestProteomicsDataForActivity)
		{
			removeShadowedProteomicData(set);
		}

		if (prioritizeActivityData)
		{
			removeOtherDataIfActivityDataIsPresent(set);
		}

		return set;
	}

	public Set<ExperimentData> getUpstreamEvidenceForSiteSpecificChange(GeneWithData gene)
	{
		return getGeneralActivationEvidence(gene);
	}

	public Set<ExperimentData> getUpstreamEvidenceForExpressionChange(GeneWithData gene)
	{
		if (mandateActivityDataUpstreamOfExpression) return gene.getData(DataType.ACTIVITY);
		else return getGeneralActivationEvidence(gene);
	}

	public Set<ExperimentData> getUpstreamEvidenceForGTPaseChange(GeneWithData gene)
	{
		return getGeneralActivationEvidence(gene);
	}

	public Set<ExperimentData> getUpstreamEvidenceForMetaboliteChange(GeneWithData gene)
	{
		return getGeneralActivationEvidence(gene);
	}

	public Set<ExperimentData> getEvidenceForPhosphoChange(GeneWithData gene)
	{
		return gene.getData(DataType.PHOSPHOPROTEIN);
	}

	public Set<ExperimentData> getEvidenceForExpressionChange(GeneWithData gene)
	{
		if (expressionEvidence != null) return gene.getData(expressionEvidence);
		return gene.getData(DataType.PROTEIN);
	}

	public Set<ExperimentData> getEvidenceForGTPaseChange(GeneWithData gene)
	{
		return gene.getData(DataType.ACTIVITY);
	}

	public Set<ExperimentData> getEvidenceForAcetylChange(GeneWithData gene)
	{
		return gene.getData(DataType.ACETYLPROTEIN);
	}

	public Set<ExperimentData> getEvidenceForMethylChange(GeneWithData gene)
	{
		return gene.getData(DataType.METHYLPROTEIN);
	}

	public Set<ExperimentData> getEvidenceForMetaboliteChange(GeneWithData gene)
	{
		return gene.getData(DataType.METABOLITE);
	}


	/**
	 * This method iterates over total protein and phosphoprotein data that has a known effect, and leaves only the one
	 * with the biggest change, removes others. This is sometimes useful for complexity management.
	 * @param data data to select from
	 */
	protected void removeShadowedProteomicData(Set<ExperimentData> data)
	{
		if (data.size() > 1)
		{
			Optional<ProteinData> opt = data.stream().filter(d -> d instanceof ProteinData && d.getEffect() != 0)
				.map(d -> (ProteinData) d).sorted((d1, d2) ->
					Double.compare(Math.abs(d2.getChangeValue()), Math.abs(d1.getChangeValue()))).findFirst();

			if (opt.isPresent())
			{
				ExperimentData ed = opt.get();
				data.removeIf(d -> d instanceof ProteinData && d != ed);
			}
		}
	}

	protected void removeOtherDataIfActivityDataIsPresent(Set<ExperimentData> data)
	{
		if (data.stream().anyMatch(d -> d instanceof ActivityData))
		{
			data.retainAll(data.stream().filter(d -> d instanceof ActivityData).collect(Collectors.toSet()));
		}
	}

	public void writeResults(String filename) throws IOException
	{
		if (pairsUsedForInference.isEmpty()) return;

		TwoDataChangeDetector relDet = pairsUsedForInference.keySet().iterator().next().chDet;
		CorrelationDetector corDet = relDet instanceof CorrelationDetector ? (CorrelationDetector) relDet : null;

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		writer.write("Source\tRelation\tTarget\tSites\t");
		if (corDet != null) writer.write("Source data ID\tTarget data ID\tCorrelation\tCorrelation pval");
		else writer.write("Source data ID\tSource change\t Source change pval\tTarget data ID\tTarget change\tTarget change pval");

		pairsUsedForInference.keySet().stream().
			sorted(Comparator.comparing(this::getRelationScore).reversed()). // Sort relations to their significance
			forEach(r -> pairsUsedForInference.get(r).stream().forEach(pair ->
		{
			Iterator<ExperimentData> iter = pair.iterator();
			ExperimentData sourceData = iter.next();
			ExperimentData targetData = iter.next();

			FileUtil.lnwrite(r.source + "\t" + r.type.getName() + "\t" + r.target + "\t" + r.getSitesInString() + "\t", writer);

			if (corDet != null)
			{
				Tuple t = corDet.calcCorrelation(sourceData, targetData);
				FileUtil.write(sourceData.getId() + "\t" + targetData.getId() + "\t" + t.v + "\t" + t.p, writer);
			}
			else
			{
				OneDataChangeDetector sDet = sourceData.getChDet();
				OneDataChangeDetector tDet = targetData.getChDet();

				FileUtil.write(sourceData.getId() + "\t" + sourceData.getChangeValue() + "\t" +
					(sDet instanceof SignificanceDetector ? ((SignificanceDetector) sDet).getPValue(sourceData) : "") + "\t", writer);

				FileUtil.write(targetData.getId() + "\t" + targetData.getChangeValue() + "\t" +
					(tDet instanceof SignificanceDetector ? ((SignificanceDetector) tDet).getPValue(targetData) : ""), writer);
			}
		}));
		writer.close();
	}

	/**
	 * We want to sort the result rows according the their significance. This method generates a score for each relation
	 * so that we can sort them using that score.
	 */
	private double getRelationScore(Relation r)
	{
		double max = -Double.MAX_VALUE;
		Set<List<ExperimentData>> pairs = pairsUsedForInference.get(r);
		for (List<ExperimentData> pair : pairs)
		{
			ExperimentData src = pair.get(0);
			ExperimentData tgt = pair.get(1);

			double val;

			if (r.chDet instanceof CorrelationDetector)
			{
				val = -((CorrelationDetector) r.chDet).calcCorrelation(src, tgt).p;
			}
			else
			{
				OneDataChangeDetector sDet = src.getChDet();
				OneDataChangeDetector tDet = tgt.getChDet();

				if (sDet instanceof SignificanceDetector || tDet instanceof SignificanceDetector)
				{
					double vS = sDet instanceof SignificanceDetector ? ((SignificanceDetector) sDet).getPValue(src) : 0;
					double vT = tDet instanceof SignificanceDetector ? ((SignificanceDetector) tDet).getPValue(tgt) : 0;
					val = -Math.max(vS, vT);
				}
				else
				{
					val = Math.min(Math.abs(sDet.getChangeValue(src)), Math.abs(tDet.getChangeValue(tgt)));
				}
			}

			if (val > max) max = val;
		}
		return max;
	}

	public void setCausal(boolean causal)
	{
		this.causal = causal ? 1 : -1;
	}

	public Set<SiteModProteinData> getDataNeedsAnnotation()
	{
		return dataNeedsAnnotation;
	}

	public Set<ExperimentData> getDataUsedForInference()
	{
		return dataUsedForInference.values().stream().flatMap(Collection::stream).collect(Collectors.toSet());
	}

	public Set<List<ExperimentData>> getPairsUsedForInference()
	{
		return pairsUsedForInference.values().stream().flatMap(Collection::stream).collect(Collectors.toSet());
	}

	//--DEBUG----------
	public void writePairsUsedForInferenceWithCorrelations(String file)
	{
		try
		{
			BufferedWriter writer = Files.newBufferedWriter(Paths.get(file));
			Set<List<ExperimentData>> pairs = getPairsUsedForInference();

			CorrelationDetector cd = new CorrelationDetector(-1, 1);
			Map<String, Double> pvals = new HashMap<>();

			for (List<ExperimentData> pair : pairs)
			{
				Iterator<ExperimentData> iter = pair.iterator();
				ExperimentData data1 = iter.next();
				ExperimentData data2 = iter.next();

				Tuple corr = cd.calcCorrelation(data1, data2);
				if (!corr.isNaN())
				{
					StringBuilder sb = new StringBuilder();
					pair.stream().sorted(Comparator.comparing(ExperimentData::getId))
						.forEach(e -> sb.append(e.getId()).append(":"));
					String id = sb.toString();

					pvals.put(id, corr.p);
				}
			}

			pvals.forEach((s, p) -> FileUtil.writeln(s + "\t" + p, writer));

			writer.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	public void printSizeOfRelationsBetweenSignificantData(Set<Relation> relations)
	{
		Set<Relation> rels = relations.stream().filter(r -> r.sourceData.getDataStream().anyMatch(d -> d.getChangeSign() != 0) &&
			r.targetData.getDataStream().anyMatch(d -> d.getChangeSign() != 0)).collect(Collectors.toSet());
		System.out.println("rels between sig data = " + rels.size());
		long protCnt = rels.stream().map(r -> Arrays.asList(r.source, r.target)).flatMap(Collection::stream).distinct().count();
		System.out.println("protCnt = " + protCnt);
	}
	//--DEBUG----------

	public Map<Relation, Set<ExperimentData>> getInferenceUnits()
	{
		return dataUsedForInference;
	}

	public void setCollectDataUsedForInference(boolean collect)
	{
		this.collectDataUsedForInference = collect;
	}

	public void setCollectDataWithMissingEffect(boolean collect)
	{
		this.collectDataWithMissingEffect = collect;
	}

	public void setMandateActivityDataUpstreamOfExpression(boolean mandate)
	{
		this.mandateActivityDataUpstreamOfExpression = mandate;
	}

	public void setUseStrongestProteomicsDataForActivity(boolean use)
	{
		this.useStrongestProteomicsDataForActivity = use;
	}


	public boolean getPrioritizeActivityData() { return this.prioritizeActivityData; }

	public void setPrioritizeActivityData(boolean prioritizeActivityData)
	{
		this.prioritizeActivityData = prioritizeActivityData;
	}

	public void addDataTypeForGeneralActivity(DataType type)
	{
		generalActivityChangeIndicators.add(type);
	}

	public boolean getForceSiteMatching() { return this.forceSiteMatching; }

	public void setForceSiteMatching(boolean forceSiteMatching)
	{
		this.forceSiteMatching = forceSiteMatching;
	}

	public int getSiteProxitmityThreshold() { return this.siteProximityThreshold; }

	public void setSiteProximityThreshold(int siteProximityThreshold)
	{
		this.siteProximityThreshold = siteProximityThreshold;
	}

	public void setGraphFilter(GraphFilter graphFilter)
	{
		this.graphFilter = graphFilter;
	}

	public boolean hasNoGraphFilter()
	{
		return graphFilter == null;
	}

	public GraphFilter getGraphFilter()
	{
		return graphFilter;
	}

	public void setExpressionEvidence(DataType... types)
	{
		expressionEvidence = types;
	}

	/**
	 * This method is experimental, only to test how good is RNA expression as a proxy to protein activity. It is not
	 * meant to be used in a regular CausalPath analysis.
	 */
	public void useExpressionForActivity()
	{
		this.generalActivityChangeIndicators = new HashSet<>(Collections.singletonList(DataType.RNA));
	}
}
