package org.panda.causalpath.run;

import com.github.jsonldjava.utils.JsonUtils;
import org.panda.causalpath.analyzer.*;
import org.panda.causalpath.data.*;
import org.panda.causalpath.network.GraphFilter;
import org.panda.causalpath.network.GraphWriter;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationType;
import org.panda.causalpath.resource.NetworkLoader;
import org.panda.causalpath.resource.ProteomicsFileReader;
import org.panda.causalpath.resource.ProteomicsLoader;
import org.panda.causalpath.resource.TCGALoader;
import org.panda.resource.ResourceDirectory;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.resource.tcga.ProteomicsFileRow;
import org.panda.utility.ArrayUtil;
import org.panda.utility.BooleanMatrixRandomizer;
import org.panda.utility.FileUtil;
import org.panda.utility.RandomizedMatrices;
import org.panda.utility.statistics.FDR;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This method is built to handle running the causality analysis using its jar file and pointing it to a related
 * resource directory.
 *
 * @author Ozgun Babur
 */
public class CausalPath
{
	/**
	 * Name of the parameters file.
	 */
	public static final String PARAMETER_FILENAME = "parameters.txt";
	public static final String CAUSATIVE_RESULT_FILE_PREFIX = "causative";
	public static final String CAUSATIVE_RESULT_FILE_DATA_CENTRIC_PREFIX = "causative-data-centric";
	public static final String CONFLICTING_RESULT_FILE_PREFIX = "conflicting";
	public static final String UNKNOWN_SITE_EFFECT_FILENAME = "unknown-site-effects.txt";
	public static final String SIGNIFICANCE_FILENAME = "significance-pvals.txt";
	public static final String VALUE_CHANGES_FILENAME = "value-changes.txt";

	/**
	 * The directory where the parameters.txt file resides in.
	 */
	private String directory;

	/**
	 * Name of the proteomics platform file. Can be the same with the values file.
	 */
	private String proteomicsPlatformFile;

	/**
	 * Name of the proteomics values file. Can be the same with the platform file.
	 */
	private String proteomicsValuesFile;

	/**
	 * Name of the column that contains IDs. This name have to be the same in both the platform and the values files.
	 */
	private String IDColumn;

	/**
	 * Name of the column that contain gene symbols related to the rows in proteomics file.
	 */
	private String symbolsColumn;

	/**
	 * Name of the column that contain protein modification sites related to the rows in proteomics file.
	 */
	private String sitesColumn;

	/**
	 * The name of the column that contains modified site effects.
	 */
	private String effectColumn;

	/**
	 * Names of the value columns when there is only one group of values.
	 */
	private List<String> valueColumn;

	/**
	 * Names of the columns that contain control values.
	 */
	private List<String> controlValueColumn;

	/**
	 * Names of the columns that contain test values.
	 */
	private List<String> testValueColumn;

	/**
	 * Whether or not to do a log transformation on the values read, before every other kinds of processing. Log base
	 * is 2.
	 */
	private boolean doLogTransfrorm = false;

	/**
	 * A threshold to determine significant changes whenever it applies. This can be a value threshold, or a p-value
	 * threshold depending on the value transformation type.
	 */
	private Map<DataType, Double> thresholdForDataSignificance;

	/**
	 * An FDR threshold for the data that provides us p-values.
	 */
	private Map<DataType, Double> fdrThresholdForDataSignificance;

	/**
	 * Set of standard deviation thresholds to apply to the data.
	 */
	private Map<DataType, Double> stDevThresholds;

	/**
	 * Whether we should pool proteomics and phosphoproteomics during FDR adjustment. Pooling makes sense for RPPA data
	 * because they are coming from single experiment.
	 */
	private boolean poolProteomicsForFDRAdjustment = false;

	/**
	 * A threshold to determine if the network size or a downstream of a protein is significantly crowded.
	 */
	private double fdrThresholdForNetworkSignificance = 0.1;

	/**
	 * A correlation threshold when the value transformation is correlation. This is optional, but if not used, then
	 * using the <code>thresholdForSignificance</code> is mandatory.
	 */
	private double pvalThresholdForCorrelation = -1;

	private double correlationValueThreshold = -1;

	/**
	 * A correlation FDR threshold when the value transformation is correlation. When this is used,
	 */
	private double fdrThresholdForCorrelation = -1;

	/**
	 * Proximity threshold to infer effects of unknown sites by the known effect of neighbor sites.
	 */
	private int siteEffectProximityThreshold = 0;

	/**
	 * A method for interpretation of values in the data file.
	 */
	private ValueTransformation transformation;

	/**
	 * A default value for the missing values in the data file.
	 */
	private Double defaultMissingValue = null;

	/**
	 * Parameter for calculating network significances. If this is true, then <code>permutationCount</code> should also
	 * be set.
	 */
	private boolean calculateNetworkSignificance = false;

	/**
	 * Number of iterations of permutations to use during network significance calculations.
	 */
	private int permutationCount = 1000;

	private String tcgaDirectory;

	/**
	 * The filename that contains mutation effects.
	 */
	private String mutationEffectFilename;

	/**
	 * The change value where the most saturated color will be used on the network. This is supposed to be a positive
	 * value and the negative saturation value will be symmetrical.
	 */
	private double colorSaturationValue = 1;

	private Double correlationUpperThreshold = null;

	private Integer minimumSampleSize = null;

	/**
	 * For providing activity changes for the comparison.
	 */
	private Map<String, Integer> activityMap;

	/**
	 * A string that can indicate the combinations of the built-in network resources to use in the analysis. The
	 * possible resources are PC, REACH, PhosphoNetworks, and TRRUST. User is supposed to provide a string
	 */
	private String networkSelection;

	private boolean generateDataCentricGraph = false;

	private String tfActivityFile;

	private boolean showAllGenesWithProteomicData = false;

	private boolean hideDataNotPartOfCausalRelations = false;

	private boolean useNetworkSignificanceForCausalReasoning = false;

	private int minimumPotentialTargetsToConsiderForDownstreamSignificance = 5;

	private boolean showInsignificantData = false;

	private boolean testMissingValues = false;

	private double missingValueTestDataSufficiencyThreshold = 0.001;

	private String randomizedMatrixDirectory;

	private RandomizedMatrices phosphoRM;
	private RandomizedMatrices totProtRM;

	/**
	 * The search engine.
	 */
	CausalitySearcher cs;

	public static void main(String[] args) throws IOException, ClassNotFoundException
	{
		if (args.length < 1)
		{
			System.out.println("Please specify the data directory that contains the \"" + PARAMETER_FILENAME + "\" " +
				"file, as the first and only program parameter. All other parameters must go into the parameters " +
				"file. Below are possible parameters in the parameters file.\n\n" + Parameter.getUsageInfo());
			return;
		}

		if (args[0].equals("params-in-json"))
		{
			Map map = Parameter.getParamsInfoAsJson();
			BufferedWriter writer = Files.newBufferedWriter(Paths.get("parameter-info.json"));
			JsonUtils.writePrettyPrint(writer, map);
			writer.close();
			return;
		}

		new CausalPath(args[0]).run();
	}

	public CausalPath(String directory) throws IOException
	{
		this.directory = directory;
		this.cs = new CausalitySearcher(true);

		readParameters();
	}

	private void readParameters() throws IOException
	{
		Files.lines(Paths.get(directory + File.separator + PARAMETER_FILENAME)).filter(l -> !l.trim().startsWith("#"))
			.filter(l -> !l.isEmpty())
			.map(l -> new String[]{l.substring(0, l.indexOf("=")).trim(), l.substring(l.indexOf("=") + 1).trim()})
			.forEach(t -> setParameter(t[0], t[1]));
	}

	private void setParameter(String key, String value)
	{
		Parameter param = Parameter.findEnum(key);

		if (param == null) throw new RuntimeException("Unknown parameter: " + key);

		param.reader.read(value, this);
	}

	/**
	 * Executes the analysis;
	 */
	public void run() throws IOException, ClassNotFoundException
	{
		System.out.println("directory = " + directory);

		// If there is no platform file, use the values file instead.
		if (proteomicsPlatformFile == null) proteomicsPlatformFile = proteomicsValuesFile;

		// Read platform file
		List<ProteomicsFileRow> rows = ProteomicsFileReader.readAnnotation(
			directory + File.separator + proteomicsPlatformFile,
			IDColumn, symbolsColumn, sitesColumn, effectColumn);

		// Marking control and test just in case needed.
		boolean[] ctrl = null;
		boolean[] test = null;

		// Read values
		List<String> vals = new ArrayList<>();
		if (transformation.isTwoGroupComparison())
		{
			vals.addAll(controlValueColumn);
			vals.addAll(testValueColumn);
			ctrl = new boolean[vals.size()];
			test = new boolean[vals.size()];
			for (int i = 0; i < controlValueColumn.size(); i++)
			{
				ctrl[i] = true;
			}
			for (int i = controlValueColumn.size(); i < vals.size(); i++)
			{
				test[i] = true;
			}
		}
		else vals.addAll(valueColumn);

		ProteomicsFileReader.addValues(rows, directory + File.separator + proteomicsValuesFile,
			IDColumn, vals, defaultMissingValue, doLogTransfrorm);

		// Add activity changes from a tf activity analysis
		readTFActivityFile(rows);

		// Add activity changes from parameters file
		addActivityChangesFromParametersFile(rows);

		// Fill-in missing effects
		SiteEffectCollective sec = new SiteEffectCollective();
		sec.fillInMissingEffect(rows, siteEffectProximityThreshold);

		ProteomicsLoader loader = new ProteomicsLoader(rows, stDevThresholds);
		if (testMissingValues) loader.initMissingDataForProteins();
//		loader.printStDevHistograms();

		// Load signed relations
		Set<Relation> relations = networkSelection == null ? NetworkLoader.load() :
			NetworkLoader.load(NetworkLoader.ResourceType.getSelectedResources(networkSelection));

		if (cs.getGraphFilter() != null) relations = cs.getGraphFilter().preAnalysisFilter(relations);

		System.out.println("Number of relations that go into analysis = " + relations.size());
		for (RelationType type : relations.stream().map(r -> r.type).collect(Collectors.toSet()))
		{
			System.out.println(type + " = " + relations.stream().filter(r -> r.type.equals(type)).count());
		}

		// Associate relations with the data
		loader.decorateRelations(relations);

		// Mark some decisions
		boolean useCorrelation = transformation == ValueTransformation.CORRELATION;
		boolean controlFDR = (transformation == ValueTransformation.SIGNIFICANT_CHANGE_OF_MEAN && fdrThresholdForDataSignificance != null) ||
			(transformation == ValueTransformation.CORRELATION && fdrThresholdForCorrelation > 0);

		Set<String> dataIDs = null;
		if (randomizedMatrixDirectory != null)
		{
			// Collect all proteomic data IDs associated with
			dataIDs = relations.stream()
				.map(r -> new GeneWithData[]{r.sourceData, r.targetData}).flatMap(Arrays::stream)
				.map(gwd -> gwd.getData(DataType.PROTEIN, DataType.PHOSPHOPROTEIN)).flatMap(Collection::stream)
				.map(ExperimentData::getId).collect(Collectors.toSet());
		}

		// Init correlation detector if needed
		CorrelationDetector corrDet = null;
		if (useCorrelation)
		{
			corrDet = new CorrelationDetector(controlFDR ? -1 : correlationValueThreshold, controlFDR ? -1 : pvalThresholdForCorrelation);
			if (minimumSampleSize != null) corrDet.setMinimumSampleSize(minimumSampleSize);
			if (correlationUpperThreshold != null) corrDet.setCorrelationUpperThreshold(correlationUpperThreshold);
			corrDet.setUseMissingData(testMissingValues);
			corrDet.setCategDataSufficiencyThreshold(missingValueTestDataSufficiencyThreshold);
			if (randomizedMatrixDirectory != null)
			{
				loadRandomMatrices(dataIDs);
				corrDet.setRandomMatrices(phosphoRM, totProtRM, vals);
			}

			for (Relation rel : relations)
			{
				rel.setChDet(corrDet);
			}
		}
		else
		{
			// Associate change detectors
//			loader.associateChangeDetector(getOneDataChangeDetector(ctrl, test), data -> data instanceof ProteinData);
			if (hasThresholdFor(DataType.PHOSPHOPROTEIN))
			{
				loader.associateChangeDetector(getOneDataChangeDetector(DataType.PHOSPHOPROTEIN, ctrl, test, dataIDs),
					data -> data instanceof PhosphoProteinData);
			}
			if (hasThresholdFor(DataType.PROTEIN))
			{
				loader.associateChangeDetector(getOneDataChangeDetector(DataType.PROTEIN, ctrl, test, dataIDs),
					data -> data instanceof ProteinData && !(data instanceof PhosphoProteinData));
			}

			// Revisit this line after adding activity file support - todo
			loader.associateChangeDetector(new ThresholdDetector(0.1, ThresholdDetector.AveragingMethod.ARITHMETIC_MEAN), data -> data instanceof ActivityData);
		}

		// Load other TCGA profiles if available
		loadOtherAvailableTCGAProfiles(ctrl, test, vals, relations);

		//---DEBUG
//		PrepareBoxPlots.runWithRelations(relations, directory, "D3", "M3");
//		ReportDifferentiallyExpressed.report(relations, directory + "/" + "differentially-expressed-proteins.txt");
		//---END OF DEBUG

		// Write down the value changes
		if (!useCorrelation)
		{
			writeValueChanges(relations);
		}

		//---DEBUG
//		RobustnessAnalysis ra = new RobustnessAnalysis(cs, relations, fdrThresholdForDataSignificance, 0.05,
//			useCorrelation, fdrThresholdForCorrelation, corrDet);
//		ra.run(1000, directory + "/robustness.txt");
//		System.exit(0);
		//---END OF DEBUG

		// Search causal or conflicting relations
		Set<Relation> causal = cs.run(relations);

		if (controlFDR)
		{
			adjustPvalThresholdToFDR(relations, useCorrelation, corrDet, cs.copy(), cs.getDataUsedForInference(),
				cs.getPairsUsedForInference(), causal);
			causal = cs.run(relations);
		}

//		loader.printStDevHistograms(cs.getDataUsedForInference());

		// Significance calculation
		NetworkSignificanceCalculator nsc = calculateNetworkSignificance(relations, useCorrelation, cs.copy());

		// Add network significance as data if opted for

		if (nsc != null && !useCorrelation && useNetworkSignificanceForCausalReasoning)
		{
			if (addNetworkSignificanceAsData(relations, (NSCForNonCorr) nsc))
			{
				// Run the inference again with new activity data
				causal = cs.run(relations);
			}
		}

		int causativeSize = causal.size();
		System.out.println("Causative relations = " + causativeSize);

		GraphWriter writer = new GraphWriter(causal, nsc);
		writer.setUseGeneBGForTotalProtein(!useCorrelation);
		writer.setColorSaturationValue(colorSaturationValue);
		writer.setShowInsignificantData(showInsignificantData);

		if (useCorrelation)
		{
			writer.setExperimentDataToDraw(cs.getPairsUsedForInference().stream().flatMap(Collection::stream)
				.collect(Collectors.toSet()));
		}
		else
		{
			if (showAllGenesWithProteomicData)
			{
				Set<GeneWithData> set = relations.stream().map(r -> r.sourceData).filter(GeneWithData::hasProteomicData).collect(Collectors.toSet());
				relations.stream().map(r -> r.targetData).filter(GeneWithData::hasProteomicData).forEach(set::add);
				writer.setOtherGenesToShow(set);
			}
			else if (hideDataNotPartOfCausalRelations)
			{
				writer.setExperimentDataToDraw(cs.getDataUsedForInference());
			}
		}

		// Generate output
		writer.writeSIFGeneCentric(directory + File.separator + CAUSATIVE_RESULT_FILE_PREFIX);
		writer.writeJSON(directory + File.separator + CAUSATIVE_RESULT_FILE_PREFIX);
		if (generateDataCentricGraph)
		{
			writer.writeSIFDataCentric(directory + File.separator + CAUSATIVE_RESULT_FILE_DATA_CENTRIC_PREFIX,
				cs.getInferenceUnits());
		}

		// Note the sites with unknown effect whose determination will improve the results
		writeSitesToCurate(cs.getDataNeedsAnnotation());

		// Do the same for conflicting relations

		cs.setCausal(false);
		Set<Relation> conflicting =  cs.run(relations);
		int conflictSize = conflicting.size();
		System.out.println("Conflicting relations = " + conflictSize);

		writer = new GraphWriter(conflicting, null);
		writer.setUseGeneBGForTotalProtein(!useCorrelation);
		writer.setColorSaturationValue(colorSaturationValue);
		if (!showInsignificantData) writer.setExperimentDataToDraw(cs.getDataUsedForInference());
		if (useCorrelation)
		{
			writer.setExperimentDataToDraw(cs.getPairsUsedForInference().stream().flatMap(Collection::stream)
				.collect(Collectors.toSet()));
		}
		writer.writeSIFGeneCentric(directory + File.separator + CONFLICTING_RESULT_FILE_PREFIX);
		writer.writeJSON(directory + File.separator + CONFLICTING_RESULT_FILE_PREFIX);

		// Report conflict/causal ratio
		if (causativeSize > 0)
		{
			System.out.println("conflict / causative ratio = " + conflictSize / (double) causativeSize);
		}

		// Estimate accuracy if did a network propagation of data
		PropagationAccuracyPredictor pred = new PropagationAccuracyPredictor();
		double accuracy = pred.run(causal, conflicting, cs.copy());
		System.out.println("accuracy = " + accuracy);
	}

	private boolean hasThresholdFor(DataType type)
	{
		if (thresholdForDataSignificance != null && thresholdForDataSignificance.containsKey(type)) return true;
		if (fdrThresholdForDataSignificance != null && fdrThresholdForDataSignificance.containsKey(type)) return true;
		return false;
	}

	/**
	 * Based on network significances, adds activation or inhibition data.
	 * @param relations relation in the analysis
	 * @param nsc the network significance calculator for non-correlation cases
	 * @return true if any data is added
	 */
	public boolean addNetworkSignificanceAsData(Set<Relation> relations, NSCForNonCorr nsc)
	{
		Map<String, Integer> geneMap = nsc.getSignificantGenes();
		Set<String> genes = geneMap.keySet();
		if (!geneMap.isEmpty())
		{
			Map<String, GeneWithData> idToGene = relations.stream().map(r -> r.sourceData).filter(g -> genes.contains(g.getId())).distinct()
				.collect(Collectors.toMap(GeneWithData::getId, g -> g));

			OneDataChangeDetector chDet = new ThresholdDetector(
				0.01, ThresholdDetector.AveragingMethod.ARITHMETIC_MEAN);

			for (String id : genes)
			{
				GeneWithData gene = idToGene.get(id);

				// If there is no data change on the gene, or if there is already an activity data associated with
				// this gene, skip it
				if (gene.getChangedData().isEmpty() || !gene.getData(DataType.ACTIVITY).isEmpty()) continue;

				if (geneMap.get(id) == 1 || geneMap.get(id) == 0)
				{
					ActivityData data = new ActivityData(id + "-active-by-network-sig", id);
					data.data = new SingleCategoricalData[]{new Activity(1)};
					data.setChDet(chDet);
					gene.add(data);
				}
				if (geneMap.get(id) == -1 || geneMap.get(id) == 0)
				{
					ActivityData data = new ActivityData(id + "-inactive-by-network-sig", id);
					data.data = new SingleCategoricalData[]{new Activity(-1)};
					data.setChDet(chDet);
					gene.add(data);
				}
			}
			return true;
		}
		return false;
	}

	public NetworkSignificanceCalculator calculateNetworkSignificance(Set<Relation> relations, boolean useCorrelation,
		CausalitySearcher cs) throws IOException
	{
		NetworkSignificanceCalculator nsc = null;

		if (calculateNetworkSignificance)
		{
			if (useCorrelation)
			{
				nsc = new NSCForCorrelation(relations, cs);
			}
			else
			{
				nsc = new NSCForNonCorr(relations, cs);
			}

			String outFile = directory + File.separator + SIGNIFICANCE_FILENAME;

			if (Files.exists(Paths.get(outFile)))
			{
				nsc.loadFromFile(outFile);
			}
			else
			{
				nsc.setMinimumPotentialTargetsToConsider(minimumPotentialTargetsToConsiderForDownstreamSignificance);
				nsc.run(permutationCount);
				nsc.writeResults(outFile);
			}

			nsc.setFDRThreshold(fdrThresholdForNetworkSignificance);
			System.out.println("Graph size pval = " + nsc.getOverallGraphSizePval());
		}
		return nsc;
	}

	public void adjustPvalThresholdToFDR(Set<Relation> relations, boolean useCorrelation, CorrelationDetector corrDet,
		CausalitySearcher cs, Set<ExperimentData> datas, Set<Set<ExperimentData>> pairs, Set<Relation> relsfromFirstRun)
		throws IOException
	{

		datas = new HashSet<>(datas);
		pairs = new HashSet<>(pairs);

		cs.setCausal(false);
		cs.setCollectDataUsedForInference(true);
		Set<Relation> testedRels = cs.run(relations);
		testedRels.addAll(relsfromFirstRun);

		// DEBUG---------------
		System.out.println("Size of relations actually tested = " + testedRels.size());
		saveRels(testedRels);
		// DEBUG---------------

		datas.addAll(cs.getDataUsedForInference());
		pairs.addAll(cs.getPairsUsedForInference());
		cs.setCausal(true);
		if (useCorrelation)
		{
			FDRAdjusterForCorrelation fad = new FDRAdjusterForCorrelation(directory, pairs, corrDet);
			fad.adjustPValueThresholdsForFDR(fdrThresholdForCorrelation);
		}
		else
		{
			FDRAdjuster fad = new FDRAdjuster(directory, poolProteomicsForFDRAdjustment);
			fad.adjustPValueThresholdsOfDatas(datas, fdrThresholdForDataSignificance);

			// fdr adjust other data types on the nodes
			Set<DataType> selectiveTypes = datas.stream().map(ExperimentData::getType).collect(Collectors.toSet());
			Set<DataType> otherTypes = fdrThresholdForDataSignificance.keySet().stream()
				.filter(t -> !selectiveTypes.contains(t)).collect(Collectors.toSet());

			if (!otherTypes.isEmpty())
			{
				Set<ExperimentData> otherData = testedRels.stream().map(Relation::getAllData).flatMap(Collection::stream)
					.collect(Collectors.toSet());
				otherData = otherData.stream().filter(d -> otherTypes.contains(d.getType())).collect(Collectors.toSet());

				fad.adjustPValueThresholdsOfDatas(otherData, fdrThresholdForDataSignificance);
			}
		}
	}

	private void saveRels(Set<Relation> rels)
	{
		try
		{
			BufferedWriter writer = Files.newBufferedWriter(Paths.get("temp.sif"));
			rels.forEach(r -> FileUtil.writeln(ArrayUtil.getString("\t", r.source, r.type.name, r.target), writer));
			writer.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public void loadOtherAvailableTCGAProfiles(boolean[] ctrl, boolean[] test, List<String> vals, Set<Relation> relations) throws IOException, ClassNotFoundException
	{
		if (tcgaDirectory != null)
		{
			TCGALoader tcga = new TCGALoader(directory + File.separator + tcgaDirectory);
			tcga.setSamples(vals.toArray(new String[vals.size()]));

			if (mutationEffectFilename != null)
			{
				tcga.loadMutationEffectMap(directory + File.separator + mutationEffectFilename);
			}

			tcga.decorateRelations(relations);

			tcga.associateChangeDetector(getOneDataChangeDetector(DataType.CNA, ctrl, test, null), data -> data instanceof CNAData);
			tcga.associateChangeDetector(getOneDataChangeDetector(DataType.RNA, ctrl, test, null), data -> data instanceof RNAData);
			tcga.associateChangeDetector(getOneDataChangeDetector(DataType.MUTATION, ctrl, test, null), data -> data instanceof MutationData);
		}
	}

	public void addActivityChangesFromParametersFile(List<ProteomicsFileRow> rows)
	{
		if (activityMap != null)
		{
			for (String gene : activityMap.keySet())
			{
				int a = activityMap.get(gene);
				ProteomicsFileRow activityRow = new ProteomicsFileRow(gene + "-" + (a < 0 ? "inactive" : "active"),
					null, Collections.singletonList(gene), null);
				activityRow.makeActivityNode(a > 0);
				rows.add(activityRow);
			}
		}
	}

	/**
	 * Reads in the transcription factor activity file if present.
	 * @param rows proteomics data loaded so far
	 * @throws IOException
	 */
	public void readTFActivityFile(List<ProteomicsFileRow> rows) throws IOException
	{
		if (tfActivityFile != null)
		{
			if (!tfActivityFile.startsWith("/")) tfActivityFile = directory + "/" + tfActivityFile;

			Files.lines(Paths.get(tfActivityFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
			{
				boolean a = t[1].startsWith("a");
				ProteomicsFileRow activityRow = new ProteomicsFileRow(t[0] + "-" + (a ? "active-tf" : "inactive-tf"),
					null, Collections.singletonList(t[0]), null);
				activityRow.makeActivityNode(a);
				rows.add(activityRow);
			});
		}
	}

	private OneDataChangeDetector getOneDataChangeDetector(DataType type, boolean[] ctrl, boolean[] test, Set<String> dataIDs) throws IOException, ClassNotFoundException
	{
		OneDataChangeDetector detector = null;

		if (transformation == ValueTransformation.ARITHMETIC_MEAN ||
			transformation == ValueTransformation.GEOMETRIC_MEAN ||
			transformation == ValueTransformation.MAX)
		{
			detector = new ThresholdDetector(thresholdForDataSignificance.get(type));
			((ThresholdDetector) detector).setAveragingMethod(transformation == ValueTransformation.ARITHMETIC_MEAN ?
				ThresholdDetector.AveragingMethod.ARITHMETIC_MEAN : transformation == ValueTransformation.GEOMETRIC_MEAN ?
				ThresholdDetector.AveragingMethod.FOLD_CHANGE_MEAN : ThresholdDetector.AveragingMethod.MAX);
		}
		else if (transformation == ValueTransformation.DIFFERENCE_OF_MEANS)
		{
			detector = new DifferenceDetector(thresholdForDataSignificance.get(type), ctrl, test);
		}
		else if (transformation == ValueTransformation.FOLD_CHANGE_OF_MEAN)
		{
			detector = new FoldChangeDetector(thresholdForDataSignificance.get(type), ctrl, test);
		}
		else if (transformation == ValueTransformation.SIGNIFICANT_CHANGE_OF_MEAN)
		{
			// if there is no fdr control, then use the given data significance threshold. Otherwise use 1, which means
			// make everything significant. And FDR correction will be applied later.
			detector = new SignificanceDetector(fdrThresholdForDataSignificance == null ?
				thresholdForDataSignificance.get(type) : 1, ctrl, test);

			if (type == DataType.PROTEIN || type == DataType.PHOSPHOPROTEIN)
			{
				((SignificanceDetector) detector).setUseMissingData(testMissingValues);
				((SignificanceDetector) detector).setCategDataSufficiencyThreshold(missingValueTestDataSufficiencyThreshold);

				if (minimumSampleSize != null)
				{
					((SignificanceDetector) detector).setMinimumSampleSize(minimumSampleSize);
				}

				if (randomizedMatrixDirectory != null)
				{
					if (phosphoRM == null)
					{
						loadRandomMatrices(dataIDs);
					}

					((SignificanceDetector) detector).setRadomizedMatrices(phosphoRM, totProtRM, new HashSet<>(testValueColumn));
				}
			}
		}
		return detector;
	}

	private void loadRandomMatrices(Set<String> dataIDs) throws IOException, ClassNotFoundException
	{
		BooleanMatrixRandomizer bmr = new BooleanMatrixRandomizer();
		totProtRM = bmr.readRandomMatrices(randomizedMatrixDirectory + "/" + BooleanMatrixRandomizer.TOTAL_PROT_DIR, dataIDs);
		phosphoRM = bmr.readRandomMatrices(randomizedMatrixDirectory + "/" + BooleanMatrixRandomizer.PHOSPHO_DIR, dataIDs);
	}

	private void writeSitesToCurate(Set<PhosphoProteinData> datas) throws IOException
	{
		Set<String> sites = new HashSet<>();
		for (PhosphoProteinData data : datas)
		{
			sites.addAll(data.getGenesWithSites().stream().collect(Collectors.toList()));
		}

		BufferedWriter writer = Files.newBufferedWriter(
			Paths.get(directory + File.separator + UNKNOWN_SITE_EFFECT_FILENAME));

		sites.stream().sorted().forEach(s -> FileUtil.writeln(s, writer));

		writer.close();
	}

	private void writeValueChanges(Set<Relation> relations) throws IOException
	{
		// collect the experiment data from relations
		Set<ExperimentData> datas = relations.stream().map(Relation::getAllData).flatMap(Collection::stream)
			.filter(ExperimentData::hasChangeDetector).collect(Collectors.toSet());

		// find the data classes
		Set<Class<? extends ExperimentData>> classes = datas.stream().map(d -> d.getClass()).collect(Collectors.toSet());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(directory + File.separator + VALUE_CHANGES_FILENAME));

		for (Class<? extends ExperimentData> clazz : classes)
		{
			writer.write("\n\nData type: " + clazz.getName());

			OneDataChangeDetector chDet = datas.stream().filter(d -> d.getClass().equals(clazz)).findAny().get().getChDet();

			if (chDet instanceof SignificanceDetector)
			{
				SignificanceDetector sigDet = (SignificanceDetector) chDet;
				if (minimumSampleSize != null)
				{
					sigDet.setMinimumSampleSize(minimumSampleSize);
				}
				if (testMissingValues)
				{
					sigDet.setUseMissingData(testMissingValues);
					sigDet.setCategDataSufficiencyThreshold(missingValueTestDataSufficiencyThreshold);
				}

				Map<ExperimentData, Double> pvals = datas.stream().filter(d -> d.getClass().equals(clazz))
					.collect(Collectors.toMap(d -> d, sigDet::getPValue));

				Map<ExperimentData, Double> qvals = FDR.getQVals(pvals, null);

				writer.write("\nRow ID\tChange amount\tP-value\tQ-value");
				datas.stream().filter(d -> d.getClass().equals(clazz))
					.sorted((d1, d2) -> pvals.get(d1).compareTo(pvals.get(d2))).forEach(d ->
					FileUtil.lnwrite(d.id + "\t" + d.getChangeValue() + "\t" + pvals.get(d) + "\t" + qvals.get(d),
						writer));
			}
			else
			{
				writer.write("\nRow ID\tChange amount");
				datas.stream().filter(d -> d.getClass().equals(clazz)).forEach(d ->
					FileUtil.lnwrite(d.id + "\t" + d.getChangeValue(), writer));
			}
		}

		writer.close();
	}

	/**
	 * Enumeration of options for how to use the value columns to determine significance of a change.
	 */
	enum ValueTransformation
	{
		ARITHMETIC_MEAN("arithmetic-mean", "The arithmetic mean value of the given values is used for significance " +
			"detection of a single change. There should only be one group of values (marked with value-column), the " +
			"values have to be distributed around zero, and a threshold value should be provided for significance " +
			"detection, using the threshold-for-data-significance.", false),

		GEOMETRIC_MEAN("geometric-mean", "The geometric mean value of the given values is used for significance " +
			"detection of a single change. This is the only case when the geometric mean is used for averaging a " +
			"group of samples, and it is appropriate if the individual values are formed of some kind of ratios. " +
			"There should only be one group of values (marked with value-column), the " +
			"values have to be distributed around zero, and a threshold value should be provided for significance " +
			"detection, using the threshold-for-data-significance.", false),

		MAX("max", "The value with maximum absolute is used for the analysis. There should only be one group of " +
			"values (marked with value-column), the values have to be distributed around " +
			"zero, and a threshold value should be provided for significance detection, using the " +
			"threshold-for-data-significance.", false),

		DIFFERENCE_OF_MEANS("difference-of-means", "There should be control and test values, whose difference would " +
			"be used for significance detection. The threshold for significance (threshold-for-data-significance) " +
			"should also be provided.", true),

		FOLD_CHANGE_OF_MEAN("fold-change-of-mean", "There should be control and test values, whose ratio will be " +
			"converted to fold change and thresholded. The fold change value will be in the range (-inf, -1] + [1, " +
			"inf). If the data file already contains a fold-change value, then please use the " + GEOMETRIC_MEAN.name +
			" as value transformation. The threshold for significance (threshold-for-data-significance) should also " +
			"be provided.", true),

		SIGNIFICANT_CHANGE_OF_MEAN("significant-change-of-mean", "There should be sufficient amount of control and " +
			"test values to detect the significance of change with a t-test. Technically there should be more than 3" +
			" controls and 3 tests, practically, they should be much more to provide statistical power. The " +
			"threshold-for-data-significance should be used for a p-value threshold, or " +
			"alternatively, fdr-threshold-for-data-significance should be used for " +
			"controlling significance at the false discovery rate level.", true),

		CORRELATION("correlation", "There should be one group of values (marked with value-column). There must be at " +
			"least 3 value columns technically, but many more " +
			"than that practically to have some statistical power for significant correlation. ", false);

		ValueTransformation(String name, String description, boolean twoGroupComparison)
		{
			this.name = name;
			this.description = description;
			this.twoGroupComparison = twoGroupComparison;
		}

		String name;
		String description;

		/**
		 * Whether there is a control and a test group in the analysis.
		 */
		private boolean twoGroupComparison;

		public static ValueTransformation fetch(String name)
		{
			for (ValueTransformation trans : values())
			{
				if (trans.name.equals(name)) return trans;
			}
			return null;
		}

		public boolean isTwoGroupComparison()
		{
			return twoGroupComparison;
		}

		public static String getUsageInfo()
		{
			StringBuilder sb = new StringBuilder();
			for (ValueTransformation trans : values())
			{
				sb.append("\t").append(trans.name).append(": ").append(trans.description).append("\n");
			}
			return sb.toString();
		}

		public static Map getValuesAsJson()
		{
			List list = new ArrayList<>();
			for (ValueTransformation transformation : values())
			{
				list.add(transformation.name);
			}

			Map map = new LinkedHashMap<>();
			map.put("name", "ValueTransformation");
			map.put("values", list);
			return map;
		}
	}

	enum Parameter
	{
		PROTEOMICS_PLATFORM_FILE((value, cp) -> cp.proteomicsPlatformFile = value,
			"Proteomics platform file",
			"Name of the proteomics platform file. Each row should belong to either a gene's total protein " +
				"measurement, or a site specific measurement. This file should contain ID, gene symbols, modification" +
				" sites, and known site effects. Platform file and values file can be the same file.",
			new EntryType(File.class), null, false, false, null),
		PROTEOMICS_VALUES_FILE((value, cp) -> cp.proteomicsValuesFile = value,
			"Proteomics values file",
			"Name of the proteomics values file. It should have at least one ID column and one or more columns for " +
				"experiment values. Platform file and values file can be the same file.",
			new EntryType(File.class), null, true, false, null),
		ID_COLUMN((value, cp) -> cp.IDColumn = value,
			"ID column in data file",
			"The name of the ID column in platform and values files.",
			new EntryType(String.class), new String[][]{{"ID"}}, true, false, null),
		SYMBOLS_COLUMN((value, cp) -> cp.symbolsColumn = value,
			"Symbols column in data file",
			"The name of the symbols column in platform file.",
			new EntryType(String.class), new String[][]{{"Symbols"}}, true, false, null),
		SITES_COLUMN((value, cp) -> cp.sitesColumn = value,
			"Sites column in data file",
			"The name of the sites column.",
			new EntryType(String.class), new String[][]{{"Sites"}}, true, false, null),
		EFFECT_COLUMN((value, cp) -> cp.effectColumn = value,
			"Site effect column in data file",
			"The name of the effect column.",
			new EntryType(String.class), new String[][]{{"Effect"}}, false, false, null),
		VALUE_TRANSFORMATION((value, cp) -> cp.transformation = ValueTransformation.fetch(value),
			"How to use values in the analysis",
			"This parameter determines how to use the values in the proteomics file. Possible values are listed " +
				"below.\n" + ValueTransformation.getUsageInfo(),
			new EntryType(ValueTransformation.class), null, true, false, null),
		VALUE_COLUMN((value, cp) ->
		{
			if (cp.valueColumn == null) cp.valueColumn = new ArrayList<>();
			cp.valueColumn.add(value);
		},
			"Value column in the data file",
			"Name of a value column. This parameter should be used when there is only one group of experiments to " +
				"consider in the analysis.",
			new EntryType(String.class), null, true, true,
			new Cond(Logical.OR,
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.ARITHMETIC_MEAN.name),
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.GEOMETRIC_MEAN.name),
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.MAX.name),
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name))),
		CONTROL_VALUE_COLUMN((value, cp) ->
		{
			if (cp.controlValueColumn == null) cp.controlValueColumn = new ArrayList<>();
			cp.controlValueColumn.add(value);
		},
			"Control value column in the data file",
			"Name of a control value column. This parameter should be used when there are control and test value " +
				"columns in the dataset.",
			new EntryType(String.class), null, true, true,
			new Cond(Logical.OR,
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.DIFFERENCE_OF_MEANS.name),
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.FOLD_CHANGE_OF_MEAN.name),
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.SIGNIFICANT_CHANGE_OF_MEAN.name))),
		TEST_VALUE_COLUMN((value, cp) ->
		{
			if (cp.testValueColumn == null) cp.testValueColumn = new ArrayList<>();
			cp.testValueColumn.add(value);
		},
			"Test value column in the data file",
			"Name of a test value column. This parameter should be used when there are control and test value " +
				"columns in the dataset.",
			new EntryType(String.class), null, true, true,
			new Cond(Logical.OR,
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.DIFFERENCE_OF_MEANS.name),
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.FOLD_CHANGE_OF_MEAN.name),
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.SIGNIFICANT_CHANGE_OF_MEAN.name))),
		DO_LOG_TRANSFORM((value, cp) -> cp.doLogTransfrorm = Boolean.valueOf(value),
			"Log transform data values",
			"Whether the proteomic values should be log transformed for the analysis. Possible values are 'true' and " +
				"'false'.",
			new EntryType(Boolean.class), new String[][]{{"false"}}, false, false, null),
		THRESHOLD_FOR_DATA_SIGNIFICANCE((value, cp) ->
		{
			if (cp.thresholdForDataSignificance == null) cp.thresholdForDataSignificance = new HashMap<>();
			String[] s = value.split("\\s+");
			cp.thresholdForDataSignificance.put(DataType.get(s[1]), Double.valueOf(s[0]));
		},
			"Threshold value for significant data",
			"A threshold value for selecting significant data. Use this parameter only when FDR controlling procedure" +
				"is already performed outside of CausalPath.",
			new EntryType(Double.class, DataType.class), null, false, true,
			new Cond(Logical.AND,
				new Cond(Logical.NOT, new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name)),
				new Cond("fdr-threshold-for-data-significance", null))), // typed the parameter out because of illegal forward reference
		FDR_THRESHOLD_FOR_DATA_SIGNIFICANCE((value, cp) ->
		{
			if (cp.fdrThresholdForDataSignificance == null) cp.fdrThresholdForDataSignificance = new HashMap<>();
			String[] s = value.split("\\s+");
			cp.fdrThresholdForDataSignificance.put(DataType.get(s[1]), Double.valueOf(s[0]));
		},
			"FDR threshold for data significance",
			"False discovery rate threshold for data significance. This parameter can be set for each " +
				"different data type separately. The parameter value has to be in the form 'fdr-val data-type', such like " +
				"'0.1 phosphoprotein'.",
			new EntryType(Double.class, DataType.class), null, false, true,
			new Cond(Logical.AND,
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.SIGNIFICANT_CHANGE_OF_MEAN.name),
				new Cond(THRESHOLD_FOR_DATA_SIGNIFICANCE.getText(), null))),
		STDEV_THRESHOLD_FOR_DATA((value, cp) ->
		{
			if (cp.stDevThresholds == null) cp.stDevThresholds = new HashMap<>();
			String[] s = value.split("\\s+");
			cp.stDevThresholds.put(DataType.get(s[1]), Double.valueOf(s[0]));
		},
			"Standard deviation threshold for data",
			"This parameter can be set for each different data type separately. The parameter value has to be in the" +
				" form 'stdev-thr data-type', such like '0.5 phosphoprotein'.",
			new EntryType(Double.class, DataType.class), null, false, true,
			new Cond(Logical.OR,
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name),
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.SIGNIFICANT_CHANGE_OF_MEAN.name))),
		POOL_PROTEOMICS_FOR_FDR_ADJUSTMENT((value, cp) -> cp.poolProteomicsForFDRAdjustment = Boolean.valueOf(value),
			"Pool proteomics data for FDR adjustment",
			"Whether to consider proteomic and phosphoproteomic data as a single dataset during FDR adjustment. This " +
				"is typically the case with RPPA data, and typically not the case with mass spectrometry data. Can be" +
				" 'true' or 'false'. Default is false.",
			new EntryType(Boolean.class), new String[][]{{"true"}}, true, false,
			new Cond(Logical.NOT, new Cond(FDR_THRESHOLD_FOR_DATA_SIGNIFICANCE.getText(), null))),
		CALCULATE_NETWORK_SIGNIFICANCE((value, cp) -> cp.calculateNetworkSignificance = Boolean.valueOf(value),
			"Calculate network significance",
			"Whether to calculate significances of the properties of the graph. When turned on, a p-value for network" +
				" size, and also downstream activity enrichment p-values for each gene on the graph are calculated.",
			new EntryType(Boolean.class), new String[][]{{"false"}}, true, false, null),
		PERMUTATIONS_FOR_SIGNIFICANCE((value, cp) -> cp.permutationCount = Integer.valueOf(value),
			"Number of permutations for calculating network significance",
			"We will do data randomization to see if the result network is large, or any protein's downstream is " +
				"enriched. This parameter indicates the number of randomizations we should perform. It should be " +
				"reasonable high, such as 1000, but not too high.",
			new EntryType(Integer.class), new String[][]{{"1000"}}, true, false,
			new Cond(CALCULATE_NETWORK_SIGNIFICANCE.getText(), "true")),
		FDR_THRESHOLD_FOR_NETWORK_SIGNIFICANCE((value, cp) ->
			cp.fdrThresholdForNetworkSignificance = Double.valueOf(value),
			"FDR threshold for network significance",
			"The false discovery rate for network significance calculations for the downstream activity enrichment of" +
				" genes.",
			new EntryType(Double.class), new String[][]{{"0.1"}}, true, false,
			new Cond(CALCULATE_NETWORK_SIGNIFICANCE.getText(), "true")),
		CORRELATION_VALUE_THRESHOLD((value, cp) ->
			cp.correlationValueThreshold = Double.valueOf(value),
			"Threshold for correlation value",
			"Option to control correlation with its value. This cannot be used with FDR control, but can be used with" +
				" p-value control.",
			new EntryType(Double.class), null, false, false,
			new Cond(Logical.AND,
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name),
				new Cond("fdr-threshold-for-correlation", null))), // typing out the parameter name due to illegal forward reference
		PVAL_THRESHOLD_FOR_CORRELATION((value, cp) -> cp.pvalThresholdForCorrelation = Double.valueOf(value),
			"P-value threshold for correlation",
			"A p-value threshold for correlation in a correlation-based causality. This parameter should only be used" +
				" when FDR control is performed outside of CausalPath.",
			new EntryType(Double.class), null, false, false,
			new Cond(Logical.AND,
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name),
				new Cond("fdr-threshold-for-correlation", null))), // typing out the parameter name due to illegal forward reference
		FDR_THRESHOLD_FOR_CORRELATION((value, cp) -> cp.fdrThresholdForCorrelation = Double.valueOf(value),
			"FDR threshold for correlation",
			"False discovery rate threshold for the correlations in a correlation-based analysis.",
			new EntryType(Double.class), new String[][]{{"0.01"}}, false, false,
			new Cond(Logical.AND,
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name),
				new Cond(PVAL_THRESHOLD_FOR_CORRELATION.getText(), null),
				new Cond(CORRELATION_VALUE_THRESHOLD.getText(), null))),
		CUSTOM_RESOURCE_DIRECTORY((value, cp) -> ResourceDirectory.set(value),
			"Custom resource directory name",
			"CausalPath downloads some data in the first run and stores in the resource directory. This directory is " +
				"'.panda' by default. If this needs to be customized, use this parameter.",
			new EntryType(String.class), null, false, false, new Cond(Logical.NOT)),
		SITE_EFFECT_PROXIMITY_THRESHOLD((value, cp) -> cp.siteEffectProximityThreshold = Integer.valueOf(value),
			"Site-effect proximity threshold",
			"CausalPath has a database of phosphorylation site effects. When not set, this parameter is 0 by default," +
				" which means exact usage of site effects. But sometimes we may see a site with unknown effect is " +
				"modified, which is very close to another site with a known effect. This parameter let's us to assume" +
				" those changing sites with unknown effect has the same effect with the neighbor site with known " +
				"effect. Use responsibly.",
			new EntryType(Double.class), new String[][]{{"0"}}, true, false, null),
		DO_SITE_MATCHING((value, cp) -> cp.cs.setForceSiteMatching(Boolean.valueOf(value)),
			"Do site matching",
			"Whether to force site matching in causality analysis. True by default.",
			new EntryType(Boolean.class), new String[][]{{"true"}}, true, false, null),
		SITE_MATCH_PROXIMITY_THRESHOLD((value, cp) -> cp.cs.setSiteProximityThreshold(Integer.valueOf(value)),
			"Site-match proximity threshold",
			"Phosphorylation relations many times know the target sites. When we observe a change in a site of the " +
				"target protein which is not targeted by the relation, but the site is very close to a known target " +
				"site of the relation, this parameter let's us to assume that the relation also applies to those " +
				"close-by sites.",
			new EntryType(Double.class), new String[][]{{"0"}}, true, false, new Cond(DO_SITE_MATCHING.getText(), "true")),
		DEFAULT_MISSING_VALUE((value, cp) -> cp.defaultMissingValue = Double.valueOf(value),
			"Default missing value in proteomics file",
			"An option to specify a default value for the missing values in the proteomics file.",
			new EntryType(Double.class), null, false, false, null),
		RELATION_FILTER_TYPE((value, cp) ->
		{
			if (!cp.cs.hasGraphFilter())
			{
				cp.cs.setGraphFilter(new GraphFilter(value));
			}
			else
			{
				cp.cs.getGraphFilter().setRelationFilterType(GraphFilter.RelationFilterType.get(value));
			}
		},
			"Network relation-type filter",
			"Use this parameter to limit the results with a specific type of relation. Possible values are below.\n" +
				GraphFilter.RelationFilterType.getUsageInfo(),
			new EntryType(GraphFilter.RelationFilterType.class),
			new String[][]{{GraphFilter.RelationFilterType.NO_FILTER.getName()}},
			true, false, null),
		GENE_FOCUS((value, cp) ->
		{
			if (!cp.cs.hasGraphFilter())
			{
				cp.cs.setGraphFilter(new GraphFilter(new HashSet<>(Arrays.asList(value.split(";")))));
			}
			else
			{
				cp.cs.getGraphFilter().setFocusGenes(new HashSet<>(Arrays.asList(value.split(";"))));
			}
		},
			"Gene to focus",
			"Use this parameter to crop the result network to the neighborhood of certain gene. You should provide " +
				"gene symbols of these genes in a row separated by a semicolon, such like 'MTOR;RPS6KB1;RPS6'",
			new EntryType(String.class), null, false, false, null),
		TCGA_DIRECTORY((value, cp) -> cp.tcgaDirectory = value,
			"TCGA data directory",
			"It is possible to add genomic data from TCGA to CausalPath analysis. This is only useful when the " +
				"proteomic data have the same sample IDs. Users can load TCGA data into a local directory from Broad " +
				"Firehose, and provide the directory here. The org.panda.resource.tcga.BroadDownloader in the project" +
				" https://github.com/PathwayAndDataAnalysis/resource is a utility that can do that.",
			new EntryType(String.class), null, false, false, new Cond(Logical.NOT)),
		MUTATION_EFFECT_FILE((value, cp) -> cp.mutationEffectFilename = value,
			"Mutation effect file",
			"When we have mutations in the analysis, users can provide mutation effects using this parameter, " +
				"otherwise all mutations are assumed to be inactivating.",
			new EntryType(File.class), null, false, false, null),
		COLOR_SATURATION_VALUE((value, cp) -> cp.colorSaturationValue = Double.valueOf(value),
			"Node color saturation value",
			"Specifies the value where node colors reach most intense color. Has to be a positive value, and used " +
				"symmetrically. In the case of value-transformation is significant-change-of-mean, the value is " +
				"-log(p) with a sign associated to it.",
			new EntryType(Double.class), new String[][]{{"1"}}, false, false,
			new Cond( Logical.NOT, new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name))),
		SHOW_ALL_GENES_WITH_PROTEOMIC_DATA((value, cp) -> cp.showAllGenesWithProteomicData = Boolean.valueOf(value),
			"Show unexplained proteomic data",
			"CausalPath generates a result graph, but what about all other significant changes that could not make " +
				"into the network? CausalPath puts those genes as disconnected nodes in the graph when the analysis " +
				"is not correlation based. This is true by default but can be turned off by setting to false.",
			new EntryType(Boolean.class), new String[][]{{"true"}}, true, false,
			new Cond( Logical.NOT, new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name))),
		SHOW_INSIGNIFICANT_DATA((value, cp) ->
			cp.showInsignificantData = Boolean.valueOf(value),
			"Show insignificant proteomic data on the graph",
			"Option to make the insignificant protein data on the result graph visible. Seeing these is good for " +
				"seeing what is being measured, but when they are too much, turning off generates a a better view.",
			new EntryType(Boolean.class), new String[][]{{"false"}}, true, false,
			new Cond(Logical.NOT, new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name))),
		HIDE_DATA_NOT_PART_OF_CAUSAL_RELATIONS((value, cp) ->
			cp.hideDataNotPartOfCausalRelations = Boolean.valueOf(value),
			"Hide data which did not contribute causal relations",
			"Limits the data drawn on the result graph to the ones that take part in the identified causal relations.",
			new EntryType(Boolean.class), new String[][]{{"false"}}, true, false,
			new Cond(Logical.AND,
				new Cond(Logical.NOT, new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name)),
				new Cond(SHOW_ALL_GENES_WITH_PROTEOMIC_DATA.getText(), "false"))),
		BUILT_IN_NETWORK_RESOURCE_SELECTION((value, cp) -> cp.networkSelection = value,
			"Built-in network resources to use",
			"Determines which network resource to use during the analysis. Multiple network resource should be " +
				"mentioned together, separated with a space or comma. Possible values are below.\n" +
				NetworkLoader.ResourceType.getUsageInfo(),
			new EntryType(NetworkLoader.ResourceType.class), new String[][]{{
				NetworkLoader.ResourceType.PC.name(),
				NetworkLoader.ResourceType.PhosphoNetworks.name(),
				NetworkLoader.ResourceType.IPTMNet.name()}},
			true, true, null),
		GENERATE_DATA_CENTRIC_GRAPH((value, cp) -> cp.generateDataCentricGraph = Boolean.valueOf(value),
			"Generate a data-centric view as well",
			"An alternative to the gene-centric graph of CausalPath is a data-centric graph where nodes are not genes" +
				" but the data. This parameter forces to generate this type of result as well. False by default.",
			new EntryType(Boolean.class), new String[][]{{"false"}}, true, false, null),
		CORRELATION_UPPER_THRESHOLD((value, cp) -> cp.correlationUpperThreshold = Double.valueOf(value),
			"An upper threshold for correlation value",
			"In some types of proteomic data, highest correlations come from errors. A way around is filtering with " +
				"an upper value.",
			new EntryType(Double.class), null, false, false,
			new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name)),
		MINIMUM_SAMPLE_SIZE((value, cp) -> cp.minimumSampleSize = Integer.valueOf(value),
			"Minimum sample size",
			"When there are missing values in proteomic file, the comparisons can have different sample sizes for " +
				"controls and tests. This parameter sets the minimum sample size of the control and test sets.",
			new EntryType(Integer.class), new String[][]{{"3"}}, true, false,
			new Cond(Logical.OR,
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name),
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.SIGNIFICANT_CHANGE_OF_MEAN.name))),
		GENE_ACTIVITY((value, cp) ->
		{
			if (cp.activityMap == null) cp.activityMap = new HashMap<>();
			String[] t = value.split(" ");
			cp.activityMap.put(t[0], ActivityLabel.getLabel(t[1]) == ActivityLabel.ACTIVATED ? 1 : -1);
		},
			"Gene activity hypotheses to include",
			"Use this parameter to assign a specific activity or inactivity to a gene in the analysis. The value has " +
				"to start with a gene name and one letter code for activity or inactivity, such as 'BRAF a', or 'PTEN" +
				" i'.",
			new EntryType(String.class, ActivityLabel.class), null, false, true,
			new Cond(Logical.NOT, new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name))),
		TF_ACTIVITY_FILE((value, cp) -> cp.tfActivityFile = value,
			"Transcription factory activity inference file",
			"CausalPath lets users to input results from an inference for transcriptional factor activities, such as" +
				" PRECEPTS, MARINa, or VIPER. For this, the results should be prepared in a file, first column " +
				"containing TF symbol and the second column whether 'activated' or 'inhibited'. The name of such file" +
				" should be provided here.",
			new EntryType(File.class), null, false, false,
			new Cond(Logical.NOT, new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name))),
		USE_STRONGEST_PROTEOMIC_DATA_PER_GENE((value, cp) ->
			cp.cs.setUseStrongestProteomicsDataForActivity(Boolean.valueOf(value)),
			"Use strongest proteomic data per gene",
			"When a proteomic experiment outputs too many phosphorylation sites with lots of changes, many proteins " +
				"have evidences for both activation and inhibition. This produces networks hard to read. A complexity" +
				" management technique is to turn on this parameter to use only the strongest proteomic feature at " +
				"the upstream of relations. This is false by default.",
			new EntryType(Boolean.class), new String[][]{{"false"}}, true, false, null),
		USE_NETWORK_SIGNIFICANCE_FOR_CAUSAL_REASONING((value, cp) ->
			cp.useNetworkSignificanceForCausalReasoning = Boolean.valueOf(value),
			"Use network significance for causal reasoning",
			"After calculation of network significances in a non-correlation-based analysis, this option introduces" +
				" the detected active and inactive proteins as data to be used in the analysis. This applies only to " +
				"the proteins that already have a changed data on them, and have no previous activity data associated.",
		new EntryType(Boolean.class), new String[][]{{"false"}}, true, false,
		new Cond(Logical.NOT, new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name))),
		MINIMUM_POTENTIAL_TARGETS_TO_CONSIDER_FOR_DOWNSTREAM_SIGNIFICANCE((value, cp) ->
			cp.minimumPotentialTargetsToConsiderForDownstreamSignificance = Integer.valueOf(value),
			"Minimum potential targets to calculate network significance",
			"While calculating downstream significance for each source gene, we may not like to include those genes " +
				"with already few qualifying targets to reduce noise in data and reduce the number of tested " +
				"hypotheses.",
			new EntryType(Integer.class), new String[][]{{"5"}}, true, false,
			new Cond(CALCULATE_NETWORK_SIGNIFICANCE.getText(), "true")),
		USE_MISSING_PROTEOMIC_DATA_FOR_TEST((value, cp) ->
			cp.testMissingValues = Boolean.valueOf(value),
			"Include missing proteomic data in tests",
			"Option to use a G-test to check unequal distribution of missing values. If opted, and data is sufficient" +
				", the G-test result is combined with t-test result with Fisher's method. But beware. This method " +
				"assumes that missing values are uniformly distributed to samples. If this is violated, then false " +
				"positives will appear. If you are not sure, stay away from this option.",
			new EntryType(Boolean.class), new String[][]{{"false"}}, true, false,
			new Cond(Logical.OR,
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.SIGNIFICANT_CHANGE_OF_MEAN.name),
				new Cond(VALUE_TRANSFORMATION.getText(), ValueTransformation.CORRELATION.name))),
		RANDOMIZED_MATRIX_DIRECTORY_FOR_MISSING_PROTEOMIC_DATA((value, cp) ->
			cp.randomizedMatrixDirectory = value,
			"Directory name for randomized boolean matrices for missing data distribution",
			"Using randomization is an alternative to using a G-test for interpreting missing data distribution. " +
				"CausalPath cannot generate those matrices, but it can use pre-generated matrices to compute " +
				"significances. This operation typically requires a lot of memory.",
			new EntryType(String.class), null, false, false, new Cond(Logical.NOT)),
		MISSING_VALUE_TEST_DATA_SUFFICIENCY_THRESHOLD((value, cp) ->
			cp.missingValueTestDataSufficiencyThreshold = Double.valueOf(value),
			"Missing value test data sufficiency threshold",
			"When we use a G-test in the analysis, we don't want to use it for every proteomic row. Some rows will " +
				"have insufficient data. To test for sufficiency, we generate an extreme case where missing data is " +
				"shifted to the smaller group and see if this can provide a p-value small enough. If this extreme " +
				"shifting cannot make the the p-value small enough (specified with this threshold), then we don't use" +
				" a G-test for that row.",
			new EntryType(Double.class), new String[][]{{"0.001"}}, true, false,
			new Cond(USE_MISSING_PROTEOMIC_DATA_FOR_TEST.getText(), "true")),
		;

		ParameterReader reader;
		String title;
		String info;
		EntryType type;
		String[][] defaultValue;
		boolean mandatory;
		boolean canBeMultiple;

		/**
		 * Condition for the parameter to be effective.
		 */
		Cond condition;

		Parameter(ParameterReader reader, String title, String info, EntryType type, String[][] defaultValue,
			boolean mandatory, boolean canBeMultiple, Cond condition)
		{
			this.reader = reader;
			this.title = title;
			this.info = info;
			this.type = type;
			this.defaultValue = defaultValue;
			this.mandatory = mandatory;
			this.canBeMultiple = canBeMultiple;
			this.condition = condition;
		}

		String getText()
		{
			return toString().toLowerCase().replaceAll("_", "-");
		}

		static Parameter findEnum(String text)
		{
			if (text == null) return null;

			text = text.trim();
			for (Parameter parameter : values())
			{
				if (parameter.getText().equals(text)) return parameter;
			}
			return null;
		}

		public String getInfo()
		{
			return info;
		}

		public EntryType getType()
		{
			return type;
		}

		public List<List<String>> getDefaultValue()
		{
			if (defaultValue == null) return null;
			List<List<String>> defList = new ArrayList<>();
			for (String[] defs : defaultValue)
			{
				defList.add(Arrays.asList(defs));
			}
			return defList;
		}

		public boolean isMandatory()
		{
			return mandatory;
		}

		public Cond getCondition()
		{
			return condition;
		}

		public static String getUsageInfo()
		{
			StringBuilder sb = new StringBuilder();
			for (Parameter param : values())
			{
				sb.append(param.getText()).append(": ").append(param.getInfo()).append("\n");
			}
			return sb.toString();
		}

		public static Map getParamsInfoAsJson()
		{
			Map map = new LinkedHashMap<>();
			List list = new ArrayList<>();
			for (Parameter param : values())
			{
				if (param.condition != null && param.condition.skipAtWebServer())
				{
					continue;
				}

				list.add(param.getInfoAsJson());
			}
			map.put("Parameters", list);

			list = new ArrayList<>();
			list.add(ValueTransformation.getValuesAsJson());
			list.add(DataType.getValuesAsJson());
			list.add(NetworkLoader.ResourceType.getValuesAsJson());
			list.add(GraphFilter.RelationFilterType.getValuesAsJson());
			list.add(ActivityLabel.getValuesAsJson());

			map.put("Enumerations", list);
			return map;
		}

		public Map getInfoAsJson()
		{
			Map map = new LinkedHashMap<>();
			map.put("ID", getText());
			map.put("Title", title);
			map.put("Description", info);
			map.put("EntryType", type.getInfoAsJson());
			if (defaultValue != null && defaultValue.length > 0) map.put("Default", getDefaultValue());
			map.put("Mandatory", Boolean.toString(mandatory));
			map.put("CanBeMultiple", Boolean.toString(canBeMultiple));
			if (condition != null) map.put("Condition", condition.getAsJson());
			return map;
		}
	}

	interface ParameterReader
	{
		void read(String value, CausalPath cp);
	}

	/**
	 * Conditions for
	 */
	static class Cond
	{
		Logical op;
		Cond[] conds;
		String param;
		String[] value;

		static final List<String> nullList = new ArrayList<>(Collections.singletonList(null));

		public Cond(String param, String... value)
		{
			if (value != null && value.length == 0)
				throw new IllegalArgumentException("Value cannot be absent.");
			this.param = param;
			this.value = value;
		}

		public Cond(Logical op, Cond... conds)
		{
			if (op == null) throw new IllegalArgumentException("Use the other constructor for basic condition.");
			if (op != Logical.NOT && conds.length == 0) throw new IllegalArgumentException("Empty conditions is not allowed.");
			if (op == Logical.NOT && conds.length > 1) throw new IllegalArgumentException("NOT can take only one condition.");
			if (op != Logical.NOT && conds.length < 2) throw new IllegalArgumentException("AND and OR have to take at least 2 conditions.");

			this.op = op;
			this.conds = conds;
		}

		public Map getAsJson()
		{
			Map map = new LinkedHashMap<>();
			if (op == null)
			{
				map.put("Parameter", param);
				map.put("Value", value == null ? nullList : Arrays.asList(value));
			}
			else
			{
				map.put("Operator", op);
				List list = new ArrayList<>();
				for (Cond cond : conds)
				{
					list.add(cond.getAsJson());
				}
				map.put("Conditions", list);
			}
			return map;
		}

		public boolean skipAtWebServer()
		{
			return op == Logical.NOT && conds.length == 0;
		}
	}

	enum Logical
	{
		AND, OR, NOT
	}

	static class EntryType
	{
		Class[] types;

		public EntryType(Class... types)
		{
			this.types = types;
		}

		public List getInfoAsJson()
		{
			List list = new ArrayList<>();
			for (Class clz : types)
			{
				String name = clz.getName();
				if (name.contains("$")) name = name.substring(name.lastIndexOf("$") + 1);
				if (name.contains(".")) name = name.substring(name.lastIndexOf(".") + 1);
				list.add(name);
			}
			return list;
		}
	}

	enum ActivityLabel
	{
		ACTIVATED,
		INHIBITED;

		public String getName()
		{
			return toString().toLowerCase();
		}

		public static ActivityLabel getLabel(String s)
		{
			s = s.trim().toLowerCase();
			if (s.startsWith("i") || s.startsWith("-")) return INHIBITED;
			return ACTIVATED;
		}

		public static Map getValuesAsJson()
		{
			List list = new ArrayList<>();
			for (ActivityLabel label : values())
			{
				list.add(label.getName());
			}

			Map map = new LinkedHashMap<>();
			map.put("name", "ActivityLabel");
			map.put("values", list);
			return map;
		}
	}
}
