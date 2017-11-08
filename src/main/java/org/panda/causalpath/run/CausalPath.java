package org.panda.causalpath.run;

import org.omg.CORBA.IDLType;
import org.panda.causalpath.analyzer.*;
import org.panda.causalpath.data.*;
import org.panda.causalpath.network.GraphFilter;
import org.panda.causalpath.network.GraphWriter;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.resource.NetworkLoader;
import org.panda.causalpath.resource.ProteomicsFileReader;
import org.panda.causalpath.resource.ProteomicsLoader;
import org.panda.causalpath.resource.TCGALoader;
import org.panda.resource.ResourceDirectory;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.resource.tcga.ProteomicsFileRow;
import org.panda.utility.FileUtil;
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
	private Set<String> valueColumn;

	/**
	 * Names of the columns that contain control values.
	 */
	private Set<String> controlValueColumn;

	/**
	 * Names of the columns that contain test values.
	 */
	private Set<String> testValueColumn;

	/**
	 * Whether or not to do a log transformation on the values read, before every other kinds of processing. Log base
	 * is 2.
	 */
	private boolean doLogTransfrorm;

	/**
	 * A threshold to determine significant changes whenever it applies. This can be a value threshold, or a p-value
	 * threshold depending on the value transformation type.
	 */
	private double thresholdForDataSignificance;

	/**
	 * An FDR threshold for the data that provides us p-values.
	 */
	private Map<DataType, Double> fdrThresholdForDataSignificance;

	/**
	 * Whether we should pool proteomics and phosphoproteomics during FDR adjustment. Pooling makes sense for RPPA data
	 * because they are coming from single experiment.
	 */
	private boolean poolProteomicsForFDRAdjustment = false;

	/**
	 * A threshold to determine if the network size or a downstream of a protein is significantly crowded.
	 */
	private double fdrThresholdForNetworkSignificance;

	/**
	 * A correlation threshold when the value transformation is correlation. This is optional, but if not used, then
	 * using the <code>thresholdForSignificance</code> is mandatory.
	 */
	private double pvalThresholdForCorrelation;

	/**
	 * A correlation FDR threshold when the value transformation is correlation. When this is used,
	 */
	private double fdrThresholdForCorrelation;

	/**
	 * Proximity threshold to infer effects of unknown sites by the known effect of neighbor sites.
	 */
	private int siteEffectProximityThreshold;

	/**
	 * A method for interpretation of values in the data file.
	 */
	private ValueTransformation transformation;

	/**
	 * A default value for the missing values in the data file.
	 */
	private double defaultMissingValue;

	/**
	 * Parameter for calculating network significances. If this is true, then <code>permutationCount</code> should also
	 * be set.
	 */
	private boolean calculateNetworkSignificance;

	/**
	 * Number of iterations of permutations to use during network significance calculations.
	 */
	private int permutationCount;

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

	private boolean showUnexplainedProteomicData = false;

	private boolean useNetworkSignificanceForCausalReasoning = false;

	/**
	 * The search engine.
	 */
	CausalitySearcher cs;

	public static void main(String[] args) throws IOException
	{
		if (args.length < 1)
		{
			System.err.println("Please specify the data directory that contains the \"" + PARAMETER_FILENAME + "\" " +
				"file, as the first and only program parameter. All other parameters must go into the parameters " +
				"file.");
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
	public void run() throws IOException
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

		// Fill-in missing effects
		SiteEffectCollective sec = new SiteEffectCollective();
		sec.fillInMissingEffect(rows, siteEffectProximityThreshold);

		ProteomicsLoader loader = new ProteomicsLoader(rows);

		// Associate change detectors
//		loader.associateChangeDetector(getOneDataChangeDetector(ctrl, test), data -> data instanceof ProteinData);
		loader.associateChangeDetector(getOneDataChangeDetector(ctrl, test), data -> data instanceof PhosphoProteinData);
		loader.associateChangeDetector(getOneDataChangeDetector(ctrl, test), data -> data instanceof ProteinData &&
			!(data instanceof PhosphoProteinData));

		// Revisit this line after adding activity file support - todo
		loader.associateChangeDetector(new ThresholdDetector(0.1, ThresholdDetector.AveragingMethod.ARITHMETIC_MEAN), data -> data instanceof ActivityData);

		// Load signed relations
		Set<Relation> relations = networkSelection == null ? NetworkLoader.load() :
			NetworkLoader.load(NetworkLoader.ResourceType.getSelectedResources(networkSelection));

		loader.decorateRelations(relations);

		boolean useCorrelation = transformation == ValueTransformation.CORRELATION;
		CorrelationDetector corrDet = null;
		if (useCorrelation)
		{
			corrDet = new CorrelationDetector(pvalThresholdForCorrelation, thresholdForDataSignificance);
			if (minimumSampleSize != null) corrDet.setMinimumSampleSize(minimumSampleSize);
			if (correlationUpperThreshold != null) corrDet.setCorrelationUpperThreshold(correlationUpperThreshold);

			for (Relation rel : relations)
			{
				rel.setChDet(corrDet);
			}
		}

		// Load other TCGA profiles if available
		if (tcgaDirectory != null)
		{
			TCGALoader tcga = new TCGALoader(directory + File.separator + tcgaDirectory);
			tcga.setSamples(vals.toArray(new String[vals.size()]));

			if (mutationEffectFilename != null)
			{
				tcga.loadMutationEffectMap(directory + File.separator + mutationEffectFilename);
			}

			tcga.decorateRelations(relations);
			OneDataChangeDetector chDet = getOneDataChangeDetector(ctrl, test);

			if (chDet != null)
			{
				tcga.associateChangeDetector(chDet.makeACopy(), data -> data instanceof CNAData);
				tcga.associateChangeDetector(chDet.makeACopy(), data -> data instanceof RNAData);
				tcga.associateChangeDetector(chDet.makeACopy(), data -> data instanceof MutationData);
			}
		}

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
//		RobustnessAnalysis ra = new RobustnessAnalysis(cs, relations, fdrThresholdForDataSignificance, 0.05);
//		ra.run(100, directory + "/robustness.txt");
//		System.exit(0);
		//---END OF DEBUG

		// Search causal or conflicting relations

		boolean controlFDR = fdrThresholdForDataSignificance != null || fdrThresholdForCorrelation > 0;

		if (!controlFDR) cs.setCollectDataWithMissingEffect(true);
		else cs.setCollectDataUsedForInference(true);

		Set<Relation> causal = cs.run(relations);

		if (!controlFDR) cs.setCollectDataWithMissingEffect(false);

		if (controlFDR)
		{
			Set<ExperimentData> datas = new HashSet<>(cs.getDataUsedForInference());
			Set<Set<ExperimentData>> pairs = new HashSet<>(cs.getPairsUsedForInference());
			cs.setCausal(false);
			cs.run(relations);
			datas.addAll(cs.getDataUsedForInference());
			pairs.addAll(cs.getPairsUsedForInference());
			cs.setCausal(true);
			if (useCorrelation)
			{
				FDRAdjusterForCorrelation fad = new FDRAdjusterForCorrelation(pairs, corrDet);
				fad.adjustPValueThresholdsForFDR(fdrThresholdForCorrelation);
			}
			else
			{
				FDRAdjuster fad = new FDRAdjuster(poolProteomicsForFDRAdjustment);
				fad.adjustPValueThresholdsForRelations(relations, fdrThresholdForDataSignificance);
				fad.adjustPValueThresholdsOfDatas(datas, fdrThresholdForDataSignificance);
			}
			cs.setCollectDataWithMissingEffect(true);
			causal = cs.run(relations);

			cs.setCollectDataWithMissingEffect(false);
			cs.setCollectDataUsedForInference(false);
		}

		// Significance calculation
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

			nsc.run(permutationCount);
			nsc.setFDRThreshold(fdrThresholdForNetworkSignificance);
			nsc.writeResults(directory + File.separator + SIGNIFICANCE_FILENAME);
			System.out.println("Graph size pval = " + nsc.getOverallGraphSizePval());
		}

		// Add network significance as data if opted for

		if (calculateNetworkSignificance && !useCorrelation && useNetworkSignificanceForCausalReasoning)
		{
			Map<String, Integer> geneMap = ((NSCForNonCorr) nsc).getSignificantGenes();
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

				// Run the inference again with new activity data

				cs.setCollectDataWithMissingEffect(true);
				cs.setCollectDataUsedForInference(true);
				causal = cs.run(relations);
				cs.setCollectDataWithMissingEffect(false);
				cs.setCollectDataUsedForInference(false);
			}
		}

		int causativeSize = causal.size();
		System.out.println("Causative relations = " + causativeSize);

		GraphWriter writer = new GraphWriter(causal, nsc);
		writer.setUseGeneBGForTotalProtein(true);
		writer.setColorSaturationValue(colorSaturationValue);

		if (!useCorrelation && showUnexplainedProteomicData)
		{
			Set<GeneWithData> set = relations.stream().map(r -> r.sourceData).filter(GeneWithData::hasChangedProteinData).collect(Collectors.toSet());
			relations.stream().map(r -> r.targetData).filter(GeneWithData::hasChangedProteinData).forEach(set::add);
			writer.setOtherGenesToShow(set);
		}
		if (useCorrelation)
		{
			writer.setExperimentDataToDraw(cs.getPairsUsedForInference().stream().flatMap(Collection::stream)
				.collect(Collectors.toSet()));
		}

		// Generate output
		writer.writeSIFGeneCentric(directory + File.separator + CAUSATIVE_RESULT_FILE_PREFIX);
		writer.writeJSON(directory + File.separator + CAUSATIVE_RESULT_FILE_PREFIX);
		if (generateDataCentricGraph)
		{
			writer.writeSIFDataCentric(directory + File.separator + CAUSATIVE_RESULT_FILE_DATA_CENTRIC_PREFIX);
		}

		// Note the sites with unknown effect whose determination will improve the results
		writeSitesToCurate(cs.getDataNeedsAnnotation());

		// Do the same for conflicting relations

		cs.setCausal(false);
		cs.setCollectDataUsedForInference(true);
		Set<Relation> conflicting =  cs.run(relations);
		cs.setCollectDataUsedForInference(false);
		int conflictSize = conflicting.size();
		System.out.println("Conflicting relations = " + conflictSize);
		writer = new GraphWriter(conflicting, null);
		writer.setUseGeneBGForTotalProtein(true);
		writer.setColorSaturationValue(colorSaturationValue);
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
		double accuracy = pred.run(causal, conflicting, cs);
		System.out.println("accuracy = " + accuracy);
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

	private OneDataChangeDetector getOneDataChangeDetector(boolean[] ctrl, boolean[] test)
	{
		OneDataChangeDetector detector = null;

		if (transformation == ValueTransformation.ARITHMETIC_MEAN ||
			transformation == ValueTransformation.GEOMETRIC_MEAN ||
			transformation == ValueTransformation.MAX)
		{
			detector = new ThresholdDetector(thresholdForDataSignificance);
			((ThresholdDetector) detector).setAveragingMethod(transformation == ValueTransformation.ARITHMETIC_MEAN ?
				ThresholdDetector.AveragingMethod.ARITHMETIC_MEAN : transformation == ValueTransformation.GEOMETRIC_MEAN ?
				ThresholdDetector.AveragingMethod.FOLD_CHANGE_MEAN : ThresholdDetector.AveragingMethod.MAX);
		}
		else if (transformation == ValueTransformation.DIFFERENCE_OF_MEANS)
		{
			detector = new DifferenceDetector(thresholdForDataSignificance, ctrl, test);
		}
		else if (transformation == ValueTransformation.FOLD_CHANGE_OF_MEAN)
		{
			detector = new FoldChangeDetector(thresholdForDataSignificance, ctrl, test);
		}
		else if (transformation == ValueTransformation.SIGNIFICANT_CHANGE_OF_MEAN)
		{
			// if there is no fdr control, then use the given data significance threshold. Otherwise use 1, which means
			// make everything significant. And FDR correction will be applied later.
			detector = new SignificanceDetector(fdrThresholdForDataSignificance == null ? thresholdForDataSignificance : 1,
				ctrl, test);
		}
		return detector;
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
			"detection of a single change. There should only be one group of values (marked with " +
			Parameter.VALUE_COLUMN.getText() + "), the values have to be distributed around zero, and a threshold " +
			"value should be provided for significance detection, using the " +
			Parameter.THRESHOLD_FOR_DATA_SIGNIFICANCE.getText() + ".", false),

		GEOMETRIC_MEAN("geometric-mean", "The geometric mean value of the given values is used for significance " +
			"detection of a single change. This is the only case when the geometric mean is used for averaging a " +
			"group of samples, and it is appropriate if the individual values are formed of some kind of ratios. " +
			"There should only be one group of values (marked with " + Parameter.VALUE_COLUMN.getText() + "), the " +
			"values have to be distributed around zero, and a threshold value should be provided for significance " +
			"detection, using the " + Parameter.THRESHOLD_FOR_DATA_SIGNIFICANCE.getText() + ".", false),

		MAX("max", "The value with maximum absolute is used for the analysis. There should only be one group of " +
			"values (marked with " + Parameter.VALUE_COLUMN.getText() + "), the values have to be distributed around " +
			"zero, and a threshold value should be provided for significance detection, using the " +
			Parameter.THRESHOLD_FOR_DATA_SIGNIFICANCE.getText() + ".", false),

		DIFFERENCE_OF_MEANS("difference-of-means", "There should be control and test values, whose difference would " +
			"be used for significance detection. The threshold for significance (" +
			Parameter.THRESHOLD_FOR_DATA_SIGNIFICANCE.getText() + ") should also be provided.", true),

		FOLD_CHANGE_OF_MEAN("fold-change-of-mean", "There should be control and test values, whose ratio will be " +
			"converted to fold change and thresholded. The fold change value will be in the range (-inf, -1] + [1, " +
			"inf). If the data file already contains a fold-change value, then please use the " + GEOMETRIC_MEAN.name +
			" as value transformation. The threshold for significance (" +
			Parameter.THRESHOLD_FOR_DATA_SIGNIFICANCE.getText() + ") should also be provided.", true),

		SIGNIFICANT_CHANGE_OF_MEAN("significant-change-of-mean", "There should be sufficient amount of control and " +
			"test values to detect the significance of change with a t-test. Technically there should be more than 3" +
			" controls and 3 tests, practically, they should be much more to provide statistical power. The " +
			Parameter.THRESHOLD_FOR_DATA_SIGNIFICANCE.getText() + " should be used for a p-value threshold, or " +
			"alternatively, " + Parameter.FDR_THRESHOLD_FOR_DATA_SIGNIFICANCE.getText() + " should be used for " +
			"controlling significance at the false discovery rate level.", true),

		CORRELATION("correlation", "There should be one group of values (marked with " +
			Parameter.VALUE_COLUMN.getText() + "). There must be at least 3 value columns technically, but many more " +
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
	}

	enum Parameter
	{
		PROTEOMICS_PLATFORM_FILE((value, cp) -> cp.proteomicsPlatformFile = value,
			"Name of the proteomics platform file. Each row should belong to either a gene's total protein " +
				"measurement, or a site specific measurement. This file should contain ID, gene symbols, modification" +
				" sites, and known site effects. Platform file and values file can be the same file."),
		PROTEOMICS_VALUES_FILE((value, cp) -> cp.proteomicsValuesFile = value,
			"Name of the proteomics values file. It should have at least one ID column and one or more columns for " +
				"experiment values. Platform file and values file can be the same file."),
		ID_COLUMN((value, cp) -> cp.IDColumn = value,
			"The name of the ID column in platform and values files."),
		SYMBOLS_COLUMN((value, cp) -> cp.symbolsColumn = value,
			"The name of the symbols column in platform file."),
		SITES_COLUMN((value, cp) -> cp.sitesColumn = value,
			"The name of the sites column."),
		EFFECT_COLUMN((value, cp) -> cp.effectColumn = value,
			"The name of the effect column."),
		VALUE_COLUMN((value, cp) ->
		{
			if (cp.valueColumn == null) cp.valueColumn = new HashSet<>();
			cp.valueColumn.add(value);
		},
			"Name of a value column. This parameter should be used when there is only one group of experiments to " +
				"consider in the analysis."),
		CONTROL_VALUE_COLUMN((value, cp) ->
		{
			if (cp.controlValueColumn == null) cp.controlValueColumn = new HashSet<>();
			cp.controlValueColumn.add(value);
		},
			"Name of a control value column. This parameter should be used when there are control and test value " +
				"columns in the dataset."),
		TEST_VALUE_COLUMN((value, cp) ->
		{
			if (cp.testValueColumn == null) cp.testValueColumn = new HashSet<>();
			cp.testValueColumn.add(value);
		},
			"Name of a test value column. This parameter should be used when there are control and test value " +
				"columns in the dataset."),
		DO_LOG_TRANSFORM((value, cp) -> cp.doLogTransfrorm = Boolean.valueOf(value),
			"Whether the proteomic values should be log transformed for the analysis. Possible values are 'true' and " +
				"'false'."),
		THRESHOLD_FOR_DATA_SIGNIFICANCE((value, cp) -> cp.thresholdForDataSignificance = Double.valueOf(value),
			"A threshold value for selecting significant data. Use this parameter only when FDR controlling procedure" +
				"is already performed outside of CausalPath."),
		FDR_THRESHOLD_FOR_DATA_SIGNIFICANCE((value, cp) ->
		{
			if (cp.fdrThresholdForDataSignificance == null) cp.fdrThresholdForDataSignificance = new HashMap<>();
			String[] s = value.split("\\s+");
			cp.fdrThresholdForDataSignificance.put(DataType.get(s[1]), Double.valueOf(s[0]));
		},
			"False discovery rate threshold for data significance. This parameter can be use multiple times for each " +
				"different data type. The parameter value has to be in the form 'fdr-val data-type', such like " +
				"'0.1 phosphoprotein'. Possible data types are below" ),
		POOL_PROTEOMICS_FOR_FDR_ADJUSTMENT((value, cp) -> cp.poolProteomicsForFDRAdjustment = Boolean.valueOf(value),
			"Whether to consider proteomic and phosphoproteomic data as a single dataset during FDR adjustment. This" +
				"is typically the case with RPPA data, and typically not the case with mass spectrometry data. Can be" +
				" 'true' or 'false'. Default is false."),
		FDR_THRESHOLD_FOR_NETWORK_SIGNIFICANCE((value, cp) ->
			cp.fdrThresholdForNetworkSignificance = Double.valueOf(value),
			"The false discovery rate for network significance calculations for the downstream activity enrichment of" +
				" genes."),
		PVAL_THRESHOLD_FOR_CORRELATION((value, cp) -> cp.pvalThresholdForCorrelation = Double.valueOf(value),
			"A p-value threshold for correlation in a correlation-based causality. This parameter should only be used" +
				" when FDR control is performed outside of CausalPath."),
		FDR_THRESHOLD_FOR_CORRELATION((value, cp) -> cp.fdrThresholdForCorrelation = Double.valueOf(value),
			"False discovery rate threshold for the correlations in a correlation-based analysis."),
		VALUE_TRANSFORMATION((value, cp) -> cp.transformation = ValueTransformation.fetch(value),
			"This parameter determines how to use the values in the proteomics file. Possible values are listed " +
				"below.\n" + ValueTransformation.getUsageInfo()),
		CUSTOM_RESOURCE_DIRECTORY((value, cp) -> ResourceDirectory.set(value),
			"CausalPath downloads some data in the first run and stores in the resource directory. This directory is " +
				"'.panda' by default. If this needs to be customized, use this parameter."),
		SITE_EFFECT_PROXIMITY_THRESHOLD((value, cp) -> cp.siteEffectProximityThreshold = Integer.valueOf(value),
			"CausalPath has a database of phosphorylation site effects. When not set, this parameter is 0 by default," +
				" which means exact usage of site effects. But sometimes we may see a site with unknown effect is " +
				"modified, which is very close to another site with a known effect. This parameter let's us to assume" +
				" those changing sites with unknown effect has the same effect with the neighbor site with known " +
				"effect. Use responsibly."),
		SITE_MATCH_PROXIMITY_THRESHOLD((value, cp) -> cp.cs.setSiteProximityThreshold(Integer.valueOf(value)),
			"Phosphorylation relations many times know the target sites. When we observe a change in a site of the " +
				"target protein which is not targeted by the relation, but the site is very close to a known target " +
				"site of the relation, this parameter let's us to assume that the relation also applies to those " +
				"close-by sites."),
		DEFAULT_MISSING_VALUE((value, cp) -> cp.defaultMissingValue = Double.valueOf(value),
			"An option to specify a default value for the missing values in the proteomics file."),
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
			"Use this parameter to limit the results with a specific type of relation. Possible values are below.\n" +
				GraphFilter.RelationFilterType.getUsageInfo()),
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
			"Use this parameter to crop the result network to the neighborhood of certain gene. You should provide " +
				"gene symbols of these genes in a row separated by a semicolon, such like 'MTOR;RPS6KB1;RPS6'"),
		CALCULATE_NETWORK_SIGNIFICANCE((value, cp) -> cp.calculateNetworkSignificance = Boolean.valueOf(value),
			"Whether to calculate significances of the properties of the graph. When turned on, a p-value for network" +
				" size, and also downstream activity enrichment p-values for each gene on the graph are calculated."),
		PERMUTATIONS_FOR_SIGNIFICANCE((value, cp) -> cp.permutationCount = Integer.valueOf(value),
			""),
		TCGA_DIRECTORY((value, cp) -> cp.tcgaDirectory = value,
			"It is possible to add genomic data from TCGA to CausalPath analysis. This is only useful when the " +
				"proteomic data have the same sample IDs. Users can load TCGA data into a local directory from Broad " +
				"Firehose, and provide the directory here. The org.panda.resource.tcga.BroadDownloader in the project" +
				" https://github.com/PathwayAndDataAnalysis/resource is a utility that can do that."),
		MUTATION_EFFECT_FILE((value, cp) -> cp.mutationEffectFilename = value,
			"When we have mutations in the analysis, users can provide mutation effects using this parameter, " +
				"otherwise all mutations are assumed to be inactivating."),
		COLOR_SATURATION_VALUE((value, cp) -> cp.colorSaturationValue = Double.valueOf(value),
			"The saturation value to use on the result graphs."),
		DO_SITE_MATCHING((value, cp) -> cp.cs.setForceSiteMatching(Boolean.valueOf(value)),
			"Whether to force site matching in causality analysis. True by default."),
		SHOW_UNEXPLAINED_PROTEOMIC_DATA((value, cp) -> cp.showUnexplainedProteomicData = Boolean.valueOf(value),
			"CausalPath generates a result graph, but what about all other significant changes that could not make " +
				"into the network? CausalPath puts those genes as disconnected nodes in the graph when the analysis " +
				"is not correlation based. This is true by default but can be turned off by setting to false."),
		BUILT_IN_NETWORK_RESOURCE_SELECTION((value, cp) -> cp.networkSelection = value,
			"Determines which network resource to use during the analysis. Multiple network resource should be " +
				"mentioned together, separated with a space or comma. Possible values are below.\n" +
				NetworkLoader.ResourceType.getUsageInfo()),
		GENERATE_DATA_CENTRIC_GRAPH((value, cp) -> cp.generateDataCentricGraph = Boolean.valueOf(value),
			"An alternative to the gene-centric graph of CausalPath is a data-centric graph where nodes are not genes" +
				" but the data. This parameter forces to generate this type of result as well. False by default."),
		CORRELATION_UPPER_THRESHOLD((value, cp) -> cp.correlationUpperThreshold = Double.valueOf(value),
			"In some types of proteomic data, highest correlations come from errors. A way around is filtering with " +
				"an upper value."),
		MINIMUM_SAMPLE_SIZE((value, cp) -> cp.minimumSampleSize = Integer.valueOf(value),
			"When there are missing values in proteomic file, the comparisons can have different sample sizes for " +
				"controls and tests. This parameter sets the minimum sample size of the control and test sets."),
		GENE_ACTIVITY((value, cp) ->
		{
			if (cp.activityMap == null) cp.activityMap = new HashMap<>();
			String[] t = value.split(" ");
			cp.activityMap.put(t[0], t[1].equals("a") ? 1 : t[1].equals("i") ? -1 : Integer.valueOf(t[1]));
		},
			"Use this parameter to assign a specific activity or inactivity to a gene in the analysis. The value has " +
				"to start with a gene name and one letter code for activity or inactivity, such as 'BRAF a', or 'PTEN" +
				" i'."),
		TF_ACTIVITY_FILE((value, cp) -> cp.tfActivityFile = value,
			"CausalPath lets users to input results from an inference for transcriptional factor activities, such as" +
				" PRECEPTS, MARINa, or VIPER. For this, the results should be prepared in a file, first column " +
				"containing TF symbol and the second column whether 'activated' or 'inhibited'. The name of such file" +
				" should be provided here."),
		USE_STRONGEST_PROTEOMIC_DATA_PER_GENE((value, cp) ->
			cp.cs.setUseStrongestProteomicsDataForActivity(Boolean.valueOf(value)),
			"When a proteomic experiment outputs too many phosphorylation sites with lots of changes, many proteins " +
				"have evidences for both activation and inhibition. This produces networks hard to read. A complexity" +
				" management technique is to turn on this parameter to use only the strongest proteomic feature at " +
				"the upstream of relations. This is false by default."),
		USE_NETWORK_SIGNIFICANCE_FOR_CAUSAL_REASONING((value, cp) ->
			cp.useNetworkSignificanceForCausalReasoning = Boolean.valueOf(value),
			"After calculation of network significances in a non-correlation-based analysis, this option introduces" +
				" the detected active and inactive proteins as data to be used in the analysis. This applies only to " +
				"the proteins that already have a changed data on them, and have no previous activity data associated."),
		;

		ParameterReader reader;
		String info;

		Parameter(ParameterReader reader, String info)
		{
			this.reader = reader;
			this.info = info;
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

		public static String getUsageInfo()
		{
			StringBuilder sb = new StringBuilder();
			for (Parameter param : values())
			{
				sb.append(param.getText()).append(": ").append(param.getInfo()).append("\n");
			}
			return sb.toString();
		}
	}

	interface ParameterReader
	{
		void read(String value, CausalPath cp);
	}
}
