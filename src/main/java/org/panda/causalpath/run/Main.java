package org.panda.causalpath.run;

import org.panda.causalpath.analyzer.*;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.PhosphoProteinData;
import org.panda.causalpath.data.ProteinData;
import org.panda.causalpath.network.GraphFilter;
import org.panda.causalpath.network.GraphWriter;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationAndSelectedData;
import org.panda.causalpath.resource.NetworkLoader;
import org.panda.causalpath.resource.ProteomicsFileReader;
import org.panda.causalpath.resource.ProteomicsLoader;
import org.panda.causalpath.resource.TCGALoader;
import org.panda.resource.PhosphoSitePlus;
import org.panda.resource.ResourceDirectory;
import org.panda.resource.tcga.ProteomicsFileRow;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * This method is built to handle running the causality analysis using its jar file and pointing it to a related
 * resource directory.
 *
 * @author Ozgun Babur
 */
public class Main
{
	/**
	 * Name of the parameters file.
	 */
	public static final String PARAMETER_FILENAME = "parameters.txt";
	public static final String CAUSATIVE_RESULT_FILE_PREFIX = "causative";
	public static final String CONFLICTING_RESULT_FILE_PREFIX = "conflicting";
	public static final String UNKNOWN_SITE_EFFECT_FILENAME = "unknown-site-effects.txt";
	public static final String SIGNIFICANCE_FILENAME = "significance-pvals.txt";

	public static final String PROTEOMICS_PLATFORM_FILE_KEY = "proteomics-platform-file";
	public static final String PROTEOMICS_VALUES_FILE_KEY = "proteomics-values-file";
	public static final String ID_COLUMN_KEY = "id-column";
	public static final String SYMBOLS_COLUMN_KEY = "symbols-column";
	public static final String SITES_COLUMN_KEY = "sites-column";
	public static final String EFFECT_COLUMN_KEY = "effect-column";
	public static final String VALUE_COLUMN_KEY = "value-column";
	public static final String CONTROL_VALUE_COLUMN_KEY = "control-value-column";
	public static final String TEST_VALUE_COLUMN_KEY = "test-value-column";
	public static final String DO_LOG_TRANSFORM_KEY = "do-log-transform";
	public static final String THRESHOLD_FOR_SIGNIFICANCE_KEY = "threshold-for-significance";
	public static final String THRESHOLD_FOR_CORRELATION_KEY = "threshold-for-correlation";
	public static final String THRESHOLD_FOR_FDR_KEY = "threshold-for-FDR";
	public static final String VALUE_TRANSFORMATION_KEY = "value-transformation";
	public static final String CUSTOM_RESOURCE_DIRECTORY_KEY = "custom-resource-directory";
	public static final String SITE_EFFECT_PROXIMITY_THRESHOLD_KEY = "site-effect-proximity-threshold";
	public static final String SITE_MATCH_PROXIMITY_THRESHOLD_KEY = "site-match-proximity-threshold";
	public static final String DEFAULT_MISSING_VALUE_KEY = "default-missing-value";
	public static final String RELATION_FILTER_TYPE_KEY = "relation-filter-type";
	public static final String CALCULATE_NETWORK_SIGNIFICANCE_KEY = "calculate-network-significance";
	public static final String PERMUTATIONS_FOR_SIGNIFICANCE_KEY = "permutations-for-significance";
	public static final String TCGA_DIRECTORY_KEY = "tcga-directory";
	public static final String MUTATION_EFFECT_FILE_KEY = "mutation-effect-file";

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
	private double thresholdForSignificance;

	/**
	 * A correlation threshold when the value transformation is correlation. This is optional, but if not used, then
	 * using the <code>thresholdForSignificance</code> is mandatory.
	 */
	private double thresholdForCorrelation;

	/**
	 * A threshold for FDR cutoff, if preferred. IF this is provided, and applicable, then the significance threshold is
	 * ignored.
	 */
	private double thresholdForFDR;

	/**
	 * Proximity threshold to infer effects of unknown sites by the known effect of neighbor sites.
	 */
	private int siteEffectProximityThreshold;

	/**
	 * Proximity threshold to infer A phosphorylates B at the new site, just because we know that A phosphorylates B at
	 * a close-by site.
	 */
	private int siteMatchProximityThreshold;

	/**
	 * A method for interpretation of values in the data file.
	 */
	private ValueTransformation transformation;

	/**
	 * A default value for the missing values in the data file.
	 */
	private double defaultMissingValue;

	/**
	 * A filter for selecting a subgraph from the results.
	 */
	private GraphFilter graphFilter;

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

	private String mutationEffectFilename;

	public static void main(String[] args) throws IOException
	{
		if (args.length < 1)
		{
			System.err.println("Please specify the data directory that contains the \"" + PARAMETER_FILENAME + "\" " +
				"file, as the first and only program parameter. All other parameters must go into the parameters " +
				"file.");
			return;
		}

		new Main(args[0]).run();
	}

	public Main(String directory) throws IOException
	{
		this.directory = directory;
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
		switch (key)
		{
			case PROTEOMICS_PLATFORM_FILE_KEY: proteomicsPlatformFile = value; break;
			case PROTEOMICS_VALUES_FILE_KEY: proteomicsValuesFile = value; break;
			case ID_COLUMN_KEY: IDColumn = value; break;
			case SYMBOLS_COLUMN_KEY: symbolsColumn = value; break;
			case SITES_COLUMN_KEY: sitesColumn = value; break;
			case EFFECT_COLUMN_KEY: effectColumn = value; break;
			case VALUE_COLUMN_KEY: if (valueColumn == null) valueColumn = new HashSet<>();
				valueColumn.add(value); break;
			case CONTROL_VALUE_COLUMN_KEY: if (controlValueColumn == null) controlValueColumn = new HashSet<>();
				controlValueColumn.add(value); break;
			case TEST_VALUE_COLUMN_KEY: if (testValueColumn == null) testValueColumn = new HashSet<>();
				testValueColumn.add(value); break;
			case DO_LOG_TRANSFORM_KEY: doLogTransfrorm = Boolean.valueOf(value); break;
			case THRESHOLD_FOR_SIGNIFICANCE_KEY: thresholdForSignificance = Double.valueOf(value); break;
			case THRESHOLD_FOR_CORRELATION_KEY: thresholdForCorrelation = Double.valueOf(value); break;
			case THRESHOLD_FOR_FDR_KEY: thresholdForFDR = Double.valueOf(value); break;
			case VALUE_TRANSFORMATION_KEY: transformation = ValueTransformation.fetch(value); break;
			case CUSTOM_RESOURCE_DIRECTORY_KEY: ResourceDirectory.set(value); break;
			case SITE_EFFECT_PROXIMITY_THRESHOLD_KEY: siteEffectProximityThreshold = Integer.valueOf(value); break;
			case SITE_MATCH_PROXIMITY_THRESHOLD_KEY: siteMatchProximityThreshold = Integer.valueOf(value); break;
			case DEFAULT_MISSING_VALUE_KEY: defaultMissingValue = Double.valueOf(value); break;
			case RELATION_FILTER_TYPE_KEY: graphFilter = new GraphFilter(value); break;
			case CALCULATE_NETWORK_SIGNIFICANCE_KEY: calculateNetworkSignificance = Boolean.valueOf(value); break;
			case PERMUTATIONS_FOR_SIGNIFICANCE_KEY: permutationCount = Integer.valueOf(value); break;
			case TCGA_DIRECTORY_KEY: tcgaDirectory = value; break;
			case MUTATION_EFFECT_FILE_KEY: mutationEffectFilename = value; break;
			default: throw new RuntimeException("Unknown parameter: " + key);
		}
	}

	/**
	 * Executes the analysis;
	 */
	public void run() throws IOException
	{
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

		// Fill-in missing effect from PhosphoSitePlus
		PhosphoSitePlus.get().fillInMissingEffect(rows, siteEffectProximityThreshold);

		ProteomicsLoader loader = new ProteomicsLoader(rows);

		// Prepare a change-detector tailored to need

		OneDataChangeDetector detector = null;

		if (transformation == ValueTransformation.ARITHMETIC_MEAN ||
			transformation == ValueTransformation.GEOMETRIC_MEAN)
		{
			detector = new ThresholdDetector(thresholdForSignificance);
			((ThresholdDetector) detector).setGeometricMean(transformation == ValueTransformation.GEOMETRIC_MEAN);
		}
		else if (transformation == ValueTransformation.DIFFERENCE_OF_MEANS)
		{
			detector = new DifferenceDetector(thresholdForSignificance, ctrl, test);
		}
		else if (transformation == ValueTransformation.FOLD_CHANGE_OF_MEAN)
		{
			detector = new FoldChangeDetector(thresholdForSignificance, ctrl, test);
		}
		else if (transformation == ValueTransformation.SIGNIFICANT_CHANGE_OF_MEAN)
		{
			detector = new SignificanceDetector(thresholdForSignificance, ctrl, test);
		}

		// Associate change detectors
		loader.associateChangeDetector(detector, data -> data instanceof ProteinData);

		// Prepare relation-target compatibility checker
		RelationTargetCompatibilityChecker rtcc = new RelationTargetCompatibilityChecker();
		rtcc.setForceSiteMatching(true);
		rtcc.setSiteProximityThreshold(siteMatchProximityThreshold);

		// Load signed relations
		Set<Relation> relations = NetworkLoader.load();
		loader.decorateRelations(relations, rtcc);

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
			tcga.associateChangeDetector(detector, data -> true);
		}

		if (transformation == ValueTransformation.CORRELATION)
		{
			CorrelationDetector corrDet = new CorrelationDetector(thresholdForCorrelation, thresholdForSignificance);

			for (Relation rel : relations)
			{
				rel.setChDet(corrDet);
			}
		}

		// Prepare causality searcher
		CausalitySearcher cs = new CausalitySearcher(rtcc);

		// Search causal or conflicting relations
		Set<RelationAndSelectedData> relDat =  cs.run(relations);
		if (graphFilter != null) relDat = graphFilter.filter(relDat);
		int causativeSize = (int) relDat.stream().map(r -> r.relation).distinct().count();
		System.out.println("Causative relations = " + causativeSize);

		// Significance calculation
		NetworkSignificanceCalculator nsc = null;

		if (calculateNetworkSignificance && transformation != ValueTransformation.CORRELATION)
		{
			nsc = new NetworkSignificanceCalculator(relations, true, siteMatchProximityThreshold, true, graphFilter);
			nsc.run(permutationCount);
			nsc.setSignificanceThreshold(thresholdForSignificance);
			nsc.writeResults(directory + File.separator + SIGNIFICANCE_FILENAME);
			System.out.println("Graph size pval = " + nsc.getOverallGraphSizePval());
		}

		GraphWriter writer = new GraphWriter(relDat, nsc);
		writer.setUseGeneBGForTotalProtein(true);

		// Generate output
		writer.writeSIFGeneCentric(directory + File.separator + CAUSATIVE_RESULT_FILE_PREFIX);
		writer.writeJSON(directory + File.separator + CAUSATIVE_RESULT_FILE_PREFIX);

		// Note the sites with unknown effect whose determination will improve the results
		writeData(cs.findDataThatNeedsAnnotation(relations, relDat));

		// Do the same for conflicting relations

		cs.setCausal(false);
		relDat =  cs.run(relations);
		if (graphFilter != null) relDat = graphFilter.filter(relDat);
		int conflictSize = (int) relDat.stream().map(r -> r.relation).distinct().count();
		System.out.println("Conflicting relations = " + conflictSize);
		writer = new GraphWriter(relDat, null);
		writer.setUseGeneBGForTotalProtein(true);
		writer.writeSIFGeneCentric(directory + File.separator + CONFLICTING_RESULT_FILE_PREFIX);
		writer.writeJSON(directory + File.separator + CONFLICTING_RESULT_FILE_PREFIX);

		// Estimate false discovery rate.
		if (causativeSize > 0)
		{
			System.out.println("conflict / causative ratio = " + conflictSize / (double) causativeSize);
		}
	}

	private void writeData(Set<ExperimentData> datas) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(
			Paths.get(directory + File.separator + UNKNOWN_SITE_EFFECT_FILENAME));

		for (ExperimentData data : datas)
		{
			if (data instanceof PhosphoProteinData)
			{
				PhosphoProteinData pdata = (PhosphoProteinData) data;
				for (String s : pdata.getGenesWithSites())
				{
					writer.write(s + "\n");
				}
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
			"detection of a single change. There should only be one group of values (marked with " + VALUE_COLUMN_KEY  +
			"), the values have to be distributed around zero, and a threshold value should be provided for " +
			"significance detection, using the " + THRESHOLD_FOR_SIGNIFICANCE_KEY + ".", false),

		GEOMETRIC_MEAN("geometric-mean", "The geometric mean value of the given values is used for significance " +
			"detection of a single change. This is the only case when the geometric mean is used for averaging a " +
			"group of samples, and it is appropriate if the individual values are formed of some kind of ratios. " +
			"There should only be one group of values (marked with " + VALUE_COLUMN_KEY  + "), the values have to be " +
			"distributed around zero, and a threshold value should be provided for significance detection, using the " +
			THRESHOLD_FOR_SIGNIFICANCE_KEY + ".", false),

		DIFFERENCE_OF_MEANS("difference-of-means", "There should be control and test values, whose difference would " +
			"be used for significance detection. The threshold for significance (" + THRESHOLD_FOR_SIGNIFICANCE_KEY +
			") should also be provided.", true),

		FOLD_CHANGE_OF_MEAN("fold-change-of-mean", "There should be control and test values, whose ratio will be " +
			"converted to fold change and thresholded. The fold change value will be in the range (-inf, -1] + [1, " +
			"inf). If the data file already contains a fold-change value, then please use the " + GEOMETRIC_MEAN.name +
			" as value transformation. The threshold for significance (" + THRESHOLD_FOR_SIGNIFICANCE_KEY +
			") should also be provided.", true),

		SIGNIFICANT_CHANGE_OF_MEAN("significant-change-of-mean", "There should be sufficient amount of control and " +
			"test values to detect the significance of change with a t-test. Technically there should be more than 3" +
			" controls and 3 tests, practically, they should be much more to provide statistical power. The " +
			THRESHOLD_FOR_SIGNIFICANCE_KEY + " should be used for a p-value threshold, or alternatively, " + THRESHOLD_FOR_FDR_KEY +
			" should be used for controlling significance at the false discovery rate level.", true),

		CORRELATION("correlation", "There should be one group of values (marked with " + VALUE_COLUMN_KEY + "). There " +
			"must be at least 3 value columns technically, but many more than that practically to have some " +
			"statistical power for significant correlation. ", false);

		ValueTransformation(String name, String description, boolean twoGroupComparison)
		{
			this.name = name;
			this.description = description;
			this.twoGroupComparison = twoGroupComparison;
		}

		String name;
		String description;
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
	}
}
