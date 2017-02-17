package org.panda.causalpath.network;

import com.github.jsonldjava.utils.JsonUtils;
import org.panda.causalpath.analyzer.NSCForNonCorr;
import org.panda.causalpath.analyzer.NetworkSignificanceCalculator;
import org.panda.causalpath.data.*;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * This class prepares the causality graph in two files, to display in ChiBE.
 */
public class GraphWriter
{
	/**
	 * Border color of activating phosphorylations or mutations.
	 */
	private Color activatingBorderColor;
	/**
	 * Border color of inhibiting phosphorylations or mutations.
	 */
	private Color inhibitingBorderColor;

	/**
	 * This color is used when the downstream of a gene suggest it is activated and inhibited at the same time.
	 */
	private Color doubleSignificanceBorderColor;

	/**
	 * Border color of nodes and phosphosites when their activating/inactivating status is not known, or their
	 * activated/inactivated status is not significant.
	 */
	private Color defaultBorderColor;

	/**
	 * Parameter to use gene background to display the total protein change. If false, then it is displayed on the gene
	 * node as a separate feature.
	 */
	private boolean useGeneBGForTotalProtein;

	/**
	 * The most intense color to show downregulation.
	 */
	private Color maxDownColor;

	/**
	 * The most intense color to show upregulation.
	 */
	private Color maxUpColor;

	/**
	 * The value where the most intense colors are reached. If a value is more extreme than this value, it will be
	 * shown with the most intense color, and would be considered "saturated".
	 */
	private double colorSaturationValue;

	/**
	 * An object that can produce a color for a given up/downregulation value.
	 */
	private ValToColor vtc;

	/**
	 * The set of relations with the associated data to draw.
	 */
	private Set<RelationAndSelectedData> relations;

	private NetworkSignificanceCalculator nsc;

	/**
	 * Constructor with the relations. Those relations are the result of the causality search.
	 */
	public GraphWriter(Set<RelationAndSelectedData> relations)
	{
		this(relations, null);
	}

	/**
	 * Constructor with the relations and significance calculation results.
	 * Those relations are the result of the causality search.
	 */
	public GraphWriter(Set<RelationAndSelectedData> relations, NetworkSignificanceCalculator nsc)
	{
		this.relations = relations;
		activatingBorderColor = new Color(0, 180, 20);
		inhibitingBorderColor = new Color(180, 0, 20);
		doubleSignificanceBorderColor = new Color(150, 150, 0);
		defaultBorderColor = new Color(50, 50, 50);

		maxDownColor = new Color(40, 80, 255);
		maxUpColor = new Color(255, 80, 40);
		colorSaturationValue = 1;

		initColorMapper();

		useGeneBGForTotalProtein = false;

		this.nsc = nsc;
	}

	/**
	 * Saturation color for downregulation.
	 */
	public void setMaxDownColor(Color maxDownColor)
	{
		this.maxDownColor = maxDownColor;
		initColorMapper();
	}

	/**
	 * Saturation color for upregulation.
	 */
	public void setMaxUpColor(Color maxUpColor)
	{
		this.maxUpColor = maxUpColor;
		initColorMapper();
	}

	/**
	 * The value where color saturation occurs. The parameter has to be a positive value. Negative saturation will be
	 * symmetrical to positive saturation.
	 */
	public void setColorSaturationValue(double colorSaturationValue)
	{
		this.colorSaturationValue = Math.abs(colorSaturationValue);
		initColorMapper();
	}

	/**
	 * Initializes the color mapping object.
	 */
	private void initColorMapper()
	{
		vtc = new ValToColor(new double[]{-colorSaturationValue, 0, colorSaturationValue},
			new Color[]{maxDownColor, Color.WHITE, maxUpColor});
	}

	public void setUseGeneBGForTotalProtein(boolean useGeneBGForTotalProtein)
	{
		this.useGeneBGForTotalProtein = useGeneBGForTotalProtein;
	}

	public void setActivatingBorderColor(Color activatingBorderColor)
	{
		this.activatingBorderColor = activatingBorderColor;
	}

	public void setInhibitingBorderColor(Color inhibitingBorderColor)
	{
		this.inhibitingBorderColor = inhibitingBorderColor;
	}

	public void setExpColorSchema(ValToColor vtc)
	{
		this.vtc = vtc;
	}

	/**
	 * Produces a causality graph where each node corresponds to a gene. In this graph, data may be displayed more than
	 * once if they map to more than one gene. The output .sif and .format files can be visualized using ChiBE. From
	 * ChiBE, do SIF -> Load SIF File ... and select the output .sif file.
	 */
	public void writeSIFGeneCentric(String filename) throws IOException
	{
		if (!filename.endsWith(".sif")) filename += ".sif";

		// write relations
		BufferedWriter writer1 = new BufferedWriter(new FileWriter(filename));
		relations.forEach(r -> FileUtil.writeln(r.relation.toString(), writer1));
		writer1.close();

		Set<String> totalProtUsedUp = new HashSet<>();

		filename = filename.substring(0, filename.lastIndexOf(".")) + ".format";
		BufferedWriter writer2 = new BufferedWriter(new FileWriter(filename));
		writer2.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer2.write("node\tall-nodes\tbordercolor\t" + inString(defaultBorderColor) + "\n");

		Set<ExperimentData> dataInGraph = getExperimentDataToDraw();

		dataInGraph.forEach(data ->
		{
			String colS = "255 255 255";
			String bor = inString(defaultBorderColor);
			String let = "x";

			if (data.hasChangeDetector())
			{
				if (data.getChangeSign() != 0) colS = vtc.getColorInString(data.getChangeValue());
			}

			if (data instanceof PhosphoProteinData)
			{
				PhosphoProteinData pd = (PhosphoProteinData) data;
				if (pd.getEffect() > 0) bor = inString(activatingBorderColor);
				else if (pd.getEffect() < 0) bor = inString(inhibitingBorderColor);

				let = "p";
			}
			else if (data instanceof ProteinData)
			{
				let = "t";
			}
			else if (data instanceof MutationData)
			{
				let = "m";
				if (data.getEffect() == 1)
					bor = inString(activatingBorderColor);
				else bor = inString(inhibitingBorderColor);
			}
			else if (data instanceof CNAData)
			{
				let = "c";
			}
			else if (data instanceof ExpressionData)
			{
				let = "e";
			}
			else if (data instanceof ActivityData)
			{
				let = "a";
				bor = inString(activatingBorderColor);
			}

			String tooltip = data.id;
			if (data.hasChangeDetector()) tooltip += ", " + data.getChangeValue();

			for (String gene : data.getGeneSymbols())
			{
				if (useGeneBGForTotalProtein && let.equals("t") && !totalProtUsedUp.contains(gene))
				{
					FileUtil.writeln("node\t" + gene + "\tcolor\t" + colS, writer2);
					FileUtil.writeln("node\t" + gene + "\ttooltip\t" + tooltip, writer2);
					if (nsc != null)
					{
						if (nsc.isDownstreamSignificant(gene))
						{
							FileUtil.writeln("node\t" + gene + "\tborderwidth\t2", writer2);
						}
						boolean act = false;
						boolean inh = false;

						if (nsc instanceof NSCForNonCorr)
						{
							act = ((NSCForNonCorr) nsc).isActivatingTargetsSignificant(gene);
							inh = ((NSCForNonCorr) nsc).isInhibitoryTargetsSignificant(gene);
						}

						if (act && !inh) FileUtil.writeln("node\t" + gene + "\tbordercolor\t" + inString(activatingBorderColor), writer2);
						else if (!act && inh) FileUtil.writeln("node\t" + gene + "\tbordercolor\t" + inString(inhibitingBorderColor), writer2);
						else if (act /** && inh **/) FileUtil.writeln("node\t" + gene + "\tbordercolor\t" + inString(doubleSignificanceBorderColor), writer2);
						//else FileUtil.writeln("node\t" + gene + "\tbordercolor\t" + inString(defaultBorderColor), writer2);
					}
					totalProtUsedUp.add(gene);
				}
				else
				{
					FileUtil.writeln("node\t" + gene + "\trppasite\t" + tooltip + "|" + let + "|" + colS + "|" + bor,
						writer2);
				}
			}
		});

		writer2.close();
	}

	/**
	 * Converts color to a string.
	 */
	private String inString(Color c)
	{
		return c.getRed() + " " + c.getGreen() + " " + c.getBlue();
	}

	/**
	 * Converts color to a JSON string.
	 */
	private String inJSONString(Color c)
	{
		return "rgb(" + c.getRed() + "," + c.getGreen() + "," + c.getBlue() + ")";
	}

	/**
	 * Generates a causality graph where each node is a measurement. In this graph, pathway relations can be displayed
	 * more than once if the same relation can explain more than one data pairs.
	 */
	public void writeSIFDataCentric(String filename) throws IOException
	{
		if (!filename.endsWith(".sif")) filename += ".sif";

		// write relations
		BufferedWriter writer1 = new BufferedWriter(new FileWriter(filename));
		relations.forEach(r -> FileUtil.writeln(r.source.id + "\t" + r.relation.type.name + "\t" + r.target.id, writer1));
		writer1.close();

		filename = filename.substring(0, filename.lastIndexOf(".")) + ".format";
		BufferedWriter writer2 = new BufferedWriter(new FileWriter(filename));
		writer2.write("node\tall-nodes\tcolor\t255 255 255\n");
		relations.stream().map(r -> new ExperimentData[]{r.source, r.target}).flatMap(Arrays::stream).distinct()
			.filter(data -> data.hasChangeDetector() && data.getChangeSign() != 0).forEach(data ->
				FileUtil.writeln("node\t" + data.id + "\tcolor\t" + vtc.getColorInString(data.getChangeValue()),
					writer2));

		writer2.close();
	}

	/**
	 * Writes the graph in JSON format, which can be viewed using the web-based proteomics analysis tool.
	 */
	public void writeJSON(String filename) throws IOException
	{
		if (!filename.endsWith(".json")) filename += ".json";

		Map<String, Object> map = new HashMap<>();
		List<Map> nodes = new ArrayList<>();
		List<Map> edges = new ArrayList<>();
		map.put("nodes", nodes);
		map.put("edges", edges);

		Map<String, Map> geneMap = new HashMap<>();
		Set<String> relMem = new HashSet<>();

		relations.forEach(rel ->
		{
			String key = rel.relation.source + "\t" + rel.relation.type.name + "\t" + rel.relation.target;
			if (relMem.contains(key)) return;
			else relMem.add(key);

			Map<String, Object> edge = new HashMap<>();
			edges.add(edge);
			Map<String, Object> dMap = new HashMap<>();
			edge.put("data", dMap);
			dMap.put("source", rel.relation.source);
			dMap.put("target", rel.relation.target);
			dMap.put("edgeType", rel.relation.type.name);
			dMap.put("tooltipText", CollectionUtil.merge(rel.relation.getTargetWithSites(0), ", "));
		});

		Set<String> totalProtUsedUp = new HashSet<>();

		Set<ExperimentData> dataInGraph = getExperimentDataToDraw();

		dataInGraph.forEach(data ->
		{
			String colS = "255 255 255";
			String bor = inJSONString(defaultBorderColor);
			String let = "x";

			if (data.hasChangeDetector())
			{
				if (data.getChangeSign() != 0) colS = vtc.getColorInJSONString(data.getChangeValue());
			}

			if (data instanceof PhosphoProteinData)
			{
				PhosphoProteinData pd = (PhosphoProteinData) data;
				if (pd.getEffect() > 0) bor = inJSONString(activatingBorderColor);
				else if (pd.getEffect() < 0) bor = inJSONString(inhibitingBorderColor);

				let = "p";
			}
			else if (data instanceof ProteinData)
			{
				let = "t";
			}
			else if (data instanceof MutationData)
			{
				let = "m";
				if (data.getEffect() == 1)
					bor = inJSONString(activatingBorderColor);
				else bor = inJSONString(inhibitingBorderColor);
			}
			else if (data instanceof CNAData)
			{
				let = "c";
			}
			else if (data instanceof ExpressionData)
			{
				let = "e";
			}
			else if (data instanceof ActivityData)
			{
				let = "a";
				bor = inJSONString(activatingBorderColor);
			}

			String tooltip = data.id;
			if (data.hasChangeDetector()) tooltip += ", " + data.getChangeValue();

			for (String sym : data.getGeneSymbols())
			{
				if (!geneMap.containsKey(sym))
				{
					HashMap<String, Object> node = new HashMap<>();
					geneMap.put(sym, node);
					nodes.add(node);
					Map<String, Object> d = new HashMap<>();
					node.put("data", d);
					d.put("id", sym);
					d.put("text", sym);
					List<Map> sites = new ArrayList<>();
					d.put("sites", sites);
				}

				Map node = geneMap.get(sym);

				if (useGeneBGForTotalProtein && let.equals("t") && !totalProtUsedUp.contains(sym))
				{
					if (!node.containsKey("css")) node.put("css", new HashMap<>());
					((Map) node.get("css")).put("backgroundColor", colS);
					((Map) node.get("data")).put("tooltipText", tooltip);
					totalProtUsedUp.add(sym);
				}
				else
				{
					Map site = new HashMap();
					((List) ((Map) node.get("data")).get("sites")).add(site);

					site.put("siteText", let);
					site.put("siteInfo", tooltip);
					site.put("siteBackgroundColor", colS);
					site.put("siteBorderColor", bor);
				}
			}
		});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		JsonUtils.writePrettyPrint(writer, map);
		writer.close();
	}

	private Set<ExperimentData> getExperimentDataToDraw()
	{
		Set<ExperimentData> dataInGraph = relations.stream().map(r -> new ExperimentData[]{r.source, r.target})
			.flatMap(Arrays::stream).collect(Collectors.toSet());

		Set<String> symInGraph = relations.stream().map(r -> new String[]{r.relation.source, r.relation.target})
			.flatMap(Arrays::stream).collect(Collectors.toSet());

		Set<ExperimentData> otherData = Stream.concat(
			relations.stream().map(r -> r.relation).map(r -> r.sourceData).flatMap(Collection::stream),
			relations.stream().map(r -> r.relation).map(r -> r.targetData).flatMap(Collection::stream)).distinct()
			.filter(d -> (d instanceof ExpressionData || d instanceof CNAData) && !dataInGraph.contains(d))
			.filter(d -> symInGraph.contains(d.getGeneSymbols().iterator().next()))
			.filter(d -> d.hasChangeDetector() && d.getChangeSign() != 0)
			.collect(Collectors.toSet());

		dataInGraph.addAll(otherData);
		return dataInGraph;
	}
}
