package org.panda.causalpath.network;

import com.github.jsonldjava.utils.JsonUtils;
import org.panda.causalpath.analyzer.NSCForComparison;
import org.panda.causalpath.analyzer.NetworkSignificanceCalculator;
import org.panda.causalpath.data.*;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
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
	private Set<Relation> relations;

	private NetworkSignificanceCalculator nsc;

	private Set<ExperimentData> experimentDataToDraw;

	private Set<GeneWithData> otherGenesToShow;

	private boolean showInsignificantData = false;

	/**
	 * Constructor with the relations. Those relations are the result of the causality search.
	 */
	public GraphWriter(Set<Relation> relations)
	{
		this(relations, null);
	}

	/**
	 * Constructor with the relations and significance calculation results.
	 * Those relations are the result of the causality search.
	 */
	public GraphWriter(Set<Relation> relations, NetworkSignificanceCalculator nsc)
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

	public void setExperimentDataToDraw(Set<ExperimentData> experimentDataToDraw)
	{
		this.experimentDataToDraw = experimentDataToDraw;
	}

	public void setOtherGenesToShow(Set<GeneWithData> set)
	{
		this.otherGenesToShow = set;
	}

	public void setShowInsignificantData(boolean showInsignificantData)
	{
		this.showInsignificantData = showInsignificantData;
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
		relations.stream().distinct().forEach(r -> FileUtil.writeln(r.toString(), writer1));
		if (otherGenesToShow != null)
		{
			Set<String> genesInGraph = getGenesInGraph();
			otherGenesToShow.stream().filter(g -> !genesInGraph.contains(g.getId())).forEach(g ->
				FileUtil.writeln(g.getId(), writer1));
		}
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

			if (data.hasChangeDetector())
			{
				if (data.getChangeSign() != 0 || showInsignificantData)
				{
					colS = vtc.getColorInString(data.getChangeValue());
				}
				else return;
			}

			String bor = inString(defaultBorderColor);
			String let = "?";

			if (data instanceof SiteModProteinData)
			{
				SiteModProteinData pd = (SiteModProteinData) data;
				if (pd.getEffect() > 0) bor = inString(activatingBorderColor);
				else if (pd.getEffect() < 0) bor = inString(inhibitingBorderColor);

				let = pd.getModification().toString().substring(0, 1).toLowerCase();
			}
			else if (data instanceof ProteinData)
			{
				let = "t";
			}
			else if (data instanceof MetaboliteData)
			{
				let = "c";
			}
			else if (data instanceof MutationData)
			{
				let = "x";
				if (data.getEffect() == 1)
					bor = inString(activatingBorderColor);
				else bor = inString(inhibitingBorderColor);
			}
			else if (data instanceof CNAData)
			{
				let = "d";
			}
			else if (data instanceof RNAData)
			{
				let = "r";
			}
			else if (data instanceof ActivityData)
			{
				let = data.getChangeSign() > 0 ? "!" : "i";
				bor = inString(activatingBorderColor);
			}

			String siteID = data.id;
			String val = data.hasChangeDetector() ? data.getChangeValue() + "" : "";

			for (String gene : data.getGeneSymbols())
			{
				if (nsc != null)
				{
					if (nsc.isDownstreamSignificant(gene))
					{
						FileUtil.writeln("node\t" + gene + "\tborderwidth\t2", writer2);
					}
					boolean act = false;
					boolean inh = false;

					if (nsc instanceof NSCForComparison)
					{
						act = ((NSCForComparison) nsc).isActivatingTargetsSignificant(gene);
						inh = ((NSCForComparison) nsc).isInhibitoryTargetsSignificant(gene);
					}

					if (act && !inh) FileUtil.writeln("node\t" + gene + "\tbordercolor\t" + inString(activatingBorderColor), writer2);
					else if (!act && inh) FileUtil.writeln("node\t" + gene + "\tbordercolor\t" + inString(inhibitingBorderColor), writer2);
					else if (act /* && inh */) FileUtil.writeln("node\t" + gene + "\tbordercolor\t" + inString(doubleSignificanceBorderColor), writer2);
					//else FileUtil.writeln("node\t" + gene + "\tbordercolor\t" + inString(defaultBorderColor), writer2);
				}

				if (useGeneBGForTotalProtein && (let.equals("t") || let.equals("c")) && !totalProtUsedUp.contains(gene))
				{
//					if (let.equals("c"))
//					{
//						FileUtil.writeln("node\t" + siteID + "\tcolor\t" + colS, writer2);
//						FileUtil.writeln("node\t" + siteID + "\ttooltip\t" + gene + ", " + val, writer2);
//						totalProtUsedUp.add(gene);
//					}
//					else
//					{
						FileUtil.writeln("node\t" + gene + "\tcolor\t" + colS, writer2);
						FileUtil.writeln("node\t" + gene + "\ttooltip\t" + siteID + ", " + val, writer2);
						totalProtUsedUp.add(gene);
//					}
				}
				else
				{
					FileUtil.writeln("node\t" + gene + "\trppasite\t" + siteID.replaceAll("\\|", "-") + "|" + let + "|" + colS + "|" + bor +
						"|" + val, writer2);
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
	 *
	 * @param filename name of the output sif file
	 * @param unitsMap The map from the inferred relations to the data that helped the inference.
	 */
	public void writeSIFDataCentric(String filename, Map<Relation, Set<ExperimentData>> unitsMap) throws IOException
	{
		if (!filename.endsWith(".sif")) filename += ".sif";

		Set<ExperimentData> used = new HashSet<>();

		// write relations
		BufferedWriter writer1 = new BufferedWriter(new FileWriter(filename));
		unitsMap.keySet().forEach(r ->
		{
			Set<ExperimentData> datas = unitsMap.get(r);
			Set<ExperimentData> sources = r.sourceData.getData().stream().filter(datas::contains)
				.collect(Collectors.toSet());
			Set<ExperimentData> targets = r.targetData.getData().stream().filter(datas::contains)
				.collect(Collectors.toSet());

			for (ExperimentData source : sources)
			{
				for (ExperimentData target : targets)
				{
					FileUtil.writeln(source.getId() + "\t" + r.type.getName() + "\t" + target.getId(), writer1);
					used.add(source);
					used.add(target);
				}
			}
		});
		writer1.close();

		filename = filename.substring(0, filename.lastIndexOf(".")) + ".format";
		BufferedWriter writer2 = new BufferedWriter(new FileWriter(filename));
		writer2.write("node\tall-nodes\tcolor\t255 255 255\n");
		used.stream().forEach(data ->
		{
			if (data.hasChangeDetector())
			{
				FileUtil.writeln("node\t" + data.getId() + "\tcolor\t" + vtc.getColorInString(data.getChangeValue()),
					writer2);
			}
			if (data.getType() == DataType.PHOSPHOPROTEIN && data.getEffect() != 0)
			{
				FileUtil.writeln("node\t" + data.getId() + "\tbordercolor\t" +
					inString(data.getEffect() == -1 ? inhibitingBorderColor : activatingBorderColor), writer2);
			}
		});
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
			String src = rel.source.startsWith("CHEBI:") ? rel.sourceData.getData().iterator().next().getId() : rel.source;
			String tgt = rel.target.startsWith("CHEBI:") ? rel.targetData.getData().iterator().next().getId() : rel.target;
			String key = src + "\t" + rel.type.getName() + "\t" + tgt;
			if (relMem.contains(key)) return;
			else relMem.add(key);

			Map<String, Object> edge = new HashMap<>();
			edges.add(edge);
			Map<String, Object> dMap = new HashMap<>();
			edge.put("data", dMap);
			dMap.put("source", src);
			dMap.put("target", tgt);
			dMap.put("edgeType", rel.type.getName());
			dMap.put("tooltipText", CollectionUtil.merge(rel.getTargetWithSites(0), ", "));

			if (rel.getMediators() != null)
			{
				List<String> medList = Arrays.asList(rel.getMediators().split(";| "));
				if (!medList.isEmpty())
				{
					dMap.put("pcLinks", medList);
				}
			}
		});

		Set<String> totalProtUsedUp = new HashSet<>();

		Set<ExperimentData> dataInGraph = getExperimentDataToDraw();

		dataInGraph.forEach(data ->
		{
			String colS = "255 255 255";

			if (data.hasChangeDetector())
			{
				if (data.getChangeSign() != 0)
				{
					colS = vtc.getColorInJSONString(data.getChangeValue());
				}
				else return;
			}

			String bor = inJSONString(defaultBorderColor);
			String let = "?";

			if (data instanceof SiteModProteinData)
			{
				SiteModProteinData pd = (SiteModProteinData) data;
				if (pd.getEffect() > 0) bor = inJSONString(activatingBorderColor);
				else if (pd.getEffect() < 0) bor = inJSONString(inhibitingBorderColor);

				let = pd.getModification().toString().substring(0, 1).toLowerCase();
			}
			else if (data instanceof ProteinData)
			{
				let = "t";
			}
			else if (data instanceof MetaboliteData)
			{
				let = "c";
			}
			else if (data instanceof MutationData)
			{
				let = "x";
				if (data.getEffect() == 1)
					bor = inJSONString(activatingBorderColor);
				else bor = inJSONString(inhibitingBorderColor);
			}
			else if (data instanceof CNAData)
			{
				let = "d";
			}
			else if (data instanceof RNAData)
			{
				let = "r";
			}
			else if (data instanceof ActivityData)
			{
				let = data.getChangeSign() > 0 ? "!" : "i";
				bor = inJSONString(activatingBorderColor);
			}

			String siteID = data.id;
			String val = data.hasChangeDetector() ? data.getChangeValue() + "" : "";

			for (String sym : data.getGeneSymbols())
			{
				String nodeText = let.equals("c") ? siteID : sym;
				initJsonNode(geneMap, nodes, nodeText);

				Map node = geneMap.get(nodeText);

				if (nsc != null)
				{
					if (nsc.isDownstreamSignificant(sym))
					{
						if (!node.containsKey("css")) node.put("css", new HashMap<>());

						((Map) node.get("css")).put("borderWidth", "2px");
					}
					boolean act = false;
					boolean inh = false;

					if (nsc instanceof NSCForComparison)
					{
						act = ((NSCForComparison) nsc).isActivatingTargetsSignificant(sym);
						inh = ((NSCForComparison) nsc).isInhibitoryTargetsSignificant(sym);
					}

					if (act || inh && !node.containsKey("css")) node.put("css", new HashMap<>());

					if (act && !inh) ((Map) node.get("css")).put("borderColor", inJSONString(activatingBorderColor));
					else if (!act && inh) ((Map) node.get("css")).put("borderColor", inJSONString(inhibitingBorderColor));
					else if (act /** && inh **/) ((Map) node.get("css")).put("borderColor", inJSONString(doubleSignificanceBorderColor));
				}

				if (useGeneBGForTotalProtein && (let.equals("t") || let.equals("c")) && !totalProtUsedUp.contains(nodeText))
				{
					if (!node.containsKey("css")) node.put("css", new HashMap<>());
					((Map) node.get("css")).put("backgroundColor", colS);
					String tooltip = (let.equals("c") ? sym : siteID) + (val.isEmpty() ? "" : ", " + val);
					((Map) node.get("data")).put("tooltipText", tooltip);
					totalProtUsedUp.add(nodeText);
				}
				else
				{
					Map site = new HashMap();
					((List) ((Map) node.get("data")).get("sites")).add(site);

					site.put("siteText", let);
					site.put("siteInfo", siteID);
					site.put("siteBackgroundColor", colS);
					site.put("siteBorderColor", bor);
				}
			}
		});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		JsonUtils.writePrettyPrint(writer, map);
		writer.close();
	}

	/**
	 * This method takes in a SIF graph defined by two files (.sif and .format), and generates a corresponding .json
	 * file that the webserver can display.
	 *
	 * @param sifFileanme SIF filename
	 * @param formatFilename Format filename
	 * @param outJasonFilename JASON filename to produce
	 */
	public static void convertSIFToJSON(String sifFileanme, String formatFilename, String outJasonFilename) throws IOException
	{
		Map<String, Object> map = new HashMap<>();
		List<Map> nodes = new ArrayList<>();
		List<Map> edges = new ArrayList<>();
		map.put("nodes", nodes);
		map.put("edges", edges);

		Map<String, Map> nodeMap = new HashMap<>();
		Set<String> relMem = new HashSet<>();

		Files.lines(Paths.get(sifFileanme)).map(l -> l.split("\t")).forEach(t ->
		{
			if (t.length > 2)
			{
				String key = t[0] + "\t" + t[1] + "\t" + t[2];
				if (relMem.contains(key)) return;
				else relMem.add(key);

				Map<String, Object> edge = new HashMap<>();
				edges.add(edge);
				Map<String, Object> dMap = new HashMap<>();
				edge.put("data", dMap);
				dMap.put("source", t[0]);
				dMap.put("target", t[2]);
				dMap.put("edgeType", t[1]);
				if (t.length > 4 && !t[4].trim().isEmpty())
				{
					dMap.put("tooltipText", t[2] + "-" + CollectionUtil.merge(Arrays.asList(t[4].split(";")), "-"));
				}

				if (t.length > 3 && !t[3].trim().isEmpty())
				{
					List<String> medList = Arrays.asList(t[3].split(";| "));
					if (!medList.isEmpty())
					{
						dMap.put("pcLinks", medList);
					}
				}
				initJsonNode(nodeMap, nodes, t[2]);
			}

			if (t.length > 0 && !t[0].isEmpty()) initJsonNode(nodeMap, nodes, t[0]);
		});

		Map<String, String> defaultColors = new HashMap<>();
		String defBGCKey = "node BG color";
		String defBorCKey = "node border color";

		Files.lines(Paths.get(formatFilename)).map(l -> l.split("\t")).filter(t -> t.length > 3).forEach(t ->
		{
			if (t[1].equals("all-nodes"))
			{
				if (t[2].equals("color")) defaultColors.put(defBGCKey, jasonizeColor(t[3]));
				else if (t[2].equals("bordercolor")) defaultColors.put(defBorCKey, jasonizeColor(t[3]));
			}

			if (t[0].equals("node"))
			{
				String name = t[1];
				Map node = nodeMap.get(name);
				if (node != null)
				{
					if (!node.containsKey("css")) node.put("css", new HashMap<>());

					switch (t[2])
					{
						case "rppasite":
							String[] x = t[3].split("\\|");
							Map site = new HashMap();
							((List) ((Map) node.get("data")).get("sites")).add(site);
							site.put("siteText", x[1]);
							site.put("siteInfo", x[0] + (x.length > 4 ? (" " + x[4]) : ""));
							site.put("siteBackgroundColor", jasonizeColor(x[2]));
							site.put("siteBorderColor", jasonizeColor(x[3]));
							break;
						case "color":
							((Map) node.get("css")).put("backgroundColor", jasonizeColor(t[3]));
							break;
						case "bordercolor":
							((Map) node.get("css")).put("borderColor", jasonizeColor(t[3]));
							break;
						case "borderwidth":
							((Map) node.get("css")).put("borderWidth", t[3] + "px");
							break;
						case "tooltip":
							((Map) node.get("data")).put("tooltipText", t[3]);
							break;
					}
				}
			}
		});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outJasonFilename));
		JsonUtils.writePrettyPrint(writer, map);
		writer.close();
	}

	private static void initJsonNode(Map<String, Map> nodeMap, List<Map> nodes, String name)
	{
		if (!nodeMap.containsKey(name))
		{
			HashMap<String, Object> node = new HashMap<>();
			nodeMap.put(name, node);
			nodes.add(node);
			Map<String, Object> d = new HashMap<>();
			node.put("data", d);
			d.put("id", name);
			d.put("text", name);
			List<Map> sites = new ArrayList<>();
			d.put("sites", sites);
		}
	}

	private static String jasonizeColor(String c)
	{
		String[] t = c.split(" ");
		return "rgb(" + ArrayUtil.getString(",", t[0], t[1], t[2]) + ")";
	}

	private Set<ExperimentData> getExperimentDataToDraw()
	{
		if (experimentDataToDraw != null) return experimentDataToDraw;

		Set<ExperimentData> datas = Stream.concat(
			relations.stream().map(r -> r.sourceData.getChangedData().keySet()).flatMap(Collection::stream),
			relations.stream().map(r -> r.targetData.getChangedData().keySet()).flatMap(Collection::stream))
			.collect(Collectors.toSet());

		if (otherGenesToShow != null)
		{
			datas.addAll(otherGenesToShow.stream().map(GeneWithData::getData).flatMap(Collection::stream)
				.collect(Collectors.toSet()));
		}

		return datas;
	}

	private Set<String> getGenesInGraph()
	{
		return relations.stream().map(r -> new String[]{r.source, r.target}).flatMap(Arrays::stream)
			.collect(Collectors.toSet());
	}
}
