package org.panda.causalpath.network;

import com.github.jsonldjava.utils.JsonUtils;
import org.panda.causalpath.data.*;
import org.panda.causalpath.resource.MutationClassifier;
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
	 * Parameter to use gene background to display the total protein change. If false, then it is displayed on the gene
	 * node as a separate feature.
	 */
	private boolean useGeneBGForTotalProtein;

	/**
	 * An object that can produce a color for a given up/downregulation value.
	 */
	private ValToColor vtc;

	/**
	 * The set of relations with the associated data to draw.
	 */
	private Set<RelationAndSelectedData> relations;

	/**
	 * Constructor with the relations. Those relations are the result of the causality search.
	 */
	public GraphWriter(Set<RelationAndSelectedData> relations)
	{
		this.relations = relations;
		activatingBorderColor = new Color(0, 180, 20);
		inhibitingBorderColor = new Color(180, 0, 20);

		vtc = new ValToColor(new double[]{-1, 0, 1},
			new Color[]{new Color(40, 80, 255), Color.WHITE, new Color(255, 80, 40)});

		useGeneBGForTotalProtein = false;
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
		relations.stream().map(r -> new ExperimentData[]{r.source, r.target}).flatMap(Arrays::stream).distinct()
			.forEach(data ->
		{
			String colS = "255 255 255";
			String bor = "0 0 0";
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
				if (data.getGeneSymbols().stream().anyMatch(MutationClassifier::isActivatedByMutation))
					bor = inString(activatingBorderColor);
				else bor = inString(inhibitingBorderColor);
			}
			else if (data instanceof CNAData)
			{
				let = "c";
			}
			else if (data instanceof ExpressionData)
			{
				let = "r";
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

		relations.stream().map(r -> new ExperimentData[]{r.source, r.target}).flatMap(Arrays::stream).distinct()
			.forEach(data ->
			{
				String colS = "255 255 255";
				String bor = "0 0 0";
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
					if (data.getGeneSymbols().stream().anyMatch(MutationClassifier::isActivatedByMutation))
						bor = inJSONString(activatingBorderColor);
					else bor = inJSONString(inhibitingBorderColor);
				}
				else if (data instanceof CNAData)
				{
					let = "c";
				}
				else if (data instanceof ExpressionData)
				{
					let = "r";
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
}
