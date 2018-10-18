package org.panda.causalpath.run;

import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class MergeResultsIntoSeriesView
{
	public static final String SIF_FILE = "causative.sif";
	public static final String PARAMETERS_FILE = "parameters.txt";
	public static final String FORMAT_FILE = "causative.format";
	public static final String SERIES_FILE = "causative.formatseries";
	public static final String DEFAULT_OUT_DIR = "AsSeries";

	List<String> inDirs;
	String outDir;

	public MergeResultsIntoSeriesView(List<String> inDirs, String outDir)
	{
		this.inDirs = inDirs;
		this.outDir = outDir;
	}

	public void run() throws IOException
	{

		Set<String> relations = new HashSet<>();
		List<Set<String>> nodeSets = new ArrayList<>();
		List<Set<String>> edgeSets = new ArrayList<>();

		for (String inDir : inDirs)
		{
			Set<String> rels = Files.lines(Paths.get(inDir + File.separator + SIF_FILE)).collect(Collectors.toSet());
			relations.addAll(rels);

			Set<String> nodes = rels.stream().map(l -> l.split("\t")).filter(t -> t.length == 1).map(t -> t[0])
				.collect(Collectors.toSet());

			Set<String> edges = new HashSet<>();

			rels.stream().map(l -> l.split("\t")).filter(t -> t.length >= 3).forEach(t ->
			{
				nodes.add(t[0]);
				nodes.add(t[2]);
				edges.add(t[0] + " " + t[1] + " " + t[2]);
			});

			nodeSets.add(nodes);
			edgeSets.add(edges);
		}

		BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(outDir + File.separator + SIF_FILE));
		relations.forEach(l -> FileUtil.writeln(l, writer1));
		writer1.close();

		List<List<String>> formatLists = new ArrayList<>();
		List<String> stepNames = new ArrayList<>();

		for (String inDir : inDirs)
		{
			formatLists.add(Files.lines(Paths.get(inDir + File.separator + FORMAT_FILE))
				.filter(l -> !l.contains("all-nodes") && !l.contains("all-edges"))
				.collect(Collectors.toList()));

			stepNames.add(inDir.substring(inDir.lastIndexOf(File.separator) + 1));
		}

		// collect representative lines from all sites
		Map<String, String> allSites = new HashMap<>();
		List<Set<String>> sitePresence = new ArrayList<>();

		for (List<String> formatList : formatLists)
		{
			Set<String> siteKeys = new HashSet<>();
			for (String line : formatList)
			{
				if (line.contains("\trppasite\t") && line.contains("|p|"))
				{
					String key = line.substring(0, line.indexOf("|"));
					siteKeys.add(key);

					if (!allSites.containsKey(key))
					{
						allSites.put(key, key + "| |255 255 255|220 220 220");
					}
				}
			}
			sitePresence.add(siteKeys);
		}

		// add ghosted out sites to the frames that don't have those sites
		for (int i = 0; i < formatLists.size(); i++)
		{
			List<String> formatList = formatLists.get(i);
			Set<String> present = sitePresence.get(i);
			for (String key : allSites.keySet())
			{
				if (!present.contains(key))
				{
					formatList.add(allSites.get(key));
				}
			}

			// order sites for consistency across frames
			Collections.sort(formatList);
		}

		List<Set<String>> borderColorSetNodesList = new ArrayList<>();
		for (List<String> formatList : formatLists)
		{
			borderColorSetNodesList.add(formatList.stream().filter(l -> l.contains("\tbordercolor\t"))
				.map(l -> l.split("\t")[1]).collect(Collectors.toSet()));
		}

		BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(outDir + File.separator + SERIES_FILE));

		for (int i = 0; i < stepNames.size(); i++)
		{
			List<String> list = formatLists.get(i);
			writer2.write("group-name\t" + stepNames.get(i) + "\n");
			writer2.write("node\tall-nodes\tcolor\t255 255 255\n");
			writer2.write("node\tall-nodes\tbordercolor\t220 220 220\n");
			writer2.write("node\tall-nodes\tborderwidth\t1\n");
			writer2.write("node\tall-nodes\ttextcolor\t220 220 220\n");
			writer2.write("edge\tall-edges\tcolor\t220 220 220\n");

			list.forEach(l -> FileUtil.writeln(l, writer2));

			for (String node : nodeSets.get(i))
			{
				writer2.write("node\t" + node + "\ttextcolor\t0 0 0\n");

				if (!borderColorSetNodesList.get(i).contains(node))
				{
					writer2.write("node\t" + node + "\tbordercolor\t50 50 50\n");
				}
			}
			for (String edge : edgeSets.get(i))
			{
				SignedType type = SignedType.typeOf(edge.split(" ")[1]);
				writer2.write("edge\t" + edge + "\tcolor\t" + getEdgeColor(type) + "\n");
			}
		}

		writer2.close();

	}

	private String getEdgeColor(SignedType type)
	{
		switch (type)
		{
			case PHOSPHORYLATES:
			case UPREGULATES_EXPRESSION: return "0 150 0";
			case DEPHOSPHORYLATES:
			case DOWNREGULATES_EXPRESSION: return "150 0 0";
			default: return null;
		}
	}

	public static void run(String dir, String... subDirNames) throws IOException
	{
		List<String> dirs;

		if (subDirNames.length > 0)
		{
			dirs = new ArrayList<>();
			for (String subDirName : subDirNames)
			{
				dirs.add(dir + File.separator + subDirName);
			}
		}
		else
		{
			dirs = Files.list(Paths.get(dir)).filter(p ->
				Files.isDirectory(p) &&
					Files.exists(Paths.get(p.toString() + File.separator + SIF_FILE)) &&
					Files.exists(Paths.get(p.toString() + File.separator + PARAMETERS_FILE)))
				.map(Path::toString).collect(Collectors.toList());
		}

		String outDir = dir + File.separator + DEFAULT_OUT_DIR;
		Files.createDirectories(Paths.get(outDir));

		MergeResultsIntoSeriesView app = new MergeResultsIntoSeriesView(dirs, outDir);
		app.run();
	}

	public static void main(String[] args) throws IOException
	{
		run("/home/ozgun/Analyses/CausalPath-paper/EGF-stimulation/run-fdr-0.1-relax",
			"2min", "4min", "8min", "16min", "32min", "64min", "128min");
	}
}
