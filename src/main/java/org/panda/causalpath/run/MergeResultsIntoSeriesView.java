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
	public static final String SIF_FILE = "oncogenic-changes-lung.sif";
	public static final String PARAMETERS_FILE = "parameters.txt";
	public static final String FORMAT_FILE = "oncogenic-changes-lung.format";
	public static final String SERIES_FILE = "oncogenic-changes-lung.formatseries";
	public static final String DEFAULT_OUT_DIR = "AsSeries";

	List<String> inDirs;
	String outDir;

	public MergeResultsIntoSeriesView(List<String> inDirs, String outDir)
	{
		if (inDirs.size() < 2)
		{
			throw new RuntimeException("Need at least two input folders. Found only " + inDirs.size() + ".");
		}

		this.inDirs = inDirs;
		this.outDir = outDir;
	}

	public void runFlatFolders() throws IOException
	{
		Files.createDirectories(Paths.get(outDir));

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
		for (String inDir : inDirs)
		{
			formatLists.add(Files.lines(Paths.get(inDir + File.separator + FORMAT_FILE))
				.filter(l -> !l.contains("all-nodes") && !l.contains("all-edges"))
				.collect(Collectors.toList()));
		}

		List<String> stepNames = getStepNames(inDirs);

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

	/**
	 * Generates step names by removing the common prefixes and suffixes of the inDirs.
	 * @param inDirs
	 * @return
	 */
	private List<String> getStepNames(List<String> inDirs)
	{
		String prefix = getCommon(inDirs, true);
		String suffix = getCommon(inDirs, false);

		return inDirs.stream().map(s -> s.substring(prefix.length(), s.length() - suffix.length()))
			.collect(Collectors.toList());
	}

	private String getCommon(List<String> inDirs, boolean prefix)
	{
		String s1 = inDirs.get(0);

		for (int i = 1; i < inDirs.size(); i++)
		{
			String s2 = inDirs.get(i);
			s1 = prefix ? getCommonPrefix(s1, s2) : getCommonSuffix(s1, s2);
		}

		return s1;
	}

	private String getCommonPrefix(String s1, String s2)
	{
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < Math.min(s1.length(), s2.length()); i++)
		{
			if (s1.charAt(i) == s2.charAt(i)) sb.append(s1.charAt(i));
			else break;
		}
		return sb.toString();
	}
	private String getCommonSuffix(String s1, String s2)
	{
		StringBuilder sb = new StringBuilder();
		for (int i = 1; i <= Math.min(s1.length(), s2.length()); i++)
		{
			if (s1.charAt(s1.length() - i) == s2.charAt(s2.length() - i)) sb.append(s1.charAt(s1.length() - i));
			else break;
		}
		return sb.reverse().toString();
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

	public static void runFlatFolders(String dir, String... subDirNames) throws IOException
	{
		List<String> inDirs;

		if (subDirNames.length == 0)
		{
			inDirs = Files.list(Paths.get(dir)).filter(p ->
				Files.isDirectory(p) &&
					Files.exists(Paths.get(p.toString() + File.separator + SIF_FILE)) &&
					Files.exists(Paths.get(p.toString() + File.separator + PARAMETERS_FILE)))
				.map(Path::toString).collect(Collectors.toList());
		}
		else
		{
			inDirs = Arrays.stream(subDirNames).map(d -> dir + File.separator + d).collect(Collectors.toList());
		}

		if (inDirs.size() < 2)
		{
			throw new RuntimeException("Need at least two input folders. Found only " + inDirs.size() + ".");
		}

		MergeResultsIntoSeriesView app = new MergeResultsIntoSeriesView(inDirs, dir + File.separator + DEFAULT_OUT_DIR);
		app.runFlatFolders();
	}

	public static void run(String outDir, String... inDirs) throws IOException
	{
		if (inDirs.length < 2)
		{
			throw new RuntimeException("Need at least two input folders. Found only " + inDirs.length + ".");
		}

		Files.createDirectories(Paths.get(outDir));

		MergeResultsIntoSeriesView app = new MergeResultsIntoSeriesView(Arrays.asList(inDirs), outDir);
		app.runFlatFolders();
	}

	public static void main(String[] args) throws IOException
	{
		runFlatFolders("/home/ozgunbabur/Analyses/CPTAC-LSCC/v3/tumors-vs-normals-per-NMF-subtype/",
			"Metabolic.proliferative", "NOTCH2.immune.cold", "TP53WT.hot.quiescent", "ROS.hot", "EMT.angio.quiescent", "TP63SOX2amp.KEAP1.prolif.cold");
	}
}
