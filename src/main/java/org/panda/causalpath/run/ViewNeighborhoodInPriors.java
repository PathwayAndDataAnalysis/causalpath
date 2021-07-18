package org.panda.causalpath.run;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * This class takes in a SIF file, a format file, and several gene names, and generates a neighborhood graph for the
 * given gene names in the given SIF file, and overlays data features in the format file. The view can be upload to
 * webserver. The purpose of this visualization is to understand what priors exist in the neighborhood of certain genes
 * and what data they have on. The view does not include any site matching or causality checks.
 */
public class ViewNeighborhoodInPriors
{
	public static void main(String[] args) throws IOException
	{
		switch (args[0])
		{
			case "neighborhood":
			{
				String priorFile = args[1];
				String formatFile = args[2];
				String outDir = args[3];
				Set<String> seed = new HashSet<>(Arrays.asList(args).subList(4, args.length));

				generateNeighborhood(priorFile, formatFile, outDir, seed);
				break;
			}
			case "paths-between":
			{
				String priorFile = args[1];
				String formatFile = args[2];
				String outDir = args[3];

				generatePathsBetween(priorFile, formatFile, outDir);
				break;
			}
		}
	}

	private static void generateNeighborhood(String priorFile, String formatFile, String outDir, Set<String> seed) throws IOException
	{
		String subsetSIF = outDir + "/neighborhood.sif";
		String subsetFormat = outDir + "/neighborhood.format";
		Files.createDirectories(Paths.get(outDir));
		Set<String> nodes = new HashSet<>();
		Set<String> edges = new HashSet<>();
		BufferedWriter sifWriter = FileUtil.newBufferedWriter(subsetSIF);
		FileUtil.lines(priorFile).forEach(l ->
		{
			String[] t = l.split("\t");
			if (seed.contains(t[0]) || seed.contains(t[2]))
			{
				FileUtil.writeln(l, sifWriter);
				nodes.add(t[0]);
				nodes.add(t[2]);
				edges.add(t[0] + " " + t[1] + " " + t[2]);
			}
		});
		sifWriter.close();

		writeFormatFile(formatFile, outDir, subsetFormat, nodes, edges);
	}

	private static void generatePathsBetween(String priorFile, String formatFile, String outDir) throws IOException
	{
		Set<String> seed = FileUtil.lines(formatFile).map(l -> l.split("\t"))
			.filter(t -> t.length > 2 && t[0].equals("node") && !t[1].equals("all-nodes"))
			.map(t -> t[1]).collect(Collectors.toSet());
		System.out.println("seed.size() = " + seed.size());
		generatePathsBetween(priorFile, formatFile, outDir, seed);
	}
	private static void generatePathsBetween(String priorFile, String formatFile, String outDir, Set<String> seed) throws IOException
	{
		String subsetSIF = outDir + "/paths-between.sif";
		String subsetFormat = outDir + "/paths-between.format";
		Files.createDirectories(Paths.get(outDir));
		Set<String> nodes = new HashSet<>();
		Set<String> edges = new HashSet<>();
		BufferedWriter sifWriter = FileUtil.newBufferedWriter(subsetSIF);
		FileUtil.lines(priorFile).forEach(l ->
		{
			String[] t = l.split("\t");
			if (seed.contains(t[0]) && seed.contains(t[2]))
			{
				FileUtil.writeln(l, sifWriter);
				nodes.add(t[0]);
				nodes.add(t[2]);
				edges.add(t[0] + " " + t[1] + " " + t[2]);
			}
		});
		sifWriter.close();

		writeFormatFile(formatFile, outDir, subsetFormat, nodes, edges);
	}

	private static void writeFormatFile(String formatFile, String outDir, String subsetFormat, Set<String> nodes, Set<String> edges) throws IOException
	{
		BufferedWriter formatWriter = FileUtil.newBufferedWriter(subsetFormat);
		FileUtil.lines(formatFile).forEach(l ->
		{
			String[] t = l.split("\t");

			if (t[1].equals("all-nodes") || t[1].equals("all-edges") ||
				(t[0].equals("node") && nodes.contains(t[1])) ||
				(t[0].equals("edge") && edges.contains(t[1])))
			{
				FileUtil.writeln(l, formatWriter);
			}
		});
		formatWriter.close();

		String graphName = subsetFormat.substring(subsetFormat.lastIndexOf(File.separator) + 1, subsetFormat.lastIndexOf("."));
		JasonizeResultGraphsRecursively.generate(outDir, outDir, Collections.singleton(graphName), outDir,
			"causative.json");
	}
}
