package org.panda.causalpath.resource;

import org.panda.resource.tcga.MutationReader;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Created by babur on 4/14/16.
 */
public class MutationClassifier
{
	public static final String FILENAME = "mutation-stats.txt";
	private static Set<String> activatedByMutation;

	public static boolean isActivatedByMutation(String gene)
	{
		return activatedByMutation.contains(gene);
	}

	private static void prepareResource() throws IOException
	{
		MutationReader reader = new MutationReader(null);

		Arrays.stream(new File("/home/babur/Documents/TCGA").listFiles())
			.filter(file -> new File(file.getPath() + "/mutation.maf").exists())
			.filter(file -> !file.getName().equals("PAAD"))
			.forEach(file ->
			{
				try{
					reader.load(file.getPath() + "/mutation.maf", null);
				} catch (IOException e){throw new RuntimeException(e);}
			});

		Map<String, Integer> recCnt = reader.getHighestRecurrenceCounts();
		Map<String, Double> delRat = reader.getRatiosOfDeleteriousMutations();

		BufferedWriter writer = new BufferedWriter(
			new FileWriter("src/main/resources/org/babur/causalpath/resource/" + FILENAME));

		recCnt.keySet().stream().filter(gene -> recCnt.get(gene) > 1)
			.sorted(Comparator.comparing(recCnt::get).reversed()).forEach(gene ->
			{
				try{
					writer.write(gene + "\t" + recCnt.get(gene) + "\t" + delRat.get(gene) + "\n");
				} catch (IOException e){throw new RuntimeException(e);}
			});

		writer.close();
	}

	static
	{
		try
		{
			activatedByMutation = Files.lines(Paths.get(MutationClassifier.class.getResource(FILENAME).toURI()))
				.map(line -> line.split("\t"))
				.filter(token -> Integer.parseInt(token[1]) > 5 && Double.parseDouble(token[2]) < 0.1)
				.map(token -> token[0]).collect(Collectors.toSet());
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	public static void main(String[] args) throws IOException
	{
		prepareResource();
	}
}
