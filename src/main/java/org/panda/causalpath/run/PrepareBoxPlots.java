package org.panda.causalpath.run;

import org.panda.causalpath.analyzer.ComparisonDetector;
import org.panda.causalpath.analyzer.OneDataChangeDetector;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.NumericData;
import org.panda.causalpath.network.Relation;
import org.panda.utility.statistics.BoxPlot;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * This is a utility class for visualizing the data distribution in a box plot.
 *
 * @author Ozgun Babur
 */
public class PrepareBoxPlots
{
	public static final String DIR = "/plots/";

	public static void runWithRelations(Set<Relation> relations, String directory, String controlName, String testName) throws IOException
	{
		Set<ExperimentData> datas = Stream.concat(relations.stream().map(r -> r.sourceData),
			relations.stream().map(r -> r.targetData)).flatMap(Collection::stream).collect(Collectors.toSet());

		run(datas, directory, controlName, testName);
	}

	public static void run(Set<ExperimentData> datas, String directory, String controlName, String testName) throws IOException
	{
		for (ExperimentData data : datas)
		{
			if (data instanceof NumericData)
			{
				OneDataChangeDetector chDet = data.getChDet();

				if (chDet instanceof ComparisonDetector)
				{
					String clazz = data.getClass().getName();
					NumericData numDat = (NumericData) data;
					ComparisonDetector comDet = (ComparisonDetector) chDet;

					String dir = directory + DIR + clazz;
					Files.createDirectories(Paths.get(dir));

					String filename = dir + File.separator + numDat.getId() + ".txt";
					List<Double> ctrl = getSubset(numDat.vals, comDet.getControl());
					List<Double> test = getSubset(numDat.vals, comDet.getTest());

					BoxPlot.write(filename, new String[]{controlName, testName}, new List[]{ctrl, test});
				}
			}
		}
	}

	private static List<Double> getSubset(double[] vals, boolean[] select)
	{
		List<Double> list = new ArrayList<>();
		for (int i = 0; i < vals.length; i++)
		{
			if (select[i]) list.add(vals[i]);
		}
		return list;
	}
}
