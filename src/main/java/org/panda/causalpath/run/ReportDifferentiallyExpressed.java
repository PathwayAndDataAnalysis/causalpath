package org.panda.causalpath.run;

import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.ProteinData;
import org.panda.causalpath.network.Relation;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Ozgun Babur
 */
public class ReportDifferentiallyExpressed
{
	public static void report(Set<Relation> relations, String outFile) throws IOException
	{
		Set<ExperimentData> datas = relations.stream().map(Relation::getAllData).flatMap(Collection::stream)
			.collect(Collectors.toSet());

		Set<String> genes = new HashSet<>();

		for (ExperimentData data : datas)
		{
			if (data instanceof ProteinData)
			{
				if (data.getChangeSign() != 0) genes.addAll(data.getGeneSymbols());
			}
		}

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		genes.stream().sorted().forEach(g -> FileUtil.writeln(g, writer));
		writer.close();
	}
}
