package org.panda.causalpath.run;

import org.panda.causalpath.analyzer.CausalityHelper;
import org.panda.causalpath.analyzer.CausalitySearcher;
import org.panda.causalpath.analyzer.SignificanceDetector;
import org.panda.causalpath.data.ExperimentData;
import org.panda.causalpath.data.MutationData;
import org.panda.causalpath.network.GraphWriter;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.resource.NetworkLoader;
import org.panda.causalpath.resource.TCGALoader;
import org.panda.resource.tcga.MutSigReader;
import org.panda.utility.Kronometre;
import org.panda.utility.statistics.FDR;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Compares samples that are separated by a gene mutation.
 */
public class TCGAEffectOfMutationRun
{
	public static void main(String[] args) throws IOException
	{
		Kronometre k = new Kronometre();
		String base = "/home/babur/Documents/RPPA/TCGA/effect-of-mutation/";
		String tcgaDataDir = "/home/babur/Documents/TCGA";

		for (File dir : new File(tcgaDataDir).listFiles())
		{
			// DEBUG: run only for one TCGA study
			if (!dir.getName().equals("BRCA")) continue;

			if (Files.exists(Paths.get(dir.getPath() + "/rppa.txt")) &&
				Files.exists(Paths.get(dir.getPath() + "/scores-mutsig.txt")) &&
				Files.exists(Paths.get(dir.getPath() + "/mutation.maf")))
			{
				Map<String, Double> pvals = MutSigReader.readPValues(dir.getPath());
				List<String> select = FDR.select(pvals, null, 0.01);
				select.retainAll(Collections.singleton("PIK3CA"));
				System.out.println("select.size() = " + select.size());
				if (select.isEmpty()) continue;

				File outDir = new File(base + dir.getName());

				Set<Relation> rels = NetworkLoader.load();
				System.out.println("rels.size() = " + rels.size());

				System.out.println("dir = " + dir.getName());

				TCGALoader loader = new TCGALoader(dir.getPath());
				loader.decorateRelations(rels);

				for (String gene : select)
				{
					// DEBUG: run only for one gene mutation
//					if (!gene.equals("PIK3CA")) continue;

					String outFile = outDir.getPath() + "/" + gene + ".sif";

					Set<ExperimentData> data = loader.getData(gene);
					Optional<ExperimentData> opt = data.stream().filter(d -> d instanceof MutationData).findAny();
					if (!opt.isPresent()) continue;

					MutationData md = (MutationData) opt.get();

					boolean[] test = md.getMutated();
					boolean[] control = md.getNotMutated();

					SignificanceDetector det = new SignificanceDetector(0.5, control, test);
					CausalityHelper ch = new CausalityHelper();
					for (Relation rel : rels)
					{
						if (!rel.sourceData.isEmpty() && !rel.targetData.isEmpty())
						{
							rel.sourceData.setChangeDet(det);
							rel.targetData.setChangeDet(det);
							rel.chDet = ch;
						}
					}

					CausalitySearcher searcher = new CausalitySearcher(true);
					Set<Relation> causal = searcher.run(rels);

					System.out.println("causal.size() = " + causal.size());

					if (!causal.isEmpty())
					{
						if (!outDir.exists()) outDir.mkdirs();

						GraphWriter writer = new GraphWriter(causal);
						writer.setUseGeneBGForTotalProtein(true);
						writer.writeSIFGeneCentric(outFile);
					}
				}
			}
		}

		k.stop();
		k.print();
	}
}
