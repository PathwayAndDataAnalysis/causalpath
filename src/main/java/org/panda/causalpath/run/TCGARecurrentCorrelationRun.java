package org.panda.causalpath.run;

import org.panda.causalpath.analyzer.CausalitySearcher;
import org.panda.causalpath.analyzer.CorrelationDetector;
import org.panda.causalpath.network.GraphWriter;
import org.panda.causalpath.network.Relation;
import org.panda.causalpath.network.RelationAndSelectedData;
import org.panda.causalpath.resource.NetworkLoader;
import org.panda.causalpath.resource.TCGALoader;
import org.panda.utility.Kronometre;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * For finding recurrent causative correlations in TCGA RPPA data.
 */
public class TCGARecurrentCorrelationRun
{
	public static void main(String[] args) throws IOException
	{
		Kronometre k = new Kronometre();
		String base = "/home/babur/Documents/RPPA/TCGA/basic-correlation-rppa-mut/";
		String single = base + "single/";
		String recurrent = base + "recurrent/";
		String tcgaDataDir = "/home/babur/Documents/TCGA";

		if (!(new File(single).exists())) new File(single).mkdirs();

		Map<RelationAndSelectedData, Integer> relCnt = new HashMap<>();

		for (File dir : new File(tcgaDataDir).listFiles())
		{
			if (dir.getName().equals("PanCan")) continue;

			if (new File(dir.getPath() + "/rppa.txt").exists())
			{
				String outFile = single + dir.getName() + ".sif";

				Set<Relation> rels = NetworkLoader.load();
				System.out.println("rels.size() = " + rels.size());

				System.out.println("dir = " + dir.getName());

				TCGALoader loader = new TCGALoader(dir.getPath());
				loader.decorateRelations(rels);

				CorrelationDetector det = new CorrelationDetector(0.01, 0.01);
				for (Relation rel : rels)
				{
					if (!rel.sourceData.isEmpty() && !rel.targetData.isEmpty())
					{
						rel.chDet = det;
					}
				}
				CausalitySearcher searcher = new CausalitySearcher();
//				searcher.setGenesWithTotalProteinData(loader.getGenesWithTotalProteinData());
				Set<RelationAndSelectedData> causal = searcher.run(rels);

				for (RelationAndSelectedData r : causal)
				{
					relCnt.put(r, relCnt.containsKey(r) ? relCnt.get(r) + 1 : 1);
				}

				System.out.println("causal.size() = " + causal.size());
				GraphWriter writer = new GraphWriter(causal);
				writer.setUseGeneBGForTotalProtein(false);
				writer.writeSIFGeneCentric(outFile);
			}
		}

		if (!(new File(recurrent).exists())) new File(recurrent).mkdirs();
		Map<Integer, Set<RelationAndSelectedData>> grouped = new HashMap<>();
		Set<Integer> rec = new HashSet<>(relCnt.values());
		for (Integer i : rec)
		{
			Set<RelationAndSelectedData> set = new HashSet<>();
			for (RelationAndSelectedData r : relCnt.keySet())
			{
				if (relCnt.get(r) >= i) set.add(r);
			}
			grouped.put(i, set);
		}

		for (Integer i : grouped.keySet())
		{
			GraphWriter writer = new GraphWriter(grouped.get(i));
			writer.setUseGeneBGForTotalProtein(false);
			writer.writeSIFGeneCentric(recurrent + i + ".sif");
		}

		k.stop();
		k.print();
	}
}
