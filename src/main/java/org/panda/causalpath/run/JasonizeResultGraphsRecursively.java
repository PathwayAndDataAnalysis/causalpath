package org.panda.causalpath.run;

import org.panda.causalpath.network.GraphWriter;
import org.panda.utility.FileUtil;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Set;

/**
 * This class is used for converting SIF and format files to json files that can be uploaded to the webserver for
 * visualization.
 *
 * Usage: java -jar causalpath.jar <in-dir> <out-dir> <sif-name-no-ext>
 *
 */
public class JasonizeResultGraphsRecursively
{
	public static void main(String[] args) throws IOException
	{
		String inBase = new File(args[0]).getPath();
		String outBase = new File(args[1]).getPath();
		Set<String> sifNames = new HashSet<>();
		for (int i = 2; i < args.length; i++)
		{
			sifNames.add(args[i]);
		}
		String jsonName = "causative.json";

		generate(inBase, inBase, sifNames, outBase, jsonName);
	}

	public static void generate(String inBase, String inDir, Set<String> sifNames, String outBase, String jsonName) throws IOException
	{
		for (String sifName : sifNames)
		{
			String sifPath = inDir + File.separator + sifName + ".sif";
			String formatPath = inDir + File.separator + sifName + ".format";

			if (Files.exists(Paths.get(sifPath)) && !FileUtil.isEmpty(sifPath) && Files.exists(Paths.get(formatPath)))
			{
				String outDir = inDir.replace(inBase, outBase);
				if (sifNames.size() > 1) outDir += File.separator + sifName;
				Files.createDirectories(Paths.get(outDir));

				GraphWriter.convertSIFToJSON(sifPath, formatPath, outDir + File.separator + jsonName);
			}
		}

		for (File sub : new File(inDir).listFiles())
		{
			if (sub.isDirectory())
			{
				//--- Temporary hack
//				if (!sub.getPath().contains("correlation") && (sub.getName().equals("sitespec") || sub.getName().equals("rnaseq") || sub.getName().equals("totprot"))) continue;
//				if (!sub.getPath().contains("correlation") && sub.getName().equals("all")) generate(inBase, inDir, sifNames, outBase, jsonName);
//				else
				//--- Temporary hack

				generate(inBase, sub.getPath(), sifNames, outBase, jsonName);
			}
		}
	}
}
