package org.panda.causalpath.run;

import org.panda.causalpath.network.GraphWriter;
import org.panda.utility.FileUtil;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * This class is used for converting SIF and format files to json files that can be uploaded to the webserver for
 * visualization.
 *
 * Usage: java -jar causalpath.jar <in-dir> <sif-name> <out-dir> <json-name>
 *
 */
public class JasonizeResultGraphsRecursively
{
	public static void main(String[] args) throws IOException
	{
		String inBase = new File(args[0]).getPath();
		String sifNameWOExtension = args[1];
		String outBase = new File(args[2]).getPath();
		String jsonName = args.length > 3 ? args[3] : "causative.json";

		generate(inBase, inBase, sifNameWOExtension, outBase, jsonName);
	}

	private static void generate(String inBase, String inDir, String sifName, String outBase, String jsonName) throws IOException
	{
		String sifPath = inDir + File.separator + sifName + ".sif";
		String formatPath = inDir + File.separator + sifName + ".format";

		if (Files.exists(Paths.get(sifPath)) && !FileUtil.isEmpty(sifPath) && Files.exists(Paths.get(formatPath)))
		{
			String outDir = inDir.replace(inBase, outBase);
			Files.createDirectories(Paths.get(outDir));

			GraphWriter.convertSIFToJSON(sifPath, formatPath, outDir + File.separator + jsonName);
		}

		for (File sub : new File(inDir).listFiles())
		{
			if (sub.isDirectory()) generate(inBase, sub.getPath(), sifName, outBase, jsonName);
		}
	}
}
