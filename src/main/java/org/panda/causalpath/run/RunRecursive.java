package org.panda.causalpath.run;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Executes CausalPath recursively starting from the given directory, and navigating through subdirectories.
 *
 * @author Ozgun Babur
 */
public class RunRecursive
{
	public static void main(String[] args) throws IOException, ClassNotFoundException
	{
		run(args[0]);
	}

	private static void run(String dir) throws IOException, ClassNotFoundException
	{
		if (Files.isDirectory(Paths.get(dir)))
		{
			if (Files.exists(Paths.get(dir + File.separator + CausalPath.PARAMETER_FILENAME))
				&&
				!Files.exists(Paths.get(dir + File.separator + CausalPath.CAUSATIVE_RESULT_FILE_PREFIX + ".sif")))
			{
				System.gc();
				CausalPath.main(new String[]{dir});
				System.out.println();
			}
			else
			{
				for (File file : new File(dir).listFiles())
				{
					if (file.isDirectory())
					{
						run(file.getPath());
					}
				}
			}
		}
	}
}
