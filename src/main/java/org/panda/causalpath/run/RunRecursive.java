package org.panda.causalpath.run;

import org.panda.utility.FileUtil;
import org.panda.utility.RunUtil;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 * Executes CausalPath recursively starting from the given directory, and navigating through subdirectories.
 *
 * @author Ozgun Babur
 */
public class RunRecursive
{
	public static void main(String[] args) throws IOException, ClassNotFoundException
	{
		if (args.length > 1)
		{
			runMultiCore(args[0], Integer.parseInt(args[1]));
		}
		else
		{
			run(args[0]);
		}
	}

	private static void run(String dir) throws IOException, ClassNotFoundException
	{
		System.out.println("runRecursive dir = " + dir);
		if (Files.isDirectory(Paths.get(dir)))
		{
			if (Files.exists(Paths.get(dir + File.separator + CausalPath.PARAMETER_FILENAME)))
//				&& !Files.exists(Paths.get(dir + File.separator + CausalPath.CAUSATIVE_RESULT_FILE_PREFIX + ".sif")))
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


	private static void runMultiCore(String dir, int cores) throws IOException, ClassNotFoundException
	{
		List<Runnable> runnables = new ArrayList<>();

		FileUtil.processDirsRecursive(new File(dir), curDir ->
		{
			if (FileUtil.exists(curDir.getPath() + File.separator + CausalPath.PARAMETER_FILENAME))
			{
				runnables.add(new Runnable()
				{
					@Override
					public void run()
					{
						System.gc();
						try
						{
							CausalPath.main(new String[]{curDir.getPath()});
						}
						catch (IOException | ClassNotFoundException e)
						{
							e.printStackTrace();
						}
					}
				});
			}
		});

		RunUtil.run(runnables, cores);
	}



}
