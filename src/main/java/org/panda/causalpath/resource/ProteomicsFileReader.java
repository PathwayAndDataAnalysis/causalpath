package org.panda.causalpath.resource;

import org.panda.resource.tcga.ProteomicsFileRow;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * This class reads a proteomics file that is specially formatted to use in this framework. This is a tab-delimited data
 * file with first row having column headers, and with following columns.
 *
 * ID: An ID for the row. A descriptive ID is preferable since it is displayed on the final interactive graph.
 * Symbols: Gene symbols associated with the row. There should be a single space character between each symbol if there
 *          is more than one.
 * Sites: The sites on the protein related to this measurement. If there are more than one site, they should be
 *          separated by a pipe (|), and if there is more than one genes for the row, sites of each gene should be
 *          separated with a space character. Each site is a string starting with a one letter amino acid code and
 *          continuing with the integer position, such like S473.
 * Effect: The effect that the site(s) make on the protein. It can be activating (a), inactivating (i), or more complex
 *          (c). Since the causality framework cannot use complex effects, it is preferable to choose either a or i to
 *          use that row in analysis. If this column is left blank, then PhosphoSitePlus data is used to fill it.
 * Value(s): This potentially corresponds to more than one columns, and holds the measurement value(s) for the row.
 *
 * The column names can be customized. The header names are passed as parameter to the related methods.
 */
public class ProteomicsFileReader
{
	protected final static double LOG2 = Math.log(2);

	/**
	 * Reads the annotation in a given proteomics file.
	 *
	 * @param filename name of the file
	 * @param idname name of the ID column
	 * @param symbolname name of the symbols column
	 * @param psitename name of the sites column
	 * @param effectName name of the effect column
	 */
	public static List<ProteomicsFileRow> readAnnotation(String filename, String idname, String symbolname,
		String psitename, String effectName) throws FileNotFoundException
	{
		List<ProteomicsFileRow> datas = new ArrayList<>();

		Scanner sc = new Scanner(new File(filename));
		String s = sc.nextLine();
		List<String> cols = Arrays.asList(s.split("\t"));

		int colInd = cols.indexOf(idname);
		int symbolInd = cols.indexOf(symbolname);
		int siteInd = cols.indexOf(psitename);
		int effectInd = effectName == null ? -1 : cols.indexOf(effectName);

		while (sc.hasNextLine())
		{
			String[] row = sc.nextLine().split("\t");
			String id = row[colInd];
			String syms = row[symbolInd];
			String sites = row.length > siteInd ? row[siteInd] : "";
			String effect = effectInd >= 0 && row.length > effectInd ? row[effectInd] : null;

			List<String> genes = Arrays.asList(syms.split("\\s+"));
			Map<String, List<String>> siteMap = sites.isEmpty() ? null : new HashMap<>();
			if (!sites.isEmpty())
			{
				String[] perGene = sites.split("\\s+");
				for (int i = 0; i < perGene.length; i++)
				{
					siteMap.put(genes.get(i), Arrays.asList(perGene[i].split("\\|")));
				}
				if (siteMap.size() < genes.size())
				{
					for (int i = siteMap.size(); i < genes.size(); i++)
					{
						siteMap.put(genes.get(i), siteMap.get(genes.get(0)));
					}
				}
			}

			ProteomicsFileRow data = new ProteomicsFileRow(id, null, genes, siteMap);

			if (effect != null)
			{
				data.effect = effect.equals("c") ? ProteomicsFileRow.SiteEffect.COMPLEX :
					effect.equals("a") ? ProteomicsFileRow.SiteEffect.ACTIVATING : effect.equals("i") ?
						ProteomicsFileRow.SiteEffect.INHIBITING : null;
			}

			datas.add(data);
		}
		return datas;
	}

	/**
	 * Reads the measurement values in the given proteomics file.
	 *
	 * @param filename File name
	 * @param idName Name of the ID column
	 * @param colname Array of the names of the value columns
	 * @return Array of values for each row ID
	 */
	public static Map<String, Double>[] readVals(String filename, String idName, String... colname)
	{
		try
		{
			Map<String, Double>[] valMaps = new Map[colname.length];
			for (int i = 0; i < valMaps.length; i++)
			{
				valMaps[i] = new HashMap<>();
			}

			Scanner sc = new Scanner(new File(filename));
			List<String> header = Arrays.asList(sc.nextLine().split("\t"));

			int idInd = header.indexOf(idName);
			int[] valInd = new int[colname.length];
			for (int i = 0; i < colname.length; i++)
			{
				valInd[i] = header.indexOf(colname[i]);
			}

			while (sc.hasNextLine())
			{
				String[] row = sc.nextLine().split("\t");

				for (int i = 0; i < colname.length; i++)
				{
					double val;
					try { val = Double.parseDouble(row[valInd[i]]); }
					catch (NumberFormatException e){val = Double.NaN;}

					valMaps[i].put(row[idInd], val);
				}
			}

			return valMaps;
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * After reading the annotations from a proteomics file, this method is used to fill in the measurement values.
	 * These two methods are not in the same method because annotation can be in a separate file or it can appear in
	 * with the values. Presence of two methods satisfies this flexibility.
	 * @param datas rows with annotations
	 * @param filename file name
	 * @param idColName name of the id column
	 * @param vals list of value columns
	 * @param missingVal what to use if a value is missing
	 */
	public static void addValues(List<ProteomicsFileRow> datas, String filename, String idColName,
		List<String> vals, Double missingVal, boolean logTransform)
	{
		Map<String, Double>[] v = readVals(
			filename, idColName, vals.toArray(new String[vals.size()]));

		for (ProteomicsFileRow data : datas)
		{
			data.vals = new double[v.length];
			for (int i = 0; i < v.length; i++)
			{
				Double doubVal = v[i].get(data.id);
				if (doubVal != null)
				{
					data.vals[i] = logTransform ? Math.log(doubVal) / LOG2 : doubVal;
				}
				else if (missingVal == null) data.vals[i] = Double.NaN;
				else data.vals[i] = missingVal;
			}
		}
	}

	public static int getPotentialIDColIndex(String[] header)
	{
		return getPotentialColIndex(header, "id");
	}

	public static int getPotentialSymbolColIndex(String[] header)
	{
		return getPotentialColIndex(header, "symbol");
	}

	public static int getPotentialSiteColIndex(String[] header)
	{
		return getPotentialColIndex(header, "site");
	}

	public static int getPotentialEffectColIndex(String[] header)
	{
		return getPotentialColIndex(header, "effect");
	}

	/**
	 * Searches certain words in column headers to infer which column contains which type of information.
	 */
	private static int getPotentialColIndex(String[] header, String find)
	{
		for (int i = 0; i < header.length; i++)
		{
			if (header[i].toLowerCase().contains(find)) return i;
		}
		return -1;
	}

	/**
	 * Checks which columns contain numbers to infer value columns.
	 */
	public static List<String> getNamesOfNumberColumns(String filename)
	{
		String[] header = getARow(filename, 1);
		String[] row = getARow(filename, 2);

		List<String> names = new ArrayList<>();

		for (int i = 0; i < row.length; i++)
		{
			try
			{
				Double.parseDouble(row[i]);
				names.add(header[i]);
			}
			catch (NumberFormatException e){}
		}

		return names;
	}

	/**
	 * Reads the requested row from the file.
	 */
	public static String[] getARow(String filename, int rowNum)
	{
		try
		{
			return Files.lines(Paths.get(filename)).skip(rowNum - 1).findFirst().get().split("\t");
		}
		catch (IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}
}
