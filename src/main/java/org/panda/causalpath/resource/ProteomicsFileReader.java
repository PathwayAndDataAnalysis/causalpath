package org.panda.causalpath.resource;

import htsjdk.tribble.index.Index;
import org.panda.causalpath.log.CPLogger;
import org.panda.resource.siteeffect.Feature;
import org.panda.resource.tcga.ProteomicsFileRow;
import org.slf4j.Logger;

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
	 * @param idCol name of the ID column
	 * @param symbolCol name of the symbols column
	 * @param siteCol name of the sites column
	 * @param effectCol name of the effect column
	 */

	/**
	 * Reads the annotation in a given proteomics file.
	 *
	 * @param filename name of the file
	 * @param idCol name of the ID column
	 * @param symbolCol name of the symbols column
	 * @param siteCol name of the sites column
	 * @param effectCol name of the effect column
	 */
	public static List<ProteomicsFileRow> readAnnotation(String filename, String idCol, String symbolCol,
		String siteCol, String modCol, String effectCol) throws FileNotFoundException
	{
		List<ProteomicsFileRow> datas = new ArrayList<>();

		Scanner sc = new Scanner(new File(filename));
		String s = sc.nextLine();
		List<String> cols = Arrays.asList(s.split("\t"));

		int colInd = cols.indexOf(idCol);
		int symbolInd = cols.indexOf(symbolCol);
		int siteInd = cols.indexOf(siteCol);
		int modInd = modCol == null ? -1 : cols.indexOf(modCol);
		int effectInd = effectCol == null ? -1 : cols.indexOf(effectCol);

		String fileNameForError = filename.substring(filename.lastIndexOf("/") + 1);

		if (colInd == -1 && CPLogger.isInitialized) CPLogger.dataError.error("ID column '" + idCol + "' not found in data file '" + fileNameForError + "'.");
		if (symbolInd == -1 && CPLogger.isInitialized) CPLogger.dataError.error("Symbol column '" + symbolCol + "' not found in data file '" + fileNameForError + "'.");
		if (siteInd == -1 && CPLogger.isInitialized) CPLogger.dataError.error("Site column '" + siteCol + "' not found in data file '" + fileNameForError + "'.");
		if (modInd == - 1 && modCol != null && CPLogger.isInitialized)
			CPLogger.dataError.error("Provided feature column '" + modCol + "' but column was not found in the data file '" + fileNameForError + "'.");
		if (effectInd == - 1 && effectCol != null && CPLogger.isInitialized)
			CPLogger.dataError.error("Provided effect column '" + effectCol + "' but column was not found in the data file '" + fileNameForError + "'.");

		int line = 1;
		while (sc.hasNextLine())
		{
			line++;
			String[] row = sc.nextLine().split("\t");
			String id = row[colInd];
			String syms = row[symbolInd];
			String sites = row.length > siteInd ? row[siteInd] : "";
			Feature feature = null;
			String effect = null;

			if (modInd >= 0 && row.length > modInd)
			{
				Feature f = Feature.getFeat(row[modInd]);
				if (CPLogger.isInitialized && f == null)
				{
					CPLogger.dataWarning.warn("Line " + line + ": Feature '" + row[modInd] + "' is not recognized. Defaulting to no feature." +
							"Possible features are: 'P', 'A', 'M', 'U', 'G', 'R', 'C'");
				}
				else
				{
					feature = f;
				}
			}

			if (effectInd >= 0 && row.length > effectInd)
			{
				if (row[effectInd].equals("i") || row[effectInd].equals("a") || row[effectInd].equals("c") || row[effectInd].equals(""))
				{
					effect = row[effectInd];
				}
				else if (CPLogger.isInitialized)
				{
					CPLogger.dataWarning.warn("Line " + line + ": Effect '" + row[effectInd] + "' is not recognized. Defaulting to no effect." +
							"Possible effects are: 'i', 'a', 'c', or leaving the field blank");
				}
			}

			List<String> genes = Arrays.asList(syms.split("\\s+"));
			Map<String, List<String>> siteMap = sites.isEmpty() ? null : new HashMap<>();
			if (!sites.isEmpty())
			{
				String[] perGene = sites.split("\\s+");
				for (int i = 0; i < perGene.length; i++)
				{
					if (genes.size() <= i)
					{
						System.out.println("id: " + id + " genes: " + genes);
					}
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

			ProteomicsFileRow data = new ProteomicsFileRow(id, null, genes, siteMap, feature);

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

			if (idInd < 0)
			{
				if (CPLogger.isInitialized) CPLogger.dataError.error("ID column '" + idName + "' not found in data file '" + filename + "'.");
				throw new RuntimeException("Cannot find \"" + idName + "\" column in values file.");
			}

			int[] valInd = new int[colname.length];
			for (int i = 0; i < colname.length; i++)
			{
				valInd[i] = header.indexOf(colname[i]);
				if (valInd[i] == -1 && CPLogger.isInitialized)
				{
					CPLogger.dataError.error("Value column '" + colname[i] + "' not found in data file '" + filename + "'.");
				}
			}

			int line = 1;
			while (sc.hasNextLine())
			{
				line++;
				String[] row = sc.nextLine().split("\t");

				for (int i = 0; i < colname.length; i++)
				{
					double val;
					try { val = Double.parseDouble(row[valInd[i]]); }
					catch (NumberFormatException e)
					{
						if (CPLogger.isInitialized) CPLogger.dataWarning.warn("Line " + line + ": Expected numerical input, instead value was '" + row[valInd[i]] + "'. Using default value instead.");
						val = Double.NaN;
					}
					catch (IndexOutOfBoundsException e)
					{
						val = Double.NaN;
					}

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
