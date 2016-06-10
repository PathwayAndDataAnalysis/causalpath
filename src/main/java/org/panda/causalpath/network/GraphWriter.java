package org.panda.causalpath.network;

import org.panda.causalpath.data.*;
import org.panda.causalpath.resource.MutationClassifier;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Created by babur on 4/6/16.
 */
public class GraphWriter
{
	private Color activatingBorderColor;
	private Color inhibitingBorderColor;
	private boolean useGeneBGForTotalProtein;

	/**
	 * An object that can produce a color for a given up/downregulation value.
	 */
	private ValToColor vtc;

	private Set<RelationAndSelectedData> relations;

	public GraphWriter(Set<RelationAndSelectedData> relations)
	{
		this.relations = relations;
		activatingBorderColor = new Color(0, 180, 20);
		inhibitingBorderColor = new Color(180, 0, 20);

		vtc = new ValToColor(new double[]{-1, 0, 1},
			new Color[]{new Color(40, 80, 255), Color.WHITE, new Color(255, 80, 40)});

		useGeneBGForTotalProtein = false;
	}

	public void setUseGeneBGForTotalProtein(boolean useGeneBGForTotalProtein)
	{
		this.useGeneBGForTotalProtein = useGeneBGForTotalProtein;
	}

	public void setActivatingBorderColor(Color activatingBorderColor)
	{
		this.activatingBorderColor = activatingBorderColor;
	}

	public void setInhibitingBorderColor(Color inhibitingBorderColor)
	{
		this.inhibitingBorderColor = inhibitingBorderColor;
	}

	public void setExpColorSchema(ValToColor vtc)
	{
		this.vtc = vtc;
	}

	public void writeGeneCentric(String filename) throws IOException
	{
		if (!filename.endsWith(".sif")) filename += ".sif";

		// write relations
		BufferedWriter writer1 = new BufferedWriter(new FileWriter(filename));
		relations.forEach(r -> FileUtil.writeln(r.relation.toString(), writer1));
		writer1.close();

		Set<String> totalProtUsedUp = new HashSet<>();

		filename = filename.substring(0, filename.lastIndexOf(".")) + ".format";
		BufferedWriter writer2 = new BufferedWriter(new FileWriter(filename));
		writer2.write("node\tall-nodes\tcolor\t255 255 255\n");
		relations.stream().map(r -> new ExperimentData[]{r.source, r.target}).flatMap(Arrays::stream).distinct()
			.forEach(data ->
		{
			String colS = "255 255 255";
			String bor = "0 0 0";
			String let = "x";

			if (data.hasChangeDetector())
			{
				if (data.getChangeSign() != 0) colS = vtc.getColorInString(data.getChangeValue());
			}

			if (data instanceof PhosphoProteinData)
			{
				PhosphoProteinData pd = (PhosphoProteinData) data;
				if (pd.getEffect() > 0) bor = inString(activatingBorderColor);
				else if (pd.getEffect() < 0) bor = inString(inhibitingBorderColor);

				let = "p";
			}
			else if (data instanceof ProteinData)
			{
				let = "t";
			}
			else if (data instanceof MutationData)
			{
				let = "m";
				if (data.getGeneSymbols().stream().anyMatch(MutationClassifier::isActivatedByMutation))
					bor = inString(activatingBorderColor);
				else bor = inString(inhibitingBorderColor);
			}
			else if (data instanceof CNAData)
			{
				let = "c";
			}
			else if (data instanceof ExpressionData)
			{
				let = "r";
			}
			else if (data instanceof ActivityData)
			{
				let = "a";
				bor = inString(activatingBorderColor);
			}

			String tooltip = data.id;
			if (data.hasChangeDetector()) tooltip += ", " + data.getChangeValue();

			for (String gene : data.getGeneSymbols())
			{
				if (useGeneBGForTotalProtein && let.equals("t") && !totalProtUsedUp.contains(gene))
				{
					FileUtil.writeln("node\t" + gene + "\tcolor\t" + colS, writer2);
					FileUtil.writeln("node\t" + gene + "\ttooltip\t" + tooltip, writer2);
					totalProtUsedUp.add(gene);
				}
				else
				{
					FileUtil.writeln("node\t" + gene + "\trppasite\t" + tooltip + "|" + let + "|" + colS + "|" + bor,
						writer2);
				}
			}
		});

		writer2.close();
	}

	private String inString(Color c)
	{
		return c.getRed() + " " + c.getGreen() + " " + c.getBlue();
	}

	public void writeDataCentric(String filename) throws IOException
	{
		if (!filename.endsWith(".sif")) filename += ".sif";

		// write relations
		BufferedWriter writer1 = new BufferedWriter(new FileWriter(filename));
		relations.forEach(r -> FileUtil.writeln(r.source.id + "\t" + r.relation.type.name + "\t" + r.target.id, writer1));
		writer1.close();

		filename = filename.substring(0, filename.lastIndexOf(".")) + ".format";
		BufferedWriter writer2 = new BufferedWriter(new FileWriter(filename));
		writer2.write("node\tall-nodes\tcolor\t255 255 255\n");
		relations.stream().map(r -> new ExperimentData[]{r.source, r.target}).flatMap(Arrays::stream).distinct()
			.filter(data -> data.hasChangeDetector() && data.getChangeSign() != 0).forEach(data ->
				FileUtil.writeln("node\t" + data.id + "\tcolor\t" + vtc.getColorInString(data.getChangeValue()),
					writer2));

		writer2.close();
	}
}
