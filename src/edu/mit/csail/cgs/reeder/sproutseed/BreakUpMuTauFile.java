package edu.mit.csail.cgs.reeder.sproutseed;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.TreeSet;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.reeder.sprout.SproutUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class BreakUpMuTauFile {

	/**
	 * @param args
	 * @throws NotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		Genome g = SproutUtils.parseGenome(args);
		int buffer = Args.parseInteger(args, "buffer", 0);
		int numregions = Args.parseInteger(args, "numregions", 1);
		String mutaufile = Args.parseString(args, "mutaufile", "");
		String outbase = Args.parseString(args, "outbase", "");
		String chrom = Args.parseString(args, "chrom", "");
		TreeSet<Point> mu = new TreeSet<Point>();
		BufferedReader r = new BufferedReader(new FileReader(mutaufile));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point tmp = Point.fromString(g, split[0]);
			if (chrom.equals("") || chrom.equals(tmp.getChrom())) {
				mu.add(tmp);
			}
		}
		int perRegion = mu.size() / numregions;
		int filecount = 0;
		int count = 0;
		PrintStream out = new PrintStream(outbase+filecount+".txt");
		int index = 0;
		Point endPoint;
		Point lastPoint = mu.first();
		for (Point mup : mu) {
			if (!mup.getChrom().equals(lastPoint.getChrom())) {
				out.flush();
				out.close();
				filecount++;
				out = new PrintStream(outbase+filecount+".txt");
				out.println(mup);
				lastPoint = mup;
				count = 1;
			} else if (count==perRegion) {
				endPoint = mup;
				out.println(mup);
				count++;
				lastPoint = mup;
			} else if (count>perRegion) {
				if (lastPoint.distance(mup) > buffer) {
					out.flush();
					out.close();
					filecount++;
					out = new PrintStream(outbase+filecount+".txt");
					out.println(mup);
					lastPoint = mup;
					count = 1;
				} else {
					endPoint = mup;
					out.println(mup);
					count++;
					lastPoint = mup;
				}
			} else {
				out.println(mup);
				count++;
				lastPoint = mup;
			}
			index++;
		}
		out.flush();
		out.close();
	}

}
