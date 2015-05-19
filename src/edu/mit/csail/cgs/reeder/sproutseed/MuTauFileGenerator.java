package edu.mit.csail.cgs.reeder.sproutseed;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.reeder.sprout.SproutUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class MuTauFileGenerator {

	/**
	 * @param args
	 * @throws NotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		Genome g = SproutUtils.parseGenome(args);
		Region region = Region.fromString(g, Args.parseString(args, "region", ""));
		int spacing = Args.parseInteger(args, "spacing", 0);
		int buffer = Args.parseInteger(args, "buffer", 0);
		String readfile = Args.parseString(args, "readfile", "");
		String outfile = Args.parseString(args, "outfile", "");
		PrintStream out = new PrintStream(outfile);
		TreeSet<StrandedPoint> points = new TreeSet<StrandedPoint>();
		SortedSet<StrandedPoint> mu = new TreeSet<StrandedPoint>();
		Set<Region> regionset = new HashSet<Region>();
		if (region==null) {
			List<String> chromlist = g.getChromList();
			for (String c : chromlist) {
				regionset.add(SproutUtils.chromRegion(c, g));
			}
		} else {
			regionset.add(region);
		}
		for (Region reg : regionset) {
			points = new TreeSet<StrandedPoint>();
			mu = new TreeSet<StrandedPoint>();
			BufferedReader r = new BufferedReader(new FileReader(readfile));
			String s;
			String[] split;
			while ((s = r.readLine()) != null) {
				split = s.split("\t");
				StrandedPoint p1 = StrandedPoint.fromString(g, split[0]);
				StrandedPoint p2 = StrandedPoint.fromString(g, split[1]);
				if (reg==null || reg.contains(p1)) {
					points.add(p1);
				}
				if (reg==null || reg.contains(p2)) {
					points.add(p2);
				}
			}
			if (!points.isEmpty()) {
				mu.add(points.first());
				while(mu.last().compareTo(points.last())<0) {
					StrandedPoint tmp = new StrandedPoint(g,mu.last().getChrom(),mu.last().getLocation()+spacing,'+');
					Pair<StrandedPoint,StrandedPoint> minmax = minmax(tmp,buffer);
					SortedSet<StrandedPoint> subset = points.subSet(minmax.car(), minmax.cdr());
					if (subset.isEmpty()) {
						StrandedPoint htmp = points.higher(tmp);
						if (htmp==null) {
							break;
						} else {
							mu.add(htmp);
						}
					} else {
						mu.add(tmp);
					}
				}
				for (Point p : mu) {
					out.println(p+"\t"+0);
					//out.println(p+"\t"+1);
					//out.println(p+"\t"+2);
				}
			}
			r.close();
		}
		out.flush();
		out.close();
	}

	public static Pair<StrandedPoint,StrandedPoint> minmax(StrandedPoint p, int buffer) {
		return new Pair<StrandedPoint,StrandedPoint>(new StrandedPoint(p.getGenome(),p.getChrom(),p.getLocation()-buffer,p.getStrand()),
				new StrandedPoint(p.getGenome(),p.getChrom(),p.getLocation()+buffer,p.getStrand()));
	}

}
