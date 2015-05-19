package edu.mit.csail.cgs.reeder.sproutseed;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.reeder.sprout.SproutUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class ConsolidateMuFileStage1Results {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws NotFoundException 
	 */
	public static void main(String[] args) throws IOException, NotFoundException {
		Genome g = SproutUtils.parseGenome(args);
		String filebase = Args.parseString(args, "filebase", "");
		String outbase = Args.parseString(args, "outbase", "");
		int numfiles = Args.parseInteger(args, "numfiles", 0);
		String readfile = Args.parseString(args, "readfile", "");
		int numreads = Args.parseInteger(args, "numreads", 0);
		Map<String,Map<Point,Float>> chromPoints = new HashMap<String,Map<Point,Float>>();
		Map<String,Float> chromNoise = new HashMap<String,Float>();
		SproutStorage storage = SproutStorage.fromFile(g, readfile, numreads);
		for (int i=0; i<numfiles; i++) {
			try {
				List<Point> pointlist = new ArrayList<Point>();
				BufferedReader r = new BufferedReader(new FileReader(filebase+i+".txt"));
				String s;
				String[] split;
				while (!(s = r.readLine()).equals("")) {
					split = s.split("\t");
					pointlist.add(Point.fromString(g, split[0]));
				}
				Region tmpreg = pointlist.get(0).expand(2000).combine(pointlist.get(pointlist.size()-1).expand(2000));
				Pair<Integer,Integer> minmaxn = storage.getLeftMinMaxn(tmpreg);
				float numpoints = (float)(minmaxn.cdr()-minmaxn.car());
				r.readLine();//rho
				s = r.readLine();//pi
				split = s.split("\t");
				if (i==0) {
					System.err.println(split.length+"\t"+pointlist.size());
				}
				String chrom = pointlist.get(0).getChrom();
				if (!chromPoints.containsKey(chrom)) {
					chromPoints.put(chrom, new HashMap<Point,Float>());
				}
				Map<Point,Float> map = chromPoints.get(chrom);
				for (int j=0; j<split.length-1; j++) {
					float tmp = Float.valueOf(split[j]);
					if (tmp>0) {
						Point tmpp = pointlist.get(j);
						
						map.put(tmpp, tmp*numpoints);
					}
				}
				if (!chromNoise.containsKey(chrom)) {
					chromNoise.put(chrom, 0f);
				}
				chromNoise.put(chrom, chromNoise.get(chrom)+(Float.valueOf(split[split.length-1])*numpoints));
				r.close();
				System.err.println("done with "+i);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		for (String c : chromPoints.keySet()) {
			Map<Point,Float> map = chromPoints.get(c);
			float norm = chromNoise.get(c);
			for (Point p : map.keySet()) {
				norm += map.get(p);
			}
			PrintStream out = new PrintStream(outbase+c+".txt");
			for (Point p : map.keySet()) {
				out.println(p+"\t"+0);
			}
			out.println();
			out.println("0.7");
			for (Point p : map.keySet()) {
				out.print((map.get(p)/norm)+"\t");
			}
			out.print((chromNoise.get(c)/norm));
			out.println();
			float chival = 1f / 3f;
			out.println(chival+"\t"+chival+"\t"+chival);
			out.flush();
			out.close();
		}
	}

}
