package edu.mit.csail.cgs.reeder.sproutseed;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Array;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.reeder.sprout.SproutUtils;
//import edu.mit.csail.cgs.shaun.deepseq.utilities.EMStepPlotter;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class SproutSeed {

	//private static final boolean DUMP = true;
	private static final int MINMAXCOUNT = 20;
	private static final int NO_ALPHA_ITER = 100;
	private static final int ALPHA_ANNEAL_ITER = 150;
	private static final int X_INCREMENT = 100;
	private static final int HISTORY = 40;
	private static final float RHO_ZONE = 0f;
	private static final float LIKELIHOOD_ZONE = 0f;
	private DateFormat dfm = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	private Genome g;
	private float alpha, beta, a, b;
	private Map<Integer,float[]>[] gamma;
	//private float[][][] gamma; // PET index x event index x (event index or -)
	private float[][] gamma2; // event index x event type
	private Map<Pair<Integer,Integer>,Float>[] gamma3;
	private float[] chi; // component weights for event types (stage 2)
	private float[] oldchi;
	private float rho;
	private float[] pi; // event index
	private Map<Pair<Integer,Integer>,Float> psi; // event index x event index
	private int[] tau; // event index -> event type
	private Map<Point,Integer> indexmap;
	private SortedSet<Point> muset;
	private Point[] mu; // event index -> genomic coordinate
	private Point[] oldmu;
	private Point[] originalmu;
	private double[] ll;
	private Map<Pair<Integer,Integer>,Double> interll;
	//private PairedReadDistribution[][] distros; // event index x (event index or -)
	//private PairedReadDistribution[] distTypes; // TSS+, TSS-, and nonTSS
	private PairedReadDistribution genericDist;
	private PairedReadDistribution dist;
	//private Pair<StrandedPoint,StrandedPoint>[] X; //the PETs
	private int nullIndex; // index associated with -
	private SproutStorage storage;
	private Set<Integer> deadComponents = new HashSet<Integer>();
	private int maxIters;
	private String dumpfile;
	private String imagename;
	private Region region;
	private double totalLikelihood;
	private double[] likelihoodHistory = new double[HISTORY];
	private float[] rhoHistory = new float[HISTORY];
	private List<Double> likelihoodList = new ArrayList<Double>();
	private boolean DUMP = false;

	/**
	 * @param args
	 * @throws IOException 
	 * @throws NotFoundException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		SproutSeed seed = new SproutSeed();
		int stage = Args.parseInteger(args, "stage", 0);
		if (stage==1) {
			seed.parseArgs(args);
			System.err.println("initialization complete "+seed.dfm.format(new Date()));
			seed.runStage1(-1);
			seed.dumpStage1(Args.parseString(args, "outfile", ""));
		} else if (stage==2) {
			seed.parseArgsStage2(args);
			System.err.println("initialization complete "+seed.dfm.format(new Date()));
			seed.runStage2(seed.maxIters);
			seed.dumpStage2(Args.parseString(args, "outfile", ""));
		} else if (stage==3) {
			seed.parseArgsStage3(args);
			System.err.println("initialization complete "+seed.dfm.format(new Date()));
			seed.runStage3();
			seed.dumpStage3(Args.parseString(args, "outfile", ""));
			seed.eventReport(Args.parseString(args, "eventout", ""));
			seed.interactionReport(Args.parseString(args, "interactionout", ""));
		} else if (stage==4) {
			seed.parseArgsStage2Hybrid(args);
			System.err.println("initialization complete "+seed.dfm.format(new Date()));
			seed.runStage2Hybrid(Args.parseInteger(args, "numcycles", 10));
			seed.dumpStage2(Args.parseString(args, "outfile", ""));
			seed.stage2HybridEventReport(Args.parseString(args, "eventout", ""));
		} else if (stage==5) {
			seed.parseArgsll(args);
			System.err.println("initialization complete "+seed.dfm.format(new Date()));
			seed.runInterLL();
			seed.llReport(Args.parseString(args, "llout", ""));
		}
	}

	public void dumpStage3(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		for (int i=0; i<mu.length; i++) {
			out.println(mu[i]);
		}
		out.println();
		out.println(rho);
		for (int i=0; i<pi.length; i++) {
			out.print(pi[i]+"\t");
		}
		out.println();
		for (int i=0; i<chi.length; i++) {
			out.print(chi[i]+"\t");
		}
		out.println();
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			out.println(comp.car()+"\t"+comp.cdr()+"\t"+psi.get(comp));
		}
		out.flush();
		out.close();
	}

	public void dumpStage2(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		for (int i=0; i<mu.length; i++) {
			out.println(mu[i]);
		}
		out.println();
		out.println(rho);
		for (int i=0; i<pi.length; i++) {
			out.print(pi[i]+"\t");
		}
		out.println();
		for (int i=0; i<chi.length; i++) {
			out.print(chi[i]+"\t");
		}
		out.println();
		out.flush();
		out.close();
	}

	public void stage2HybridEventReport(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		float[] count = new float[pi.length];
		for (int j=0; j<mu.length; j++) {
			Point currentmu = mu[j];
			Region currentmuregion = currentmu.expand(dist.radius);
			Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(currentmuregion);
			int leftminn = leftminmaxn.car();
			int leftmaxn = leftminmaxn.cdr();
			for (int n = leftminn; n < leftmaxn; n++) {
				if (gamma[n].containsKey(j)) {
					count[j] += gamma[n].get(j)[0];
				}
			}
		}
		float pinorm = 0;
		for (int j=0; j<mu.length; j++) {
			count[j] = Math.max(count[j]-alpha, 0f);
			pinorm += count[j];
		}
		count[count.length-1]= 0; 
		for (int j=0; j<storage.getLeftN(); j++) {
			if (gamma[j].containsKey(nullIndex)) {
				count[count.length-1] += gamma[j].get(nullIndex)[0];
			}
		}
		pinorm += count[count.length-1];
		out.println("Attachment\tPETs (out of "+pinorm+")\tPositive TSS\tNegative TSS\tNon-TSS");
		for (int i=0; i<mu.length; i++) {
			out.print(mu[i]+"\t"+count[i]);
			for (int j=0; j<chi.length; j++) {
				out.print("\t"+gamma2[i][j]);
			}
			out.print("\t"+originalmu[i]+"\t"+originalmu[i].distance(mu[i]));
			out.println();
		}
		out.println("noise\t"+count[count.length-1]);
		out.flush();
		out.close();
	}

	public void dumpStage1(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		for (int i=0; i<mu.length; i++) {
			out.println(mu[i]+"\t"+tau[i]);
		}
		out.println();
		out.println(rho);
		for (int i=0; i<pi.length; i++) {
			out.print(pi[i]+"\t");
		}
		out.println();
		out.flush();
		out.close();
	}

	public void dumpMu(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		for (int i=0; i<mu.length; i++) {
			out.println(mu[i]+"\t"+tau[i]);
		}
		out.println();
		out.flush();
		out.close();
	}

	public void dumpMuStage2(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		for (int i=0; i<mu.length; i++) {
			out.println(mu[i]);
		}
		out.println();
		out.flush();
		out.close();
	}

	public void dumpRho(String file ) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		out.println(rho);
		out.flush();
		out.close();
	}

	public void dumpPi(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		for (int i=0; i<pi.length; i++) {
			out.print(pi[i]+"\t");
		}
		out.println();
		out.flush();
		out.close();
	}

	public void dumpLikelihood(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		out.println(totalLikelihood);
		out.flush();
		out.close();
	}

	public void dumpChi(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		for (int i=0; i<chi.length; i++) {
			out.print(chi[i]+"\t");
		}
		out.println();
		out.flush();
		out.close();
	}

	private void dumpGamma2(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		for (int i=0; i<gamma2.length; i++) {
			for (int j=0; j<gamma2[i].length; j++) {
				out.print(gamma2[i][j]+"\t");
			}
			out.println();
		}
		out.println();
		out.flush();
		out.close();
	}


	private void dumpGamma(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		for (int i=0; i<gamma.length; i++) {
			for (int j : gamma[i].keySet()) {
				float[] tmp = gamma[i].get(j);
				out.print(tmp[0]+"\t");
			}
			out.println();
		}
		out.println();
		out.flush();
		out.close();
	}


	public void dump(String file) throws FileNotFoundException {
		FileOutputStream fout = new FileOutputStream(file,true);
		PrintStream out = new PrintStream(fout);
		for (int i=0; i<mu.length; i++) {
			out.println(mu[i]+"\t"+tau[i]);
		}
		out.println();
		out.println(rho);
		for (int i=0; i<pi.length; i++) {
			out.print(pi[i]+"\t");
		}
		out.println();
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			out.println(comp.car()+"\t"+comp.cdr()+"\t"+psi.get(comp));
		}
		out.flush();
		out.close();
	}

	private void dumpStage1() throws FileNotFoundException {
		if (DUMP) {
			dumpMu(dumpfile+"mu.txt");
			dumpRho(dumpfile+"rho.txt");
			dumpPi(dumpfile+"pi.txt");
		}
		dumpLikelihood(dumpfile+"likelihood.txt");
	}

	private void dumpStage2() throws FileNotFoundException {
		if (DUMP) {
			dumpMuStage2(dumpfile+"mu.txt");
			dumpChi(dumpfile+"chi.txt");
			dumpGamma2(dumpfile+"gamma2.txt");
		}
		dumpLikelihood(dumpfile+"likelihood.txt");
	}

	private void dumpStage3() throws FileNotFoundException {
		dumpLikelihood(dumpfile+"likelihood.txt");
	}

	private void dump() throws FileNotFoundException {
		dump(dumpfile);
	}

	public void parseArgs(String[] args) throws NotFoundException, IOException {
		maxIters = Args.parseInteger(args, "maxiters", 1000);
		g = SproutUtils.parseGenome(args);
		alpha = Args.parseFloat(args, "alpha", 0);
		beta = Args.parseFloat(args, "beta", 0);
		a = Args.parseFloat(args, "a", 0);
		b = Args.parseFloat(args, "b", 0);
		region = Region.fromString(g, Args.parseString(args, "region", ""));
		String mutaufile = Args.parseString(args, "mutaufile", "");
		if (region==null) {
			initializeMuTauAndRegion(mutaufile);
		} else {
			initializeMuTau(mutaufile);
		}
		String readfile = Args.parseString(args, "readfile", "");
		if (region==null) {
			storage = SproutStorage.fromFile(g, readfile);
		} else {
			storage = SproutStorage.fromFile(g, readfile, region);
		}
		String directory = Args.parseString(args, "directory", "");
		dist = new PairedReadDistribution(0);
		dist.initializeFromDirectoryStage1(directory);
		//gamma = new float[storage.getLeftN()][mu.length+1][2];
		gamma = (Map<Integer,float[]>[]) Array.newInstance(Map.class, storage.getLeftN());
		initializeHashGamma();
		rho = Args.parseFloat(args, "rho", 0);
		pi = new float[mu.length+1];
		for (int i=0; i<pi.length; i++) {
			pi[i] = 1f / ((float)pi.length);
		}
		/*
		psi = new float[mu.length+1][mu.length+1];
		float fmu = (float)(mu.length+1);
		float psiinit = 1f / (fmu*(fmu+1f)/2f);
		for (int i=0; i<psi.length; i++) {
			for (int j=i; j<psi[i].length; j++) {
				psi[i][j] = psiinit;
			}
		}
		 */
		if (Args.parseFlags(args).contains("dump")) DUMP = true;
		dumpfile = Args.parseString(args, "dumpfile", "");
		imagename = Args.parseString(args, "imagename", "");
	}

	public void parseArgsStage2Hybrid(String[] args) throws NotFoundException, IOException {
		maxIters = Args.parseInteger(args, "maxiters", 1000);
		g = SproutUtils.parseGenome(args);
		alpha = Args.parseFloat(args, "alpha", 0);
		beta = Args.parseFloat(args, "beta", 0);
		a = Args.parseFloat(args, "a", 0);
		b = Args.parseFloat(args, "b", 0);
		region = Region.fromString(g, Args.parseString(args, "region", ""));
		String mutaufile = Args.parseString(args, "mutaufile", "");
		if (region==null) {
			initializeMuTauAndRegion(mutaufile);
		} else {
			initializeMuTau(mutaufile);
		}
		String readfile = Args.parseString(args, "readfile", "");
		if (region==null) {
			storage = SproutStorage.fromFile(g, readfile);
		} else {
			storage = SproutStorage.fromFile(g, readfile, region);
		}
		String directory = Args.parseString(args, "directory", "");
		dist = new PairedReadDistribution(0);
		dist.initializeFromDirectory(directory);
		//gamma = new float[storage.getLeftN()][mu.length+1][2];
		gamma = (Map<Integer,float[]>[]) Array.newInstance(Map.class, storage.getLeftN());

		String stage1file = Args.parseString(args, "stage1file", "");
		parseStage1ResultsRegion(stage1file);
		initializeHashGamma();
		gamma2 = new float[mu.length][3];
		chi = new float[3];
		oldchi = new float[3];
		float chiInit = 1f / ((float)chi.length);
		for (int i=0; i<chi.length; i++) {
			chi[i] = chiInit;
			oldchi[i] = chiInit;
		}
		for (int i=0; i<gamma2.length; i++) {
			for (int j=0; j<gamma2[i].length; j++) {
				gamma2[i][j] = chiInit;
			}
		}
		if (Args.parseFlags(args).contains("dump")) DUMP = true;
		dumpfile = Args.parseString(args, "dumpfile", "");
		imagename = Args.parseString(args, "imagename", "");
	}

	private void initializeHashGamma() {
		muset = new TreeSet<Point>();
		indexmap = new HashMap<Point,Integer>();
		for (int i=0; i<mu.length; i++) {
			muset.add(mu[i]);
			indexmap.put(mu[i], i);
		}
		for (int i=0; i<storage.getLeftN(); i++) {
			gamma[i] = new HashMap<Integer,float[]>();
			Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
			Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius-1);
			Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius+1);
			SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
			for (Point left : leftset) {
				float[] tmparr = {0f, 0f};
				gamma[i].put(indexmap.get(left),tmparr);
			}
			gamma[i].put(nullIndex, new float[1]);
		}
	}

	private void updateHashGamma(Point newmu, int newindex) {
		Region newmuregion = newmu.expand(dist.radius);
		Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(newmuregion);
		int leftminn = leftminmaxn.car();
		int leftmaxn = leftminmaxn.cdr();
		for (int i=leftminn; i<leftmaxn; i++) {
			if (!gamma[i].containsKey(newindex)) {
				gamma[i].put(newindex, new float[2]);
			}
		}
	}


	public void parseArgsStage3(String[] args) throws NotFoundException, IOException {
		maxIters = Args.parseInteger(args, "maxiters", 1000);
		g = SproutUtils.parseGenome(args);
		alpha = Args.parseFloat(args, "alpha", 0);
		beta = Args.parseFloat(args, "beta", 0);
		a = Args.parseFloat(args, "a", 0);
		b = Args.parseFloat(args, "b", 0);
		boolean hybrid = Args.parseFlags(args).contains("hybrid");
		String stage2file = Args.parseString(args, "stage2file", "");
		String gamma2file = Args.parseString(args, "gamma2file", "");
		parseStage2Results(stage2file);
		String readfile = Args.parseString(args, "readfile", "");
		region = Region.fromString(g, Args.parseString(args, "region", ""));
		if (region==null) {
			storage = SproutStorage.fromFile(g, readfile);
		} else {
			storage = SproutStorage.fromFile(g, readfile, region);
		}
		String genericdirectory = Args.parseString(args, "genericdirectory", "");
		genericDist = new PairedReadDistribution(0);
		genericDist.initializeFromDirectoryStage1(genericdirectory);
		String directory = Args.parseString(args, "directory", "");
		if (!hybrid && directory.equals(genericdirectory)) {
			dist = new PairedReadDistribution(0);
			dist.initializeFromDirectoryStage1(directory);
		} else {
			dist = new PairedReadDistribution(0);
			dist.initializeFromDirectory(directory);
		}
		gamma = (Map<Integer,float[]>[]) Array.newInstance(Map.class, storage.getLeftN());
		if (hybrid) {
			computeGammaForStage3Hybrid();
		} else {
			computeGammaForStage3();
		}
		if (gamma2file.equals("")) {
			gamma2 = new float[mu.length][1];
			for (int i=0; i<mu.length; i++) {
				gamma2[i][0] = 1.0f;
			}
		} else {
			gamma2 = new float[mu.length][3];
			if (gamma2file.equals("compute")) {
				computeGamma2ForStage3();
			} else {
				parseGamma2File(gamma2file);
			}
		}
		gamma3 = (Map<Pair<Integer,Integer>,Float>[]) Array.newInstance(Map.class, storage.getLeftN());
		initializeGamma3();
		psi = new HashMap<Pair<Integer,Integer>,Float>();
		for (int i=0; i<gamma3.length; i++) {
			for (Pair<Integer,Integer> comp : gamma3[i].keySet()) {
				psi.put(comp, 0f);
			}
		}
		float psiinit = 1f/((float)psi.size());
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			psi.put(comp, psiinit);
		}
		if (Args.parseFlags(args).contains("dump")) DUMP = true;
		dumpfile = Args.parseString(args, "dumpfile", "");
		dumpGamma(dumpfile+"gamma.txt");
	}

	public void parseArgsll(String[] args) throws NotFoundException, IOException {
		maxIters = Args.parseInteger(args, "maxiters", 1000);
		g = SproutUtils.parseGenome(args);
		alpha = Args.parseFloat(args, "alpha", 0);
		beta = Args.parseFloat(args, "beta", 0);
		a = Args.parseFloat(args, "a", 0);
		b = Args.parseFloat(args, "b", 0);
		boolean hybrid = Args.parseFlags(args).contains("hybrid");
		String stage3file = Args.parseString(args, "stage3file", "");
		String gamma2file = Args.parseString(args, "gamma2file", "");
		parseStage3Results(stage3file);
		String readfile = Args.parseString(args, "readfile", "");
		region = Region.fromString(g, Args.parseString(args, "region", ""));
		if (region==null) {
			storage = SproutStorage.fromFile(g, readfile);
		} else {
			storage = SproutStorage.fromFile(g, readfile, region);
		}
		String genericdirectory = Args.parseString(args, "genericdirectory", "");
		genericDist = new PairedReadDistribution(0);
		//genericDist.initializeFromDirectoryStage1(genericdirectory);
		String directory = Args.parseString(args, "directory", "");
		if (!hybrid && directory.equals(genericdirectory)) {
			dist = new PairedReadDistribution(0);
			dist.initializeFromDirectoryStage1(directory);
		} else {
			dist = new PairedReadDistribution(0);
			dist.initializeFromDirectory(directory);
		}
		gamma = (Map<Integer,float[]>[]) Array.newInstance(Map.class, storage.getLeftN());
		if (hybrid) {
			computeGammaForStage3Hybrid();
		} else {
			computeGammaForStage3Hybrid();
		}
		if (gamma2file.equals("")) {
			gamma2 = new float[mu.length][1];
			for (int i=0; i<mu.length; i++) {
				gamma2[i][0] = 1.0f;
			}
		} else {
			gamma2 = new float[mu.length][3];
			if (gamma2file.equals("compute")) {
				computeGamma2ForStage3();
			} else {
				parseGamma2File(gamma2file);
			}
		}
		gamma3 = (Map<Pair<Integer,Integer>,Float>[]) Array.newInstance(Map.class, storage.getLeftN());
		initializeGamma3ll();
		if (Args.parseFlags(args).contains("dump")) DUMP = true;
		dumpfile = Args.parseString(args, "dumpfile", "");
		//dumpGamma(dumpfile+"gamma.txt");
		interll = new HashMap<Pair<Integer,Integer>,Double>();
	}


	private void initializeGamma3() {
		muset = new TreeSet<Point>();
		indexmap = new HashMap<Point,Integer>();
		for (int i=0; i<mu.length; i++) {
			muset.add(mu[i]);
			indexmap.put(mu[i], i);
		}
		for (int i=0; i<storage.getLeftN(); i++) {
			gamma3[i] = new HashMap<Pair<Integer,Integer>,Float>();
			Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
			Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius);
			Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius);
			SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
			tmpleft = new Point(g,tmp.cdr().getChrom(), tmp.cdr().getLocation()-dist.radius);
			tmpright = new Point(g, tmp.cdr().getChrom(), tmp.cdr().getLocation()+dist.radius);
			SortedSet<Point> rightset = muset.subSet(tmpleft, tmpright);
			for (Point left : leftset) {
				for (Point right : rightset) {
					gamma3[i].put(new Pair<Integer,Integer>(indexmap.get(left),indexmap.get(right)),0f);
					gamma3[i].put(new Pair<Integer,Integer>(nullIndex,indexmap.get(right)), 0f);
				}
				gamma3[i].put(new Pair<Integer,Integer>(indexmap.get(left),nullIndex), 0f);
			}
			gamma3[i].put(new Pair<Integer,Integer>(nullIndex,nullIndex), 0f);
		}
	}

	private void initializeGamma3ll() {
		muset = new TreeSet<Point>();
		indexmap = new HashMap<Point,Integer>();
		for (int i=0; i<mu.length; i++) {
			muset.add(mu[i]);
			indexmap.put(mu[i], i);
		}
		float rsqr = (float)(region.getWidth()*region.getWidth());
		//float minval = 10f*Float.MIN_VALUE*rsqr / (1 - rho - 20f * nullIndex * Float.MIN_VALUE * rsqr);
		float minval = 2.5E-28f;
		System.err.println("minval "+minval);
		System.err.println("Float.MIN_VALUE "+Float.MIN_VALUE);
		Pair<Integer,Integer> comp;
		for (int i=0; i<storage.getLeftN(); i++) {
			gamma3[i] = new HashMap<Pair<Integer,Integer>,Float>();
			Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
			Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius);
			Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius);
			SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
			tmpleft = new Point(g,tmp.cdr().getChrom(), tmp.cdr().getLocation()-dist.radius);
			tmpright = new Point(g, tmp.cdr().getChrom(), tmp.cdr().getLocation()+dist.radius);
			SortedSet<Point> rightset = muset.subSet(tmpleft, tmpright);
			for (Point left : leftset) {
				for (Point right : rightset) {
					comp = new Pair<Integer,Integer>(indexmap.get(left),indexmap.get(right));
					if (psi.containsKey(comp)) {
						gamma3[i].put(comp,0f);
					}
					comp = new Pair<Integer,Integer>(nullIndex,indexmap.get(right));
					if (psi.containsKey(comp)) {
						gamma3[i].put(comp, 0f);
						psi.put(comp, minval);
					} else {
						gamma3[i].put(comp, 0f);
						psi.put(comp, minval);
					}
				}
				comp = new Pair<Integer,Integer>(indexmap.get(left),nullIndex);
				if (psi.containsKey(comp)) {
					gamma3[i].put(comp, 0f);
					psi.put(comp, minval);
				} else {
					gamma3[i].put(comp, 0f);
					psi.put(comp, minval);
				}
			}
			comp = new Pair<Integer,Integer>(nullIndex,nullIndex);
			if (psi.containsKey(comp)) {
				gamma3[i].put(comp, 0f);
				psi.put(comp, minval);
			} else {
				gamma3[i].put(comp, 0f);
				psi.put(comp, minval);
			}
		}
		float sum = 0f;
		for (Pair<Integer,Integer> tmpcomp : psi.keySet()) {
			sum += psi.get(tmpcomp);
		}
		for (Pair<Integer,Integer> tmpcomp : psi.keySet()) {
			psi.put(tmpcomp, psi.get(tmpcomp) / sum);
		}
	}


	public void parseArgsStage2(String[] args) throws NotFoundException, IOException {
		maxIters = Args.parseInteger(args, "maxiters", 1000);
		g = SproutUtils.parseGenome(args);
		alpha = Args.parseFloat(args, "alpha", 0);
		beta = Args.parseFloat(args, "beta", 0);
		a = Args.parseFloat(args, "a", 0);
		b = Args.parseFloat(args, "b", 0);
		String stage1file = Args.parseString(args, "stage1file", "");
		parseStage1Results(stage1file);
		String readfile = Args.parseString(args, "readfile", "");
		region = Region.fromString(g, Args.parseString(args, "region", ""));
		if (region==null) {
			storage = SproutStorage.fromFile(g, readfile);
		} else {
			storage = SproutStorage.fromFile(g, readfile, region);
		}
		String genericdirectory = Args.parseString(args, "genericdirectory", "");
		genericDist = new PairedReadDistribution(0);
		genericDist.initializeFromDirectoryStage1(genericdirectory);
		String directory = Args.parseString(args, "directory", "");
		dist = new PairedReadDistribution(0);
		dist.initializeFromDirectory(directory);
		gamma = (Map<Integer,float[]>[]) Array.newInstance(Map.class, storage.getLeftN());
		computeStage1Gamma();
		gamma2 = new float[mu.length][3];
		if (Args.parseFlags(args).contains("dump")) DUMP = true;
		dumpfile = Args.parseString(args, "dumpfile", "");
		dumpGamma(dumpfile+"gamma.txt");
		chi = new float[3];
		oldchi = new float[3];
		float chiInit = 1f / ((float)chi.length);
		for (int i=0; i<chi.length; i++) {
			chi[i] = chiInit;
			oldchi[i] = chiInit;
		}
	}


	public void parseStage2Results(String file) throws IOException {
		List<Point> mu = new ArrayList<Point>();
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while (!(s = r.readLine()).equals("")) {
			split = s.split("\t");
			mu.add(Point.fromString(g, split[0]));
		}
		rho = Float.valueOf(r.readLine());
		System.err.println("rho: "+rho);
		s = r.readLine();
		split = s.split("\t");
		System.err.println(split.length);
		List<Integer> nonzero = new ArrayList<Integer>();
		for (int i=0; i<split.length-1; i++) {
			float tmp = Float.valueOf(split[i]);
			if (tmp>0) {
				nonzero.add(i);
			}
		}
		this.mu = new Point[nonzero.size()];
		nullIndex = this.mu.length;
		this.oldmu = new Point[nonzero.size()];
		this.pi = new float[nonzero.size()+1];
		for (int i=0; i<nonzero.size(); i++) {
			this.mu[i] = mu.get(nonzero.get(i));
			this.oldmu[i] = mu.get(nonzero.get(i));
			this.pi[i] = Float.valueOf(split[nonzero.get(i)]);
		}
		this.pi[this.pi.length-1] = Float.valueOf(split[split.length-1]);
		float sum = 0;
		for (int i=0; i<this.pi.length; i++) {
			sum += this.pi[i];
		}
		System.err.println("number of components: "+this.mu.length);
		System.err.println("stage 1 pi sum: "+sum);
		this.chi = new float[3];
		this.oldchi = new float[3];
		s = r.readLine();
		if (s==null) {
			this.chi = new float[1];
			this.oldchi = new float[1];
			this.chi[0] = 1.0f;
			this.chi[0] = 1.0f;
		} else {
			split = s.split("\t");
			for (int i=0; i<chi.length; i++) {
				chi[i] = Float.valueOf(split[i]);
				oldchi[i] = chi[i];
			}
		}
		r.close();
	}

	public void parseStage3Results(String file) throws IOException {
		List<Point> mu = new ArrayList<Point>();
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while (!(s = r.readLine()).equals("")) {
			split = s.split("\t");
			mu.add(Point.fromString(g, split[0]));
		}
		rho = Float.valueOf(r.readLine());
		System.err.println("rho: "+rho);
		s = r.readLine();
		split = s.split("\t");
		System.err.println(split.length);
		List<Integer> nonzero = new ArrayList<Integer>();
		for (int i=0; i<split.length-1; i++) {
			float tmp = Float.valueOf(split[i]);
			if (tmp>=0) { //accepting zero values
				nonzero.add(i);
			}
		}
		this.mu = new Point[nonzero.size()];
		nullIndex = this.mu.length;
		this.oldmu = new Point[nonzero.size()];
		this.pi = new float[nonzero.size()+1];
		for (int i=0; i<nonzero.size(); i++) {
			this.mu[i] = mu.get(nonzero.get(i));
			this.oldmu[i] = mu.get(nonzero.get(i));
			this.pi[i] = Float.valueOf(split[nonzero.get(i)]);
		}
		this.pi[this.pi.length-1] = Float.valueOf(split[split.length-1]);
		float sum = 0;
		for (int i=0; i<this.pi.length; i++) {
			sum += this.pi[i];
		}
		System.err.println("number of components: "+this.mu.length);
		System.err.println("stage 1 pi sum: "+sum);
		this.chi = new float[3];
		this.oldchi = new float[3];
		s = r.readLine();
		if (s==null) {
			this.chi = new float[1];
			this.oldchi = new float[1];
			this.chi[0] = 1.0f;
			this.chi[0] = 1.0f;
		} else {
			split = s.split("\t");
			this.chi = new float[split.length];
			this.oldchi = new float[split.length];
			for (int i=0; i<split.length; i++) {
				chi[i] = Float.valueOf(split[i]);
				oldchi[i] = chi[i];
			}
		}
		this.psi = new HashMap<Pair<Integer,Integer>,Float>();
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			if (Float.valueOf(split[2])>0) {
				this.psi.put(new Pair<Integer,Integer>(Integer.valueOf(split[0]),Integer.valueOf(split[1])), Float.valueOf(split[2]));
			}
		}
		r.close();
	}

	public void parseStage1Results(String file) throws IOException {
		List<Point> mu = new ArrayList<Point>();
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while (!(s = r.readLine()).equals("")) {
			split = s.split("\t");
			mu.add(Point.fromString(g, split[0]));
		}
		rho = Float.valueOf(r.readLine());
		System.err.println("rho: "+rho);
		s = r.readLine();
		split = s.split("\t");
		System.err.println(split.length);
		List<Integer> nonzero = new ArrayList<Integer>();
		for (int i=0; i<split.length-1; i++) {
			float tmp = Float.valueOf(split[i]);
			if (tmp>0) {
				nonzero.add(i);
			}
		}
		this.mu = new Point[nonzero.size()];
		nullIndex = this.mu.length;
		this.oldmu = new Point[nonzero.size()];
		this.originalmu = new Point[nonzero.size()];
		this.pi = new float[nonzero.size()+1];
		for (int i=0; i<nonzero.size(); i++) {
			this.mu[i] = mu.get(nonzero.get(i));
			this.oldmu[i] = mu.get(nonzero.get(i));
			this.originalmu[i] = mu.get(nonzero.get(i));
			this.pi[i] = Float.valueOf(split[nonzero.get(i)]);
		}
		this.pi[this.pi.length-1] = Float.valueOf(split[split.length-1]);
		float sum = 0;
		for (int i=0; i<this.pi.length; i++) {
			sum += this.pi[i];
		}
		r.close();
		System.err.println("number of components: "+this.mu.length);
		System.err.println("stage 1 pi sum: "+sum);
	}

	public void parseStage1ResultsRegion(String file) throws IOException {
		List<Integer> contained = new ArrayList<Integer>();
		List<Point> mu = new ArrayList<Point>();
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while (!(s = r.readLine()).equals("")) {
			split = s.split("\t");
			Point tmp = Point.fromString(g, split[0]);
			mu.add(tmp);
			if (region.contains(tmp)) {
				contained.add(mu.size()-1);
			}
		}
		rho = Float.valueOf(r.readLine());
		System.err.println("rho: "+rho);
		s = r.readLine();
		split = s.split("\t");
		System.err.println(split.length);
		List<Integer> nonzero = new ArrayList<Integer>();
		for (int i=0; i<split.length-1; i++) {
			float tmp = Float.valueOf(split[i]);
			if (tmp>0 && contained.contains(i)) {
				nonzero.add(i);
			}
		}
		this.mu = new Point[nonzero.size()];
		nullIndex = this.mu.length;
		this.oldmu = new Point[nonzero.size()];
		this.originalmu = new Point[nonzero.size()];
		this.pi = new float[nonzero.size()+1];
		for (int i=0; i<nonzero.size(); i++) {
			this.mu[i] = mu.get(nonzero.get(i));
			this.oldmu[i] = mu.get(nonzero.get(i));
			this.originalmu[i] = mu.get(nonzero.get(i));
			this.pi[i] = Float.valueOf(split[nonzero.get(i)]);
		}
		this.pi[this.pi.length-1] = Float.valueOf(split[split.length-1]);
		float sum = 0;
		for (int i=0; i<this.pi.length-1; i++) {
			sum += this.pi[i];
		}
		sum /= (1-this.pi[this.pi.length-1]);
		for (int i=0; i<this.pi.length-1; i++) {
			this.pi[i] /= sum;
		}
		r.close();
		System.err.println("number of components: "+this.mu.length);
		System.err.println("stage 1 pi sum: "+sum);
		printPi();
	}

	public void initializeMuTau(String file) throws IOException {
		List<Point> mulist = new ArrayList<Point>();
		List<Integer> taulist = new ArrayList<Integer>();
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point tmp = Point.fromString(g, split[0]);
			if (region.contains(tmp)) {
				mulist.add(tmp);
				taulist.add(Integer.valueOf(split[1]));
			}
		}
		nullIndex = mulist.size();
		mu = new Point[mulist.size()];
		oldmu = new Point[mulist.size()];
		tau = new int[mulist.size()];
		for (int i=0; i<mu.length; i++) {
			mu[i] = mulist.get(i);
			oldmu[i] = mulist.get(i);
			tau[i] = taulist.get(i);
		}
		r.close();
	}

	public void initializeMuTauAndRegion(String file) throws NumberFormatException, IOException {
		List<Point> mulist = new ArrayList<Point>();
		List<Integer> taulist = new ArrayList<Integer>();
		Point minPoint = new Point(g, "Z", Integer.MAX_VALUE);
		Point maxPoint = new Point(g, "", -1);
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point tmp = Point.fromString(g, split[0]);
			mulist.add(tmp);
			if (split.length>1) {
				taulist.add(Integer.valueOf(split[1]));
			} else {
				taulist.add(0);
			}
			if (tmp.compareTo(minPoint)<0) {
				minPoint = tmp;
			}
			if (tmp.compareTo(maxPoint)>0) {
				maxPoint = tmp;
			}
		}
		region = minPoint.expand(2000);
		region = region.combine(maxPoint.expand(2000));
		System.err.println("region: "+region);
		nullIndex = mulist.size();
		mu = new Point[mulist.size()];
		oldmu = new Point[mulist.size()];
		tau = new int[mulist.size()];
		for (int i=0; i<mu.length; i++) {
			mu[i] = mulist.get(i);
			oldmu[i] = mulist.get(i);
			tau[i] = taulist.get(i);
		}
		r.close();
	}


	public void eventReport(String file) throws FileNotFoundException {
		PrintStream out = new PrintStream(file);
		float[] count = new float[pi.length];
		for (int j=0; j<mu.length; j++) {
			Point currentmu = mu[j];
			Region currentmuregion = currentmu.expand(dist.radius);
			Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(currentmuregion);
			int leftminn = leftminmaxn.car();
			int leftmaxn = leftminmaxn.cdr();
			for (int n = leftminn; n < leftmaxn; n++) {
				if (gamma[n].containsKey(j)) {
					count[j] += gamma[n].get(j)[0];
				}
			}
		}
		float pinorm = 0;
		for (int j=0; j<mu.length; j++) {
			count[j] = Math.max(count[j]-alpha, 0f);
			pinorm += count[j];
		}
		count[count.length-1]= 0; 
		for (int j=0; j<storage.getLeftN(); j++) {
			if (gamma[j].containsKey(nullIndex)) {
				count[count.length-1] += gamma[j].get(nullIndex)[0];
			}
		}
		pinorm += count[count.length-1];
		out.println("Attachment\tPETs (out of "+pinorm+")\tPositive TSS\tNegative TSS\tNon-TSS");
		for (int i=0; i<mu.length; i++) {
			out.print(mu[i]+"\t"+count[i]);
			for (int j=0; j<chi.length; j++) {
				out.print("\t"+gamma2[i][j]);
			}
			out.println();
		}
		out.println("noise\t"+count[count.length-1]);
		out.flush();
		out.close();
	}


	public void interactionReport(String file) throws FileNotFoundException {
		PrintStream out = new PrintStream(file);
		Map<Pair<Integer,Integer>,Float> count = new HashMap<Pair<Integer,Integer>,Float>();
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			count.put(comp, 0f);
		}
		for (int i=0; i<storage.getLeftN(); i++) {
			Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
			Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius);
			Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius);
			SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
			tmpleft = new Point(g,tmp.cdr().getChrom(), tmp.cdr().getLocation()-dist.radius);
			tmpright = new Point(g, tmp.cdr().getChrom(), tmp.cdr().getLocation()+dist.radius);
			SortedSet<Point> rightset = muset.subSet(tmpleft, tmpright);
			Pair<Integer,Integer> tmpcomp;
			for (Point left : leftset) {
				for (Point right : rightset) {
					tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),indexmap.get(right));
					count.put(tmpcomp, count.get(tmpcomp)+gamma3[i].get(tmpcomp));
					tmpcomp = new Pair<Integer,Integer>(nullIndex,indexmap.get(right));
					count.put(tmpcomp, count.get(tmpcomp)+gamma3[i].get(tmpcomp));
				}
				tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),nullIndex);
				count.put(tmpcomp, count.get(tmpcomp)+gamma3[i].get(tmpcomp));
			}
			tmpcomp = new Pair<Integer,Integer>(nullIndex,nullIndex);
			count.put(tmpcomp, count.get(tmpcomp)+gamma3[i].get(tmpcomp));
		}
		float psinorm = 0;
		int psialive = 0;
		for (Pair<Integer,Integer> comp : count.keySet()) {
			float tmp = count.get(comp);
			float newtmp = Math.max(0,tmp-beta);
			psinorm += newtmp;
			if (newtmp>0) psialive++;
			count.put(comp, newtmp);

		}
		out.println("Left Anchor\tRight Anchor\tPETs (out of "+psinorm+")\tDistance\tWarpdrive Input");
		for (Pair<Integer,Integer> comp : count.keySet()) {
			float tmpcount = count.get(comp);
			if (tmpcount>0) {
				String left, right;
				if (comp.car().equals(nullIndex)) {
					left = "noise";
				} else {
					left = mu[comp.car()].toString();
				}
				if (comp.cdr().equals(nullIndex)) {
					right = "noise";
				} else {
					right = mu[comp.cdr()].toString();
				}
				String warpstring = "";
				int distance = -1;
				if (!comp.car().equals(nullIndex) && !comp.cdr().equals(nullIndex) && mu[comp.car()].getChrom().equals(mu[comp.cdr()].getChrom())) {
					Region warpregion = mu[comp.car()].expand(1).combine(mu[comp.cdr()].expand(1));
					warpstring = warpregion.toString();
					distance = warpregion.getWidth();
				}
				out.println(left+"\t"+right+"\t"+tmpcount+"\t"+distance+"\t"+warpstring);
			}
		}
		out.flush();
		out.close();
	}

	public void llReport(String file) throws FileNotFoundException {
		PrintStream out = new PrintStream(file);
		Map<Pair<Integer,Integer>,Float> count = new HashMap<Pair<Integer,Integer>,Float>();
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			count.put(comp, 0f);
		}
		for (int i=0; i<storage.getLeftN(); i++) {
			Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
			Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius);
			Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius);
			SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
			tmpleft = new Point(g,tmp.cdr().getChrom(), tmp.cdr().getLocation()-dist.radius);
			tmpright = new Point(g, tmp.cdr().getChrom(), tmp.cdr().getLocation()+dist.radius);
			SortedSet<Point> rightset = muset.subSet(tmpleft, tmpright);
			Pair<Integer,Integer> tmpcomp;
			for (Point left : leftset) {
				for (Point right : rightset) {
					tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),indexmap.get(right));
					if (gamma3[i].containsKey(tmpcomp)) {
						count.put(tmpcomp, count.get(tmpcomp)+gamma3[i].get(tmpcomp));
					}
					tmpcomp = new Pair<Integer,Integer>(nullIndex,indexmap.get(right));
					if (gamma3[i].containsKey(tmpcomp)) {
						count.put(tmpcomp, count.get(tmpcomp)+gamma3[i].get(tmpcomp));
					}
				}
				tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),nullIndex);
				if (gamma3[i].containsKey(tmpcomp)) {
					count.put(tmpcomp, count.get(tmpcomp)+gamma3[i].get(tmpcomp));
				}
			}
			tmpcomp = new Pair<Integer,Integer>(nullIndex,nullIndex);
			if (gamma3[i].containsKey(tmpcomp)) {
				count.put(tmpcomp, count.get(tmpcomp)+gamma3[i].get(tmpcomp));
			}
		}
		float psinorm = 0;
		int psialive = 0;
		for (Pair<Integer,Integer> comp : count.keySet()) {
			float tmp = count.get(comp);
			float newtmp = Math.max(0,tmp-beta);
			psinorm += newtmp;
			if (newtmp>0) psialive++;
			count.put(comp, newtmp);

		}
		out.println("Left Anchor\tRight Anchor\tPETs (out of "+psinorm+")\tDistance\tWarpdrive Input\tLikelihood Difference");
		for (Pair<Integer,Integer> comp : count.keySet()) {
			float tmpcount = count.get(comp);
			if (tmpcount>0) {
				String left, right;
				if (comp.car().equals(nullIndex)) {
					left = "noise";
				} else {
					left = mu[comp.car()].toString();
				}
				if (comp.cdr().equals(nullIndex)) {
					right = "noise";
				} else {
					right = mu[comp.cdr()].toString();
				}
				String warpstring = "";
				int distance = -1;
				if (!comp.car().equals(nullIndex) && !comp.cdr().equals(nullIndex) && mu[comp.car()].getChrom().equals(mu[comp.cdr()].getChrom())) {
					Region warpregion = mu[comp.car()].expand(1).combine(mu[comp.cdr()].expand(1));
					warpstring = warpregion.toString();
					distance = warpregion.getWidth();
				}
				out.println(left+"\t"+right+"\t"+tmpcount+"\t"+distance+"\t"+warpstring+"\t"+interll.get(comp));
			}
		}
		out.flush();
		out.close();
	}


	// inference with latent Z and fixed Chi 
	public void runStage3() throws FileNotFoundException {
		int numPairs = storage.getLeftN();
		float regionSize = region.getWidth();
		float backgroundProb = 1f/(regionSize*(regionSize+1f)/2f);
		float singleBackProb = 1f/regionSize;
		double logSingleBackProb = Math.log(singleBackProb);
		float localAlpha = 0;
		float localBeta = 0;
		int iter = 0;
		int subiter = 250;
		boolean converged = false;
		while (!converged) {
			totalLikelihood = 0;
			if (iter>=NO_ALPHA_ITER) {
				localAlpha = (((float)Math.min(iter-NO_ALPHA_ITER,ALPHA_ANNEAL_ITER)) / ((float)ALPHA_ANNEAL_ITER))*alpha;
				localBeta = (((float)Math.min(iter-NO_ALPHA_ITER,ALPHA_ANNEAL_ITER)) / ((float)ALPHA_ANNEAL_ITER))*beta;
			}
			System.err.println("alpha: "+localAlpha);
			System.err.println("beta: "+localBeta);
			//E Step
			for (int i=0; i<gamma3.length; i++) {
				float norm = 0;
				Map<Pair<Integer,Integer>,Float> tmpmap = gamma3[i];
				Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(i);
				for (Pair<Integer,Integer> comp : tmpmap.keySet()) { //iterating over components
					float compsum = 0;
					float psivalue = psi.get(comp);
					if (comp.car()==nullIndex || comp.cdr()==nullIndex) {
						float leftProb=0, rightProb=0;
						if (comp.car()==nullIndex) {
							leftProb = singleBackProb;
						} else {
							for (int j=0; j<chi.length; j++) {
								leftProb += gamma2[comp.car()][j]*dist.probability(currentPair.car(), mu[comp.car()], j);
							}
						}
						if (comp.cdr()==nullIndex) {
							rightProb = singleBackProb;
						} else {
							for (int j=0; j<chi.length; j++) {
								rightProb += gamma2[comp.cdr()][j]*dist.probability(currentPair.cdr(), mu[comp.cdr()], j);
							}
						}
						compsum += (1-rho)*psivalue*leftProb*rightProb;
					} else {
						for (int j=0; j<chi.length; j++) {
							for (int k=0; k<chi.length; k++) {
								float tmpprob = dist.probability(currentPair, mu[comp.car()], mu[comp.cdr()], j, k);
								compsum +=  (1-rho)*psivalue*gamma2[comp.car()][j]*gamma2[comp.cdr()][k]*tmpprob;
							}
						}
					}
					tmpmap.put(comp, compsum);
					norm += compsum;
				}
				for (Pair<Integer,Integer> comp : tmpmap.keySet()) {
					tmpmap.put(comp, tmpmap.get(comp)/norm);
				}
				if (dist.isSelfLigation(currentPair)) {
					float selfsum = 0;
					norm = 0;
					for (int j : gamma[i].keySet()) {
						if (j!=nullIndex) {
							selfsum = 0;
							for (int k=0; k<chi.length; k++) {
								selfsum += rho*pi[j]*gamma2[j][k]*dist.probability(currentPair, mu[j], null, k, nullIndex);
							}
							gamma[i].get(j)[0] = selfsum;
							norm += selfsum;
						}
					}
					float tmpback = rho*pi[nullIndex]*backgroundProb;
					gamma[i].get(nullIndex)[0] = tmpback;
					norm += tmpback;
					if (norm>0) {
						for (int j : gamma[i].keySet()) {
							gamma[i].get(j)[0] /= norm;
						}
					}
				}
			}
			//M Step
			System.err.println(deadComponents.size()+" dead components");
			for (int j=0; j<pi.length; j++) {
				pi[j] = 0;
			}
			for (Pair<Integer,Integer> comp : psi.keySet()) {
				psi.put(comp, 0f);
			}
			for (int i=0; i<numPairs; i++) {
				Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
				Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius);
				Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius);
				SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
				tmpleft = new Point(g,tmp.cdr().getChrom(), tmp.cdr().getLocation()-dist.radius);
				tmpright = new Point(g, tmp.cdr().getChrom(), tmp.cdr().getLocation()+dist.radius);
				SortedSet<Point> rightset = muset.subSet(tmpleft, tmpright);
				Pair<Integer,Integer> tmpcomp;
				for (Point left : leftset) {
					for (Point right : rightset) {
						tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),indexmap.get(right));
						psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
						tmpcomp = new Pair<Integer,Integer>(nullIndex,indexmap.get(right));
						psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
					}
					tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),nullIndex);
					psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
				}
				tmpcomp = new Pair<Integer,Integer>(nullIndex,nullIndex);
				psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
			}
			float psinorm = 0;
			int psialive = 0;
			for (Pair<Integer,Integer> comp : psi.keySet()) {
				float tmp = psi.get(comp);
				float newtmp = Math.max(0,tmp-localBeta);
				psinorm += newtmp;
				if (newtmp>0) psialive++;
				psi.put(comp, newtmp);
			}
			System.err.println(psialive+" psi components alive");
			for (Pair<Integer,Integer> comp : psi.keySet()) {
				psi.put(comp, psi.get(comp)/psinorm);
			}
			for (int n=0; n<numPairs; n++) {
				for (int j : gamma[n].keySet()) {
					pi[j] += gamma[n].get(j)[0];
				}
			}
			float pinorm = 0;
			for (int j=0; j<mu.length; j++) {
				pi[j] = Math.max(pi[j]-localAlpha, 0f);
				pinorm += pi[j];
			}
			pi[pi.length-1]= 0; 
			for (int j=0; j<numPairs; j++) {
				if (gamma[j].containsKey(nullIndex)) {
					pi[pi.length-1] += gamma[j].get(nullIndex)[0];
				}
			}

			pinorm += pi[pi.length-1];
			System.err.println("pinorm: "+pinorm);
			for (int j=0; j<pi.length; j++) {
				pi[j] /= pinorm;
			}
			rho = (pinorm+a) / (((float)numPairs)+a+b);

			double likelihood = 0;
			for (int j=0; j<mu.length; j++) {
				Point currentmu = mu[j];
				Region currentmuregion = currentmu.expand(dist.radius);
				Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(currentmuregion);
				int leftminn = leftminmaxn.car();
				int leftmaxn = leftminmaxn.cdr();
				for (int n = leftminn; n < leftmaxn; n++) {
					if (gamma[n].containsKey(j)) {
						Pair<StrandedPoint,StrandedPoint> r_n = storage.getLeftPair(n);
						for (int k=0; k<chi.length; k++) {
							float tmpprob = dist.probability(r_n, currentmu, null, k, nullIndex);
							if (tmpprob==0) {
								tmpprob = Float.MIN_VALUE;
							}
							likelihood += gamma[n].get(j)[0]*gamma2[j][k]*Math.log(tmpprob);
							if (Double.isNaN(likelihood)) {
								throw new FileNotFoundException(n+" "+j+" "+k+" "+gamma[n].get(j)[0]+" "+gamma2[j][k]+" "+tmpprob);
								//System.err.println(n+" "+j+" "+k+" "+gamma[n].get(j)[0]+" "+gamma2[j][k]+" "+tmpprob);
							}
						}
					}
				}
			}
			for (int i=0; i<numPairs; i++) {
				Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
				Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius);
				Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius);
				SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
				tmpleft = new Point(g,tmp.cdr().getChrom(), tmp.cdr().getLocation()-dist.radius);
				tmpright = new Point(g, tmp.cdr().getChrom(), tmp.cdr().getLocation()+dist.radius);
				SortedSet<Point> rightset = muset.subSet(tmpleft, tmpright);
				Pair<Integer,Integer> tmpcomp;
				float tmpgamma3;
				for (Point left : leftset) {
					int leftIndex = indexmap.get(left);
					for (Point right : rightset) {
						int rightIndex = indexmap.get(right);
						tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),indexmap.get(right));
						tmpgamma3 = gamma3[i].get(tmpcomp);
						for (int k1=0; k1<chi.length; k1++) {
							for (int k2=0; k2<chi.length; k2++) {
								float tmpprob = dist.probability(tmp, left, right, k1, k2);
								if (tmpprob==0) {
									tmpprob = Float.MIN_VALUE;
								}
								likelihood += tmpgamma3*gamma2[leftIndex][k1]*gamma2[rightIndex][k2]*Math.log(tmpprob);
								if (Double.isNaN(likelihood)) {
									throw new FileNotFoundException(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+rightIndex+" "+k2+" "+tmpgamma3+" "+gamma2[leftIndex][k1]+" "+gamma2[rightIndex][k2]+" "+tmpprob);
									//System.err.println(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+rightIndex+" "+k2+" "+tmpgamma3+" "+gamma2[leftIndex][k1]+" "+gamma2[rightIndex][k2]+" "+tmpprob);
								}

							}

						}
						tmpcomp = new Pair<Integer,Integer>(nullIndex,indexmap.get(right));
						tmpgamma3 = gamma3[i].get(tmpcomp);
						for (int k2=0; k2<chi.length; k2++) {
							float tmpprob = dist.probability(tmp.cdr(), right, k2);
							if (tmpprob==0) {
								tmpprob = Float.MIN_VALUE;
							}
							double prevlik = likelihood;
							likelihood += tmpgamma3*gamma2[rightIndex][k2]*(Math.log(tmpprob)+logSingleBackProb);
							if (Double.isNaN(likelihood)) {
								System.err.println(Math.log(tmpprob*singleBackProb));
								System.err.println(tmpgamma3*gamma2[rightIndex][k2]*Math.log(tmpprob*singleBackProb));
								throw new FileNotFoundException(i+" "+tmpcomp+" "+rightIndex+" "+k2+" "+tmpgamma3+" "+gamma2[rightIndex][k2]+" "+tmpprob+" "+singleBackProb+" "+prevlik);
								//System.err.println(i+" "+tmpcomp+" "+rightIndex+" "+k2+" "+tmpgamma3+" "+gamma2[rightIndex][k2]+" "+tmpprob+" "+singleBackProb);
							}
						}
					}
					tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),nullIndex);
					tmpgamma3 = gamma3[i].get(tmpcomp);
					for (int k1=0; k1<chi.length; k1++) {
						float tmpprob = dist.probability(tmp.car(), left, k1);
						if (tmpprob==0) {
							tmpprob = Float.MIN_VALUE;
						}
						likelihood += tmpgamma3*gamma2[leftIndex][k1]*(Math.log(tmpprob)+logSingleBackProb);
						if (Double.isNaN(likelihood)) {
							throw new FileNotFoundException(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+tmpgamma3+" "+gamma2[leftIndex][k1]+" "+tmpprob+" "+singleBackProb);
							//System.err.println(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+tmpgamma3+" "+gamma2[leftIndex][k1]+" "+tmpprob+" "+singleBackProb);
						}
					}
				}
				tmpcomp = new Pair<Integer,Integer>(nullIndex,nullIndex);
				tmpgamma3 = gamma3[i].get(tmpcomp);
				likelihood += tmpgamma3*2*logSingleBackProb;
				if (Double.isNaN(likelihood)) {
					throw new FileNotFoundException(i+" "+tmpcomp+" "+tmpgamma3+" "+singleBackProb);
					//System.err.println(i+" "+tmpcomp+" "+tmpgamma3+" "+singleBackProb);
				}
			}
			totalLikelihood = likelihood;

			for (int i=0; i<rhoHistory.length-1; i++) {
				rhoHistory[i] = rhoHistory[i+1];
			}
			rhoHistory[rhoHistory.length-1] = rho;
			for (int i=0; i<likelihoodHistory.length-1; i++) {
				likelihoodHistory[i] = likelihoodHistory[i+1];
			}
			likelihoodHistory[likelihoodHistory.length-1] = totalLikelihood;
			// check for convergence
			converged = true;
			if (iter==maxIters) {
				converged = false;
				break;
			}
			if (iter<NO_ALPHA_ITER+ALPHA_ANNEAL_ITER) {
				converged = false;
			} else {
				/*
				float[] minMaxRho = minMax(rhoHistory);
				double[] minMaxLikelihood = minMax(likelihoodHistory);
				if (minMaxRho[1]-minMaxRho[0]>RHO_ZONE || minMaxLikelihood[1]-minMaxLikelihood[0]>LIKELIHOOD_ZONE) {
					converged = false;
				}
				 */
				if (subiter>0) {
					subiter--;
					converged = false;
				} else {
					converged = true;
				}
			}
			iter++;
			System.err.println(iter+" iterations completed "+dfm.format(new Date()));
			dumpStage3();
		}
		if (converged) {
			System.err.println("converged after "+iter+" iterations "+dfm.format(new Date()));
		} else {
			System.err.println("did not converge "+dfm.format(new Date()));
		}
	}

	public void runLL() throws FileNotFoundException {
		int numPairs = storage.getLeftN();
		float regionSize = region.getWidth();
		float backgroundProb = 1f/(regionSize*(regionSize+1f)/2f);
		float singleBackProb = 1f/regionSize;
		double logSingleBackProb = Math.log(singleBackProb);
		float localAlpha = 0;
		float localBeta = alpha;
		totalLikelihood = beta;
		System.err.println("alpha: "+localAlpha);
		System.err.println("beta: "+localBeta);
		//E Step
		for (int i=0; i<gamma3.length; i++) {
			float norm = 0;
			Map<Pair<Integer,Integer>,Float> tmpmap = gamma3[i];
			Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(i);
			for (Pair<Integer,Integer> comp : tmpmap.keySet()) { //iterating over components
				float compsum = 0;
				float psivalue = psi.get(comp);
				if (comp.car()==nullIndex || comp.cdr()==nullIndex) {
					float leftProb=0, rightProb=0;
					if (comp.car()==nullIndex) {
						leftProb = singleBackProb;
					} else {
						for (int j=0; j<chi.length; j++) {
							leftProb += gamma2[comp.car()][j]*dist.probability(currentPair.car(), mu[comp.car()], j);
						}
					}
					if (comp.cdr()==nullIndex) {
						rightProb = singleBackProb;
					} else {
						for (int j=0; j<chi.length; j++) {
							rightProb += gamma2[comp.cdr()][j]*dist.probability(currentPair.cdr(), mu[comp.cdr()], j);
						}
					}
					compsum += (1-rho)*psivalue*leftProb*rightProb;
				} else {
					for (int j=0; j<chi.length; j++) {
						for (int k=0; k<chi.length; k++) {
							float tmpprob = dist.probability(currentPair, mu[comp.car()], mu[comp.cdr()], j, k);
							compsum +=  (1-rho)*psivalue*gamma2[comp.car()][j]*gamma2[comp.cdr()][k]*tmpprob;
						}
					}
				}
				tmpmap.put(comp, compsum);
				norm += compsum;
			}
			for (Pair<Integer,Integer> comp : tmpmap.keySet()) {
				tmpmap.put(comp, tmpmap.get(comp)/norm);
			}
			if (dist.isSelfLigation(currentPair)) {
				float selfsum = 0;
				norm = 0;
				for (int j : gamma[i].keySet()) {
					if (j!=nullIndex) {
						selfsum = 0;
						for (int k=0; k<chi.length; k++) {
							selfsum += rho*pi[j]*gamma2[j][k]*dist.probability(currentPair, mu[j], null, k, nullIndex);
						}
						gamma[i].get(j)[0] = selfsum;
						norm += selfsum;
					}
				}
				float tmpback = rho*pi[nullIndex]*backgroundProb;
				gamma[i].get(nullIndex)[0] = tmpback;
				norm += tmpback;
				if (norm>0) {
					for (int j : gamma[i].keySet()) {
						gamma[i].get(j)[0] /= norm;
					}
				}
			}
		}
		//M Step
		System.err.println(deadComponents.size()+" dead components");
		for (int j=0; j<pi.length; j++) {
			pi[j] = 0;
		}
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			psi.put(comp, 0f);
		}
		for (int i=0; i<numPairs; i++) {
			Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
			Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius);
			Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius);
			SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
			tmpleft = new Point(g,tmp.cdr().getChrom(), tmp.cdr().getLocation()-dist.radius);
			tmpright = new Point(g, tmp.cdr().getChrom(), tmp.cdr().getLocation()+dist.radius);
			SortedSet<Point> rightset = muset.subSet(tmpleft, tmpright);
			Pair<Integer,Integer> tmpcomp;
			for (Point left : leftset) {
				for (Point right : rightset) {
					tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),indexmap.get(right));
					psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
					tmpcomp = new Pair<Integer,Integer>(nullIndex,indexmap.get(right));
					psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
				}
				tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),nullIndex);
				psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
			}
			tmpcomp = new Pair<Integer,Integer>(nullIndex,nullIndex);
			psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
		}
		float psinorm = 0;
		int psialive = 0;
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			float tmp = psi.get(comp);
			float newtmp = Math.max(0,tmp-localBeta);
			psinorm += newtmp;
			if (newtmp>0) psialive++;
			psi.put(comp, newtmp);
		}
		System.err.println(psialive+" psi components alive");
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			psi.put(comp, psi.get(comp)/psinorm);
		}
		for (int n=0; n<numPairs; n++) {
			for (int j : gamma[n].keySet()) {
				pi[j] += gamma[n].get(j)[0];
			}
		}
		float pinorm = 0;
		for (int j=0; j<mu.length; j++) {
			pi[j] = Math.max(pi[j]-localAlpha, 0f);
			pinorm += pi[j];
		}
		pi[pi.length-1]= 0; 
		for (int j=0; j<numPairs; j++) {
			if (gamma[j].containsKey(nullIndex)) {
				pi[pi.length-1] += gamma[j].get(nullIndex)[0];
			}
		}

		pinorm += pi[pi.length-1];
		System.err.println("pinorm: "+pinorm);
		for (int j=0; j<pi.length; j++) {
			pi[j] /= pinorm;
		}
		rho = (pinorm+a) / (((float)numPairs)+a+b);

		for (int j=0; j<mu.length; j++) {
			Map<Integer,float[]>[] newgamma = (Map<Integer,float[]>[]) Array.newInstance(Map.class, storage.getLeftN());
			Map<Pair<Integer,Integer>,Float>[] newgamma3 = (Map<Pair<Integer,Integer>,Float>[]) Array.newInstance(Map.class, storage.getLeftN());
			for (int i = 0; i<gamma.length; i++) {
				float sum = 0;
				for (Integer jp : gamma[i].keySet()) {
					if (jp!=j) {
						float[] tmp = gamma[i].get(jp);
						sum += tmp[0];
						newgamma[i].put(jp,tmp);
					}
				}
				for (Integer jp : newgamma[i].keySet()) {
					newgamma[i].get(jp)[0] /= sum;
				}
				sum = 0;
				for (Pair<Integer,Integer> pair : gamma3[i].keySet()) {
					if (pair.car()!=j && pair.cdr()!=j) {
						float tmp = gamma3[i].get(pair);
						sum += tmp;
						newgamma3[i].put(pair,tmp);
					}
				}
				for (Pair<Integer,Integer> pair : newgamma3[i].keySet()) {
					newgamma3[i].put(pair, newgamma3[i].get(pair)/sum);
				}

			}
			double likelihood = 0;
			for (int jp=0; jp<mu.length; jp++) {
				Point currentmu = mu[jp];
				Region currentmuregion = currentmu.expand(dist.radius);
				Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(currentmuregion);
				int leftminn = leftminmaxn.car();
				int leftmaxn = leftminmaxn.cdr();
				for (int n = leftminn; n < leftmaxn; n++) {
					if (newgamma[n].containsKey(jp)) {
						Pair<StrandedPoint,StrandedPoint> r_n = storage.getLeftPair(n);
						for (int k=0; k<chi.length; k++) {
							float tmpprob = dist.probability(r_n, currentmu, null, k, nullIndex);
							if (tmpprob==0) {
								tmpprob = Float.MIN_VALUE;
							}
							likelihood += newgamma[n].get(jp)[0]*gamma2[jp][k]*Math.log(tmpprob);
							if (Double.isNaN(likelihood)) {
								throw new FileNotFoundException(n+" "+jp+" "+k+" "+newgamma[n].get(jp)[0]+" "+gamma2[jp][k]+" "+tmpprob);
								//System.err.println(n+" "+j+" "+k+" "+gamma[n].get(j)[0]+" "+gamma2[j][k]+" "+tmpprob);
							}
						}
					}
				}
			}
			for (int i=0; i<numPairs; i++) {
				Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
				Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius);
				Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius);
				SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
				tmpleft = new Point(g,tmp.cdr().getChrom(), tmp.cdr().getLocation()-dist.radius);
				tmpright = new Point(g, tmp.cdr().getChrom(), tmp.cdr().getLocation()+dist.radius);
				SortedSet<Point> rightset = muset.subSet(tmpleft, tmpright);
				Pair<Integer,Integer> tmpcomp;
				float tmpgamma3;
				for (Point left : leftset) {
					int leftIndex = indexmap.get(left);
					for (Point right : rightset) {
						int rightIndex = indexmap.get(right);
						tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),indexmap.get(right));
						tmpgamma3 = newgamma3[i].get(tmpcomp);
						for (int k1=0; k1<chi.length; k1++) {
							for (int k2=0; k2<chi.length; k2++) {
								float tmpprob = dist.probability(tmp, left, right, k1, k2);
								if (tmpprob==0) {
									tmpprob = Float.MIN_VALUE;
								}
								likelihood += tmpgamma3*gamma2[leftIndex][k1]*gamma2[rightIndex][k2]*Math.log(tmpprob);
								if (Double.isNaN(likelihood)) {
									throw new FileNotFoundException(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+rightIndex+" "+k2+" "+tmpgamma3+" "+gamma2[leftIndex][k1]+" "+gamma2[rightIndex][k2]+" "+tmpprob);
									//System.err.println(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+rightIndex+" "+k2+" "+tmpgamma3+" "+gamma2[leftIndex][k1]+" "+gamma2[rightIndex][k2]+" "+tmpprob);
								}

							}

						}
						tmpcomp = new Pair<Integer,Integer>(nullIndex,indexmap.get(right));
						tmpgamma3 = newgamma3[i].get(tmpcomp);
						for (int k2=0; k2<chi.length; k2++) {
							float tmpprob = dist.probability(tmp.cdr(), right, k2);
							if (tmpprob==0) {
								tmpprob = Float.MIN_VALUE;
							}
							double prevlik = likelihood;
							likelihood += tmpgamma3*gamma2[rightIndex][k2]*(Math.log(tmpprob)+logSingleBackProb);
							if (Double.isNaN(likelihood)) {
								System.err.println(Math.log(tmpprob*singleBackProb));
								System.err.println(tmpgamma3*gamma2[rightIndex][k2]*Math.log(tmpprob*singleBackProb));
								throw new FileNotFoundException(i+" "+tmpcomp+" "+rightIndex+" "+k2+" "+tmpgamma3+" "+gamma2[rightIndex][k2]+" "+tmpprob+" "+singleBackProb+" "+prevlik);
								//System.err.println(i+" "+tmpcomp+" "+rightIndex+" "+k2+" "+tmpgamma3+" "+gamma2[rightIndex][k2]+" "+tmpprob+" "+singleBackProb);
							}
						}
					}
					tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),nullIndex);
					tmpgamma3 = newgamma3[i].get(tmpcomp);
					for (int k1=0; k1<chi.length; k1++) {
						float tmpprob = dist.probability(tmp.car(), left, k1);
						if (tmpprob==0) {
							tmpprob = Float.MIN_VALUE;
						}
						likelihood += tmpgamma3*gamma2[leftIndex][k1]*(Math.log(tmpprob)+logSingleBackProb);
						if (Double.isNaN(likelihood)) {
							throw new FileNotFoundException(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+tmpgamma3+" "+gamma2[leftIndex][k1]+" "+tmpprob+" "+singleBackProb);
							//System.err.println(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+tmpgamma3+" "+gamma2[leftIndex][k1]+" "+tmpprob+" "+singleBackProb);
						}
					}
				}
				tmpcomp = new Pair<Integer,Integer>(nullIndex,nullIndex);
				tmpgamma3 = newgamma3[i].get(tmpcomp);
				likelihood += tmpgamma3*2*logSingleBackProb;
				if (Double.isNaN(likelihood)) {
					throw new FileNotFoundException(i+" "+tmpcomp+" "+tmpgamma3+" "+singleBackProb);
					//System.err.println(i+" "+tmpcomp+" "+tmpgamma3+" "+singleBackProb);
				}
			}
			ll[j] = likelihood;
		}


	}

	public void runInterLL() throws FileNotFoundException {
		int numPairs = storage.getLeftN();
		float regionSize = region.getWidth();
		float backgroundProb = 1f/(regionSize*(regionSize+1f)/2f);
		float singleBackProb = 1f/regionSize;
		System.err.println("singleBackProb "+singleBackProb);
		double logSingleBackProb = Math.log(singleBackProb);
		float localAlpha = alpha;
		float localBeta = beta;
		totalLikelihood = 0;
		System.err.println("alpha: "+localAlpha);
		System.err.println("beta: "+localBeta);
		//E Step
		for (int i=0; i<gamma3.length; i++) {
			float norm = 0;
			Map<Pair<Integer,Integer>,Float> tmpmap = gamma3[i];
			Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(i);
			for (Pair<Integer,Integer> comp : tmpmap.keySet()) { //iterating over components
				float compsum = 0;
				float psivalue = psi.get(comp);
				if (psivalue==0) {
					System.err.println("0 psivalue");
				}
				if (comp.car()==nullIndex || comp.cdr()==nullIndex) {
					float leftProb=0, rightProb=0;
					if (comp.car()==nullIndex) {
						leftProb = singleBackProb;
					} else {
						for (int j=0; j<chi.length; j++) {
							leftProb += gamma2[comp.car()][j]*dist.probability(currentPair.car(), mu[comp.car()], j);
						}
					}
					if (comp.cdr()==nullIndex) {
						rightProb = singleBackProb;
					} else {
						for (int j=0; j<chi.length; j++) {
							rightProb += gamma2[comp.cdr()][j]*dist.probability(currentPair.cdr(), mu[comp.cdr()], j);
							if (rightProb==0) {
								//System.err.println("F "+gamma2[comp.cdr()][j]+" "+dist.probability(currentPair.cdr(), mu[comp.cdr()], j)+" "+comp+" "+j+" "+currentPair+" "+mu[comp.cdr()]);
							}
							if (Float.isNaN(rightProb)) {
								System.err.println("D "+gamma2[comp.cdr()][j]+" "+dist.probability(currentPair.cdr(), mu[comp.cdr()], j)+" "+comp+" "+j+" "+currentPair+" "+mu[comp.cdr()]);
							}
						}
					}
					compsum += (1-rho)*psivalue*leftProb*rightProb;
					if (compsum==0) {
						//System.err.println("C "+rho+" "+psivalue+" "+leftProb+" "+rightProb+" "+comp);
					}
					if (Float.isNaN(compsum)) {
						System.err.println("A "+rho+" "+psivalue+" "+leftProb+" "+rightProb+" "+comp);
					}
				} else {
					for (int j=0; j<chi.length; j++) {
						for (int k=0; k<chi.length; k++) {
							float tmpprob = dist.probability(currentPair, mu[comp.car()], mu[comp.cdr()], j, k);
							compsum +=  (1-rho)*psivalue*gamma2[comp.car()][j]*gamma2[comp.cdr()][k]*tmpprob;
							if (Float.isNaN(compsum)) {
								System.err.println("B "+rho+" "+psivalue+" "+gamma2[comp.car()][j]+" "+gamma2[comp.cdr()][k]+" "+tmpprob+" "+comp);
							}
						}
					}
				}
				tmpmap.put(comp, compsum);
				norm += compsum;
			}
			if (norm==0) {
				System.err.print(storage.getLeftPair(i));
				for (Pair<Integer,Integer> pair : gamma3[i].keySet()) {
					System.err.print(" "+pair);
				}
				System.err.println(" nullIndex: "+nullIndex);
			} else {
				for (Pair<Integer,Integer> comp : tmpmap.keySet()) {
					tmpmap.put(comp, tmpmap.get(comp)/norm);
				}
			}
			if (dist.isSelfLigation(currentPair)) {
				float selfsum = 0;
				norm = 0;
				for (int j : gamma[i].keySet()) {
					if (j!=nullIndex) {
						selfsum = 0;
						for (int k=0; k<chi.length; k++) {
							selfsum += rho*pi[j]*gamma2[j][k]*dist.probability(currentPair, mu[j], null, k, nullIndex);
						}
						gamma[i].get(j)[0] = selfsum;
						norm += selfsum;
					}
				}
				float tmpback = rho*pi[nullIndex]*backgroundProb;
				gamma[i].get(nullIndex)[0] = tmpback;
				norm += tmpback;
				if (norm>0) {
					for (int j : gamma[i].keySet()) {
						gamma[i].get(j)[0] /= norm;
					}
				}
			}
		}
		//M Step
		/*
		System.err.println(deadComponents.size()+" dead components");
		for (int j=0; j<pi.length; j++) {
			pi[j] = 0;
		}
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			psi.put(comp, 0f);
		}
		for (int i=0; i<numPairs; i++) {
			Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
			Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius);
			Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius);
			SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
			tmpleft = new Point(g,tmp.cdr().getChrom(), tmp.cdr().getLocation()-dist.radius);
			tmpright = new Point(g, tmp.cdr().getChrom(), tmp.cdr().getLocation()+dist.radius);
			SortedSet<Point> rightset = muset.subSet(tmpleft, tmpright);
			Pair<Integer,Integer> tmpcomp;
			for (Point left : leftset) {
				for (Point right : rightset) {
					tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),indexmap.get(right));
					if (psi.containsKey(tmpcomp)) {
						psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
					}
					tmpcomp = new Pair<Integer,Integer>(nullIndex,indexmap.get(right));
					if (psi.containsKey(tmpcomp)) {
						psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
					}
				}
				tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),nullIndex);
				if (psi.containsKey(tmpcomp)) {
					psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
				}
			}
			tmpcomp = new Pair<Integer,Integer>(nullIndex,nullIndex);
			if (psi.containsKey(tmpcomp)) {
				psi.put(tmpcomp, psi.get(tmpcomp)+gamma3[i].get(tmpcomp));
			}
		}
		float psinorm = 0;
		int psialive = 0;
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			float tmp = psi.get(comp);
			float newtmp = Math.max(0,tmp-localBeta);
			psinorm += newtmp;
			if (newtmp>0) psialive++;
			psi.put(comp, newtmp);
		}
		System.err.println(psialive+" psi components alive");
		for (Pair<Integer,Integer> comp : psi.keySet()) {
			psi.put(comp, psi.get(comp)/psinorm);
		}
		for (int n=0; n<numPairs; n++) {
			for (int j : gamma[n].keySet()) {
				pi[j] += gamma[n].get(j)[0];
			}
		}
		float pinorm = 0;
		for (int j=0; j<mu.length; j++) {
			pi[j] = Math.max(pi[j]-localAlpha, 0f);
			pinorm += pi[j];
		}
		pi[pi.length-1]= 0; 
		for (int j=0; j<numPairs; j++) {
			if (gamma[j].containsKey(nullIndex)) {
				pi[pi.length-1] += gamma[j].get(nullIndex)[0];
			}
		}

		pinorm += pi[pi.length-1];
		System.err.println("pinorm: "+pinorm);
		for (int j=0; j<pi.length; j++) {
			pi[j] /= pinorm;
		}
		rho = (pinorm+a) / (((float)numPairs)+a+b);
		 */

		int iter = 0;
		int div = Math.max(psi.size()/100,1);
		for (Pair<Integer,Integer> toremove : psi.keySet()) {
			if (!toremove.car().equals(nullIndex) && !toremove.cdr().equals(nullIndex)) {
				Point leftmu = mu[toremove.car()];
				Point rightmu = mu[toremove.cdr()];
				Point lessermu;
				if (leftmu.compareTo(rightmu)>0) {
					lessermu = rightmu;
				} else {
					lessermu = leftmu;
				}
				Region leftregion = lessermu.expand(dist.radius);
				Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(leftregion);
				int leftminn = leftminmaxn.car();
				int leftmaxn = leftminmaxn.cdr();
				Map<Pair<Integer,Integer>,Float>[] newgamma3 = (Map<Pair<Integer,Integer>,Float>[]) Array.newInstance(Map.class, storage.getLeftN());
				for (int i = leftminn; i<leftmaxn; i++) {
					newgamma3[i] = new HashMap<Pair<Integer,Integer>,Float>();
					float sum = 0;
					for (Pair<Integer,Integer> pair : gamma3[i].keySet()) {
						if (!pair.equals(toremove)) {
							float tmp = gamma3[i].get(pair);
							sum += tmp;
							newgamma3[i].put(pair,tmp);
						}
					}
					if (sum==0) {
						System.err.print(toremove+" "+storage.getLeftPair(i));
						for (Pair<Integer,Integer> pair : gamma3[i].keySet()) {
							System.err.print(" "+pair);
						}
						System.err.println();
					} else {
						for (Pair<Integer,Integer> pair : newgamma3[i].keySet()) {
							newgamma3[i].put(pair, newgamma3[i].get(pair)/sum);
						}
					}

				}
				double newlikelihood = 0;
				double oldlikelihood = 0;
				/*
			for (int jp=0; jp<mu.length; jp++) {
				Point currentmu = mu[jp];
				Region currentmuregion = currentmu.expand(dist.radius);
				Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(currentmuregion);
				int leftminn = leftminmaxn.car();
				int leftmaxn = leftminmaxn.cdr();
				for (int n = leftminn; n < leftmaxn; n++) {
					if (gamma[n].containsKey(jp)) {
						Pair<StrandedPoint,StrandedPoint> r_n = storage.getLeftPair(n);
						for (int k=0; k<chi.length; k++) {
							float tmpprob = dist.probability(r_n, currentmu, null, k, nullIndex);
							if (tmpprob==0) {
								tmpprob = Float.MIN_VALUE;
							}
							likelihood += gamma[n].get(jp)[0]*gamma2[jp][k]*Math.log(tmpprob);
							if (Double.isNaN(likelihood)) {
								throw new FileNotFoundException(n+" "+jp+" "+k+" "+gamma[n].get(jp)[0]+" "+gamma2[jp][k]+" "+tmpprob);
								//System.err.println(n+" "+j+" "+k+" "+gamma[n].get(j)[0]+" "+gamma2[j][k]+" "+tmpprob);
							}
						}
					}
				}
			}
				 */
				for (int i = leftminn; i<leftmaxn; i++) {
					Pair<StrandedPoint,StrandedPoint> tmp = storage.getLeftPair(i);
					Point tmpleft = new Point(g,tmp.car().getChrom(), tmp.car().getLocation()-dist.radius);
					Point tmpright = new Point(g, tmp.car().getChrom(), tmp.car().getLocation()+dist.radius);
					SortedSet<Point> leftset = muset.subSet(tmpleft, tmpright);
					tmpleft = new Point(g,tmp.cdr().getChrom(), tmp.cdr().getLocation()-dist.radius);
					tmpright = new Point(g, tmp.cdr().getChrom(), tmp.cdr().getLocation()+dist.radius);
					SortedSet<Point> rightset = muset.subSet(tmpleft, tmpright);
					Pair<Integer,Integer> tmpcomp;
					float tmpgamma3new;
					float tmpgamma3old;
					for (Point left : leftset) {
						int leftIndex = indexmap.get(left);
						for (Point right : rightset) {
							int rightIndex = indexmap.get(right);
							tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),indexmap.get(right));
							if (psi.containsKey(tmpcomp) && !tmpcomp.equals(toremove)) {
								tmpgamma3new = newgamma3[i].get(tmpcomp);
								tmpgamma3old = gamma3[i].get(tmpcomp);
								for (int k1=0; k1<chi.length; k1++) {
									for (int k2=0; k2<chi.length; k2++) {
										float tmpprob = dist.probability(tmp, left, right, k1, k2);
										if (tmpprob==0) {
											tmpprob = Float.MIN_VALUE;
										}
										newlikelihood += tmpgamma3new*gamma2[leftIndex][k1]*gamma2[rightIndex][k2]*Math.log(tmpprob);
										oldlikelihood += tmpgamma3old*gamma2[leftIndex][k1]*gamma2[rightIndex][k2]*Math.log(tmpprob);
										if (Double.isNaN(newlikelihood)) {
											throw new FileNotFoundException(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+rightIndex+" "+k2+" "+tmpgamma3new+" "+gamma2[leftIndex][k1]+" "+gamma2[rightIndex][k2]+" "+tmpprob);
											//System.err.println(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+rightIndex+" "+k2+" "+tmpgamma3+" "+gamma2[leftIndex][k1]+" "+gamma2[rightIndex][k2]+" "+tmpprob);
										}

									}

								}
							}
							tmpcomp = new Pair<Integer,Integer>(nullIndex,indexmap.get(right));
							if (psi.containsKey(tmpcomp) && !tmpcomp.equals(toremove)) {
								tmpgamma3new = newgamma3[i].get(tmpcomp);
								tmpgamma3old = gamma3[i].get(tmpcomp);
								for (int k2=0; k2<chi.length; k2++) {
									float tmpprob = dist.probability(tmp.cdr(), right, k2);
									if (tmpprob==0) {
										tmpprob = Float.MIN_VALUE;
									}
									double prevlik = newlikelihood;
									newlikelihood += tmpgamma3new*gamma2[rightIndex][k2]*(Math.log(tmpprob)+logSingleBackProb);
									oldlikelihood += tmpgamma3old*gamma2[rightIndex][k2]*(Math.log(tmpprob)+logSingleBackProb);
									if (Double.isNaN(newlikelihood)) {
										System.err.println(Math.log(tmpprob*singleBackProb));
										System.err.println(tmpgamma3new*gamma2[rightIndex][k2]*Math.log(tmpprob*singleBackProb));
										throw new FileNotFoundException(i+" "+tmpcomp+" "+rightIndex+" "+k2+" "+tmpgamma3new+" "+gamma2[rightIndex][k2]+" "+tmpprob+" "+singleBackProb+" "+prevlik);
										//System.err.println(i+" "+tmpcomp+" "+rightIndex+" "+k2+" "+tmpgamma3+" "+gamma2[rightIndex][k2]+" "+tmpprob+" "+singleBackProb);
									}
								}
							}
						}
						tmpcomp = new Pair<Integer,Integer>(indexmap.get(left),nullIndex);
						if (psi.containsKey(tmpcomp) && !tmpcomp.equals(toremove)) {
							tmpgamma3new = newgamma3[i].get(tmpcomp);
							tmpgamma3old = gamma3[i].get(tmpcomp);
							for (int k1=0; k1<chi.length; k1++) {
								float tmpprob = dist.probability(tmp.car(), left, k1);
								if (tmpprob==0) {
									tmpprob = Float.MIN_VALUE;
								}
								newlikelihood += tmpgamma3new*gamma2[leftIndex][k1]*(Math.log(tmpprob)+logSingleBackProb);
								oldlikelihood += tmpgamma3old*gamma2[leftIndex][k1]*(Math.log(tmpprob)+logSingleBackProb);
								if (Double.isNaN(newlikelihood)) {
									throw new FileNotFoundException(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+tmpgamma3new+" "+gamma2[leftIndex][k1]+" "+tmpprob+" "+singleBackProb);
									//System.err.println(i+" "+tmpcomp+" "+leftIndex+" "+k1+" "+tmpgamma3+" "+gamma2[leftIndex][k1]+" "+tmpprob+" "+singleBackProb);
								}
							}
						}
					}
					tmpcomp = new Pair<Integer,Integer>(nullIndex,nullIndex);
					if (psi.containsKey(tmpcomp) && !tmpcomp.equals(toremove)) {
						tmpgamma3new = newgamma3[i].get(tmpcomp);
						tmpgamma3old = gamma3[i].get(tmpcomp);
						newlikelihood += tmpgamma3new*2*logSingleBackProb;
						oldlikelihood += tmpgamma3old*2*logSingleBackProb;
						if (Double.isNaN(newlikelihood)) {
							throw new FileNotFoundException(i+" "+tmpcomp+" "+tmpgamma3new+" "+singleBackProb);
							//System.err.println(i+" "+tmpcomp+" "+tmpgamma3+" "+singleBackProb);
						}
					}
				}
				interll.put(toremove, newlikelihood-oldlikelihood);
			}
			iter++;
			if (iter % div == 0) {
				System.err.println(iter+" iterations completed "+dfm.format(new Date()));
			}
		}


	}

	public void runStage2Hybrid(int numcycles) throws FileNotFoundException {
		for (int i=0; i<numcycles; i++) {
			runStage1Hybrid(5);
			runStage2(10);
		}
	}

	// inference with latent Tau and fixed Z (or gamma?) 
	public void runStage2(int numiters) throws FileNotFoundException {
		int numPairs = storage.getLeftN();
		int iter = 0;
		boolean converged = false;
		while (!converged) {
			totalLikelihood = 0;
			for (int j=0; j<mu.length; j++) {
				oldmu[j] = mu[j];
			}
			for (int j=0; j<chi.length; j++) {
				oldchi[j] = chi[j];
			}
			// E Step
			for (int i=0; i<gamma2.length; i++) {
				for (int j=0; j<gamma2[i].length; j++) {
					gamma2[i][j] = 0;
				}
			}
			for (int n=0; n<numPairs; n++) {
				Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(n);
				for (int i : gamma[n].keySet()) {
					if (i!=nullIndex) {
						float[] tmp = gamma[n].get(i);
						for (int j=0; j<gamma2[i].length; j++) {
							float tmpprob = dist.probability(currentPair, mu[i], null, j, nullIndex);
							gamma2[i][j] += tmp[0]*chi[j]*tmpprob;
							gamma2[i][j] += tmp[1]*chi[j]*dist.probability(currentPair, mu[i], mu[i], j, j);
						}
					}
				}
			}
			for (int i=0; i<gamma2.length; i++) {
				float norm = 0;
				for (int j=0; j<gamma2[i].length; j++) {
					norm += gamma2[i][j];
				}
				if (norm>0) {
					for (int j=0; j<gamma2[i].length; j++) {
						gamma2[i][j] /= norm;
					}
				}
			}
			// M Step
			/*
			System.err.println(deadComponents.size()+" dead components");
			for (int j=0; j<mu.length; j++) {
				if (componentAliveStage1(j)) {
					Point currentmu = mu[j];
					int bestx = -1;
					float bestlikelihood = Float.NEGATIVE_INFINITY;
					Region currentmuregion = currentmu.expand(dist.radius);
					Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(currentmuregion);
					int leftminn = leftminmaxn.car();
					int leftmaxn = leftminmaxn.cdr();
					for (int x = Math.max(0,currentmu.getLocation()); x <= currentmu.getLocation(); x += 1) {
						Point xpoint = new Point(g,mu[j].getChrom(), x);
						float likelihood = 0;
						for (int n = leftminn; n < leftmaxn; n++) {
							Pair<StrandedPoint,StrandedPoint> r_n = storage.getLeftPair(n);
							for (int k=0; k<gamma2[j].length; k++) {
								float tmpprob = dist.probability(r_n, xpoint, null, k, nullIndex);
								if (tmpprob==0) {
									tmpprob = Float.MIN_VALUE;
								}
								likelihood += gamma[n][j][0]*gamma2[j][k]*Math.log(tmpprob);
								tmpprob = dist.probability(r_n, xpoint, xpoint, k, k);
								if (tmpprob==0) {
									tmpprob = Float.MIN_VALUE;
								}
								likelihood += gamma[n][j][1]*gamma2[j][k]*Math.log(tmpprob);
							}
						}
						if (likelihood>bestlikelihood) {
							bestx = x;
							bestlikelihood = likelihood;
						}
					}
					if (bestx!=-1) {
						totalLikelihood += bestlikelihood;
						Point newpoint = new Point(currentmu.getGenome(), currentmu.getChrom(), bestx);
						mu[j] = newpoint;
					}
				}
			}
			 */
			/*
			float[] N = new float[gamma2[0].length];
			float Nnorm = 0;
			for (int j=0; j<gamma2.length; j++) {
				for (int k=0; k<gamma2[j].length; k++) {
					N[k] += gamma2[j][k];
					Nnorm += gamma2[j][k];
				}
			}
			for (int k=0; k<gamma2[0].length; k++) {
				chi[k] = N[k] / Nnorm;
			}
			 */
			// convergence?
			//converged = true;
			if (iter==maxIters || iter==numiters) {
				converged = false;
				break;
			}
			iter++;
			System.err.println(iter+" iterations completed "+dfm.format(new Date()));
			dumpStage2();
		}
		if (converged) {
			System.err.println("converged after "+iter+" iterations "+dfm.format(new Date()));
		} else {
			System.err.println("did not converge "+dfm.format(new Date()));
		}
	}


	private void parseGamma2File(String file) throws IOException {
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		for (int i=0; i<mu.length; i++) {
			s = r.readLine();
			split = s.split("\t");
			for (int j=0; j<chi.length; j++) {
				gamma2[i][j] = Float.valueOf(split[j]);
			}
		}
		r.close();
	}


	private void computeGammaForStage3() {
		float regionSize = region.getWidth();
		float backgroundProb = 1f/(regionSize*(regionSize+1f)/2f);
		for (int i=0; i<gamma.length; i++) {
			gamma[i] = new HashMap<Integer,float[]>();
			float norm = 0;
			Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(i);
			for (int j=0; j<mu.length; j++) { //iterating over components
				float tmpprob = genericDist.probability(currentPair, mu[j], null, 0, nullIndex);
				if (tmpprob>0) {
					float[] tmp = {rho*pi[j]*tmpprob};
					gamma[i].put(j, tmp);
					norm += tmp[0];
				}
			}
			if (genericDist.isSelfLigation(currentPair)) {
				float[] tmp = {rho*pi[nullIndex]*backgroundProb};
				gamma[i].put(nullIndex, tmp);
				norm += tmp[0];
			}
			if (norm>0) {
				for (int j : gamma[i].keySet()) {
					gamma[i].get(j)[0] /= norm;
				}
			}
		}
	}

	private void computeGammaForStage3Hybrid() {
		float regionSize = region.getWidth();
		float backgroundProb = 1f/(regionSize*(regionSize+1f)/2f);
		for (int i=0; i<gamma.length; i++) {
			gamma[i] = new HashMap<Integer,float[]>();
			float norm = 0;
			Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(i);
			for (int j=0; j<mu.length; j++) { //iterating over components
				float tmp[] = {0};
				for (int k=0; k<chi.length; k++) {
					float tmpprob = dist.probability(currentPair, mu[j], null, k, nullIndex);
					tmp[0] += rho*chi[k]*pi[j]*tmpprob;
				}
				if (tmp[0]>0) {
					gamma[i].put(j, tmp);
					norm += tmp[0];
				}
			}
			if (dist.isSelfLigation(currentPair)) {
				float[] tmp = {rho*pi[nullIndex]*backgroundProb};
				gamma[i].put(nullIndex, tmp);
				norm += tmp[0];
			}
			if (norm>0) {
				for (int j : gamma[i].keySet()) {
					gamma[i].get(j)[0] /= norm;
				}
			}
		}
	}

	private void computeGamma2ForStage3() {
		int numPairs = storage.getLeftN();
		for (int n=0; n<numPairs; n++) {
			Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(n);
			for (int i : gamma[n].keySet()) {
				if (i!=nullIndex) {
					float[] tmp = gamma[n].get(i);
					for (int j=0; j<gamma2[i].length; j++) {
						gamma2[i][j] += tmp[0]*chi[j]*dist.probability(currentPair, mu[i], null, j, nullIndex);
						//gamma2[i][j] += tmp[1]*chi[j]*dist.probability(currentPair, mu[i], mu[i], j, j);
					}
				}
			}
		}
		for (int i=0; i<gamma2.length; i++) {
			float norm = 0;
			for (int j=0; j<gamma2[i].length; j++) {
				norm += gamma2[i][j];
			}
			if (norm==0) {
				for (int j=0; j<gamma2[i].length; j++) {
					gamma2[i][j] = 1f/3f;
				}
			} else {
				for (int j=0; j<gamma2[i].length; j++) {
					gamma2[i][j] /= norm;
				}
			}
		}
	}

	private void computeStage1Gamma() {
		float regionSize = region.getWidth();
		float backgroundProb = 1f/(regionSize*(regionSize+1f)/2f);
		for (int i=0; i<gamma.length; i++) {
			gamma[i] = new HashMap<Integer,float[]>();
			float norm = 0;
			Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(i);
			for (int j=0; j<mu.length; j++) { //iterating over components
				float tmpprob = genericDist.probability(currentPair, mu[j], null, 0, nullIndex);
				if (tmpprob>0) {
					float[] tmp = new float[2];
					tmp[0] = rho*pi[j]*tmpprob;
					norm += tmp[0];
					tmpprob = genericDist.probability(currentPair, mu[j], mu[j], 0, 0);
					tmp[1] = (1f-rho)*pi[j]*tmpprob;
					norm += tmp[1];
					gamma[i].put(j, tmp);
				}
			}
			float[] tmp = {pi[nullIndex]*backgroundProb};
			gamma[i].put(nullIndex, tmp);
			norm += tmp[0];
			if (norm>0) {
				for (int j : gamma[i].keySet()) {
					tmp = gamma[i].get(j);
					tmp[0] /= norm;
					if (j!=nullIndex) {
						tmp[1] /= norm;
					}
				}
			}
		}
	}


	/* inference with one dimensional latent Z and fixed tau with generic read distributions */
	public void runStage1(int numiters) throws FileNotFoundException {
		float regionSize = region.getWidth();
		float backgroundProb = 1f/(regionSize*(regionSize+1f)/2f);
		float localAlpha = 0;
		int localIncrement = X_INCREMENT;
		int iter = 0;
		int subiter = 250;
		if (numiters>-1) {
			iter = NO_ALPHA_ITER+ALPHA_ANNEAL_ITER;
			localAlpha = alpha;
			subiter = numiters;
			localIncrement = 1;
		}
		boolean converged = false;
		double likmax = -Float.MAX_VALUE;
		int maxcount = 0;
		int lastseen = 0;
		while (!converged) {
			totalLikelihood = 0;
			/*
			for (int j=0; j<mu.length; j++) {
				oldmu[j] = mu[j];
			}
			 */
			if (iter>=NO_ALPHA_ITER) {
				localAlpha = (((float)Math.min(iter-NO_ALPHA_ITER,ALPHA_ANNEAL_ITER)) / ((float)ALPHA_ANNEAL_ITER))*alpha;
			}
			System.err.println("alpha: "+localAlpha);
			//E Step
			for (int i=0; i<gamma.length; i++) {
				float norm = 0;
				Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(i);
				float[] tmpgamma;
				for (int j : gamma[i].keySet()) { //iterating over components
					if (j!=nullIndex) {
						tmpgamma = gamma[i].get(j);
						float tmpprob = dist.probability(currentPair, mu[j], null, tau[j], nullIndex);
						tmpgamma[0] = rho*pi[j]*tmpprob;
						norm += tmpgamma[0];
						tmpprob = dist.probability(currentPair, mu[j], mu[j], tau[j], tau[j]);
						tmpgamma[1] = (1f-rho)*pi[j]*tmpprob;
						norm += tmpgamma[1];
					}
				}
				tmpgamma = gamma[i].get(nullIndex);
				tmpgamma[0] = pi[nullIndex]*backgroundProb;
				norm += tmpgamma[0];
				if (norm>0) {
					for (int j : gamma[i].keySet()) {
						if (j==nullIndex) {
							gamma[i].get(j)[0] /= norm;
						} else {
							tmpgamma = gamma[i].get(j);
							tmpgamma[0] /= norm;
							tmpgamma[1] /= norm;
						}
					}
				}
			}
			//M Step
			float[] N_pi = new float[mu.length];
			float[] N = new float[mu.length];
			//mu (and tau?)
			System.err.println(deadComponents.size()+" dead components");
			for (int j=0; j<mu.length; j++) {
				if (componentAliveStage1(j)) {
					Point currentmu = mu[j];
					int bestx = -1;
					float bestlikelihood = Float.NEGATIVE_INFINITY;
					Region currentmuregion = currentmu.expand(dist.radius);
					Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(currentmuregion);
					int leftminn = leftminmaxn.car();
					int leftmaxn = leftminmaxn.cdr();
					for (int x = Math.max(0,currentmu.getLocation()-2*dist.radius); x <= currentmu.getLocation()+2*dist.radius; x += localIncrement) {
						Point xpoint = new Point(g,mu[j].getChrom(), x);
						float likelihood = 0;
						for (int n = leftminn; n < leftmaxn; n++) {
							Pair<StrandedPoint,StrandedPoint> r_n = storage.getLeftPair(n);
							float tmpprob = dist.probability(r_n, xpoint, null, tau[j], nullIndex);
							if (tmpprob==0) {
								tmpprob = Float.MIN_VALUE;
							}
							if (gamma[n].get(j)==null) {
								System.err.println(n+" "+r_n.car()+" "+j+" "+mu[j]);
							}
							likelihood += gamma[n].get(j)[0]*Math.log(tmpprob);
							tmpprob = dist.probability(r_n, xpoint, xpoint, tau[j], tau[j]);
							if (tmpprob==0) {
								tmpprob = Float.MIN_VALUE;
							}
							likelihood += gamma[n].get(j)[1]*Math.log(tmpprob);
						}
						if (likelihood>bestlikelihood) {
							bestx = x;
							bestlikelihood = likelihood;
						}
					}
					if (bestx!=-1) {
						totalLikelihood += bestlikelihood;
						Point newpoint = new Point(currentmu.getGenome(), currentmu.getChrom(), bestx);
						mu[j] = newpoint;
						updateHashGamma(newpoint,j);
					}
					for (int n = leftminn; n < leftmaxn; n++) {
						for (int k=0; k<2; k++) {
							N[j] += gamma[n].get(j)[k];

						}
						N_pi[j] += gamma[n].get(j)[0];
					}
				}
			}
			float N_piall = 0;
			float pinorm = 0;
			for (int j=0; j<N.length; j++) {
				N_piall += Math.max(N_pi[j]-localAlpha, 0f);
				pi[j] = Math.max(N[j]-localAlpha, 0f);
				pinorm += pi[j];
			}
			int numpairs = storage.getLeftN();
			pi[pi.length-1]= 0; 
			for (int j=0; j<numpairs; j++) {
				pi[pi.length-1] += gamma[j].get(nullIndex)[0]; 
			}
			rho = (N_piall+a) / (((float)pinorm)+a+b);
			pinorm += pi[pi.length-1];
			System.err.println("pinorm: "+pinorm);
			for (int j=0; j<pi.length; j++) {
				pi[j] /= pinorm;
			}

			for (int i=0; i<rhoHistory.length-1; i++) {
				rhoHistory[i] = rhoHistory[i+1];
			}
			rhoHistory[rhoHistory.length-1] = rho;
			for (int i=0; i<likelihoodHistory.length-1; i++) {
				likelihoodHistory[i] = likelihoodHistory[i+1];
			}
			likelihoodHistory[likelihoodHistory.length-1] = totalLikelihood;
			likelihoodList.add(totalLikelihood);
			// check for convergence
			converged = true;
			if (iter==maxIters) {
				converged = false;
				break;
			}
			if (iter<NO_ALPHA_ITER+ALPHA_ANNEAL_ITER) {
				converged = false;
			} else {
				/*
				if (totalLikelihood>likmax) {
					lastseen = iter;
					likmax = totalLikelihood;
					maxcount = 0;
					converged = false;
				} else if (totalLikelihood==likmax) {
					lastseen = iter;
					maxcount++;
					if (maxcount<MINMAXCOUNT) {
						converged = false;
					}
				} else if (iter-lastseen>MINMAXCOUNT) {
					converged = false;
					likmax = localMax(likelihoodList, MINMAXCOUNT-2);
					lastseen = iter;
					maxcount = 0;
				} else {
					converged = false;
				}
				 */
				/*
				for (int i=3; i<likelihoodHistory.length; i += 2) {
					if (likelihoodHistory[i]!=likelihoodHistory[i-2]) {
						converged = false;
					}
				}
				 */

				if (subiter>0) {
					subiter--;
					converged = false;
				} else {
					converged = true;
				}

				/*
				float[] minMaxRho = minMax(rhoHistory);
				float[] minMaxLikelihood = minMax(likelihoodHistory);
				if (minMaxRho[1]-minMaxRho[0]>RHO_ZONE || minMaxLikelihood[1]-minMaxLikelihood[0]>LIKELIHOOD_ZONE) {
					converged = false;
				}
				 */
			}
			if (converged && localIncrement>1) {
				lastseen = iter;
				likmax = -Float.MAX_VALUE;
				maxcount = 0;
				converged = false;
				subiter = 250;
				localIncrement = Math.max(localIncrement/10, 1);
				System.err.println("Increment reduced to "+localIncrement);
			}
			iter++;
			System.err.println(iter+" iterations completed "+dfm.format(new Date()));
			dumpStage1();
			/*
			int[][] muarr = new int[1][mu.length];
			for (int i=0; i<mu.length; i++) {
				muarr[0][i] = mu[i].getLocation();
			}
			double[][] piarr = new double[1][mu.length];
			for (int i=0; i<mu.length; i++) {
				piarr[0][i] = pi[i];
			}
			EMStepPlotter.execute(imagename+iter, region, muarr, piarr, null, 1, mu.length, iter, 0, 0);
			 */
		}
		if (converged) {
			System.err.println("converged after "+iter+" iterations "+dfm.format(new Date()));
		} else {
			System.err.println("did not converge "+dfm.format(new Date()));
		}
	}

	public void runStage1Hybrid(int numiters) throws FileNotFoundException {
		float regionSize = region.getWidth();
		float backgroundProb = 1f/(regionSize*(regionSize+1f)/2f);
		System.err.println("background prob: "+backgroundProb);
		float localAlpha = 0;
		int localIncrement = X_INCREMENT;
		int iter = 0;
		int subiter = 250;
		if (numiters>-1) {
			iter = NO_ALPHA_ITER+ALPHA_ANNEAL_ITER;
			localAlpha = alpha;
			subiter = numiters;
			localIncrement = 1;
		}
		boolean converged = false;
		double likmax = -Float.MAX_VALUE;
		int maxcount = 0;
		int lastseen = 0;
		while (!converged) {
			totalLikelihood = 0;
			/*
			for (int j=0; j<mu.length; j++) {
				oldmu[j] = mu[j];
			}
			 */
			if (iter>=NO_ALPHA_ITER) {
				localAlpha = (((float)Math.min(iter-NO_ALPHA_ITER,ALPHA_ANNEAL_ITER)) / ((float)ALPHA_ANNEAL_ITER))*alpha;
			}
			System.err.println("alpha: "+localAlpha);
			//E Step
			for (int i=0; i<gamma.length; i++) {
				float norm = 0;
				Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(i);
				float[] tmpgamma;
				for (int j : gamma[i].keySet()) { //iterating over components
					if (j!=nullIndex) {
						tmpgamma = gamma[i].get(j);
						tmpgamma[0] = 0;
						tmpgamma[1] = 0;
						for (int k=0; k<chi.length; k++) {

							float tmpprob = dist.probability(currentPair, mu[j], null, k, nullIndex);
							/*
							if (j==12) {
								System.err.println(tmpprob+"\t"+pi[j]+"\t"+gamma2[j][k]);
							}
							 */
							tmpgamma[0] += rho*pi[j]*gamma2[j][k]*tmpprob;
							norm += tmpgamma[0];
							tmpprob = dist.probability(currentPair, mu[j], mu[j], k, k);
							/*
							if (j==12) {
								System.err.println(tmpprob+"\t"+pi[j]+"\t"+gamma2[j][k]);
							}
							 */
							tmpgamma[1] += (1f-rho)*pi[j]*gamma2[j][k]*tmpprob;
							norm += tmpgamma[1];
						}
					}
				}
				tmpgamma = gamma[i].get(nullIndex);
				tmpgamma[0] = pi[nullIndex]*backgroundProb;
				norm += tmpgamma[0];
				if (norm>0) {
					for (int j : gamma[i].keySet()) {
						if (j==nullIndex) {
							gamma[i].get(j)[0] /= norm;
						} else {
							tmpgamma = gamma[i].get(j);
							tmpgamma[0] /= norm;
							tmpgamma[1] /= norm;
						}
					}
				}
			}

			//M Step
			float[] N_pi = new float[mu.length];
			float[] N = new float[mu.length];
			//mu (and tau?)
			System.err.println(deadComponents.size()+" dead components");
			for (int j=0; j<mu.length; j++) {
				if (componentAliveStage1(j)) {
					Point currentmu = mu[j];
					int bestx = -1;
					float bestlikelihood = Float.NEGATIVE_INFINITY;
					Region currentmuregion = currentmu.expand(dist.radius);
					Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(currentmuregion);
					int leftminn = leftminmaxn.car();
					int leftmaxn = leftminmaxn.cdr();
					for (int x = Math.max(0,currentmu.getLocation()-2*dist.radius); x <= currentmu.getLocation()+2*dist.radius; x += localIncrement) {
						Point xpoint = new Point(g,mu[j].getChrom(), x);
						float likelihood = 0;
						for (int n = leftminn; n < leftmaxn; n++) {
							Pair<StrandedPoint,StrandedPoint> r_n = storage.getLeftPair(n);
							for (int k=0; k<chi.length; k++) {
								float tmpprob = dist.probability(r_n, xpoint, null, k, nullIndex);
								if (tmpprob==0) {
									tmpprob = Float.MIN_VALUE;
								}
								if (gamma[n].get(j)==null) {
									System.err.println(n+" "+r_n.car()+" "+j+" "+mu[j]);
								}
								likelihood += gamma[n].get(j)[0]*gamma2[j][k]*Math.log(tmpprob);
								tmpprob = dist.probability(r_n, xpoint, xpoint, k, k);
								if (tmpprob==0) {
									tmpprob = Float.MIN_VALUE;
								}
								likelihood += gamma[n].get(j)[1]*gamma2[j][k]*Math.log(tmpprob);
							}
						}
						if (likelihood>bestlikelihood) {
							bestx = x;
							bestlikelihood = likelihood;
						}
					}
					if (bestx!=-1) {
						totalLikelihood += bestlikelihood;
						Point newpoint = new Point(currentmu.getGenome(), currentmu.getChrom(), bestx);
						mu[j] = newpoint;
						updateHashGamma(newpoint,j);
					}
					for (int n = leftminn; n < leftmaxn; n++) {
						for (int k=0; k<2; k++) {
							N[j] += gamma[n].get(j)[k];

						}
						N_pi[j] += gamma[n].get(j)[0];
					}
				}
			}
			float N_piall = 0;
			float pinorm = 0;
			for (int j=0; j<N.length; j++) {
				N_piall += Math.max(N_pi[j]-localAlpha, 0f);
				pi[j] = Math.max(N[j]-localAlpha, 0f);
				pinorm += pi[j];
			}
			int numpairs = storage.getLeftN();
			pi[pi.length-1]= 0; 
			for (int j=0; j<numpairs; j++) {
				pi[pi.length-1] += gamma[j].get(nullIndex)[0]; 
			}
			rho = (N_piall+a) / (((float)pinorm)+a+b);
			pinorm += pi[pi.length-1];
			System.err.println("pinorm: "+pinorm);
			for (int j=0; j<pi.length; j++) {
				pi[j] /= pinorm;
			}
			//printPi();
			for (int i=0; i<rhoHistory.length-1; i++) {
				rhoHistory[i] = rhoHistory[i+1];
			}
			rhoHistory[rhoHistory.length-1] = rho;
			for (int i=0; i<likelihoodHistory.length-1; i++) {
				likelihoodHistory[i] = likelihoodHistory[i+1];
			}
			likelihoodHistory[likelihoodHistory.length-1] = totalLikelihood;
			likelihoodList.add(totalLikelihood);
			// check for convergence
			converged = true;
			if (iter==maxIters) {
				converged = false;
				break;
			}
			if (iter<NO_ALPHA_ITER+ALPHA_ANNEAL_ITER) {
				converged = false;
			} else {
				/*
				if (totalLikelihood>likmax) {
					lastseen = iter;
					likmax = totalLikelihood;
					maxcount = 0;
					converged = false;
				} else if (totalLikelihood==likmax) {
					lastseen = iter;
					maxcount++;
					if (maxcount<MINMAXCOUNT) {
						converged = false;
					}
				} else if (iter-lastseen>MINMAXCOUNT) {
					converged = false;
					likmax = localMax(likelihoodList, MINMAXCOUNT-2);
					lastseen = iter;
					maxcount = 0;
				} else {
					converged = false;
				}
				 */
				/*
				for (int i=3; i<likelihoodHistory.length; i += 2) {
					if (likelihoodHistory[i]!=likelihoodHistory[i-2]) {
						converged = false;
					}
				}
				 */

				if (subiter>0) {
					subiter--;
					converged = false;
				} else {
					converged = true;
				}

				/*
				float[] minMaxRho = minMax(rhoHistory);
				float[] minMaxLikelihood = minMax(likelihoodHistory);
				if (minMaxRho[1]-minMaxRho[0]>RHO_ZONE || minMaxLikelihood[1]-minMaxLikelihood[0]>LIKELIHOOD_ZONE) {
					converged = false;
				}
				 */
			}
			if (converged && localIncrement>1) {
				lastseen = iter;
				likmax = -Float.MAX_VALUE;
				maxcount = 0;
				converged = false;
				subiter = 250;
				localIncrement = Math.max(localIncrement/10, 1);
				System.err.println("Increment reduced to "+localIncrement);
			}
			iter++;
			System.err.println(iter+" iterations completed "+dfm.format(new Date()));
			dumpStage1();
			/*
			int[][] muarr = new int[1][mu.length];
			for (int i=0; i<mu.length; i++) {
				muarr[0][i] = mu[i].getLocation();
			}
			double[][] piarr = new double[1][mu.length];
			for (int i=0; i<mu.length; i++) {
				piarr[0][i] = pi[i];
			}
			EMStepPlotter.execute(imagename+iter, region, muarr, piarr, null, 1, mu.length, iter, 0, 0);
			 */
		}
		if (converged) {
			System.err.println("converged after "+iter+" iterations "+dfm.format(new Date()));
		} else {
			System.err.println("did not converge "+dfm.format(new Date()));
		}
	}

	private void printPi() {
		for (int i=0; i<pi.length; i++) {
			System.err.print(pi[i]+"\t");
		}
		System.err.println();
	}

	public double localMax(List<Double> dlist, int range) {
		double max = -Double.MAX_VALUE;
		for (int i=dlist.size()-range; i<dlist.size(); i++) {
			if (dlist.get(i)>max) {
				max = dlist.get(i);
			}
		}
		return max;
	}

	private float[] minMax(float[] arr) {
		float[] tor = {Float.MAX_VALUE, -Float.MAX_VALUE};
		for (int i=0; i<arr.length; i++) {
			if (arr[i]<tor[0]) tor[0] = arr[i];
			if (arr[i]>tor[1]) tor[1] = arr[i];
		}
		return tor;
	}

	private double[] minMax(double[] arr) {
		double[] tor = {Double.MAX_VALUE, -Double.MAX_VALUE};
		for (int i=0; i<arr.length; i++) {
			if (arr[i]<tor[0]) tor[0] = arr[i];
			if (arr[i]>tor[1]) tor[1] = arr[i];
		}
		return tor;
	}

	/*
	public void run() throws FileNotFoundException {
		//variables should be initialized
		float localAlpha = 0;
		float localBeta = 0;
		int iter = 0;
		boolean converged = false;
		while (!converged) {
			for (int j=0; j<mu.length; j++) {
				oldmu[j] = mu[j];
			}
			if (iter>=NO_ALPHA_ITER) {
				localAlpha = (((float)Math.min(iter-NO_ALPHA_ITER,ALPHA_ANNEAL_ITER)) / ((float)ALPHA_ANNEAL_ITER))*alpha;
				localBeta = (((float)Math.min(iter-NO_ALPHA_ITER,ALPHA_ANNEAL_ITER)) / ((float)ALPHA_ANNEAL_ITER))*beta;
			}
			//E Step
			for (int i=0; i<gamma.length; i++) {
				float norm = 0;
				Pair<StrandedPoint,StrandedPoint> currentPair = storage.getLeftPair(i);
				for (int j=0; j<gamma[i].length; j++) {
					for (int k=j; k<gamma[i][j].length; k++) {
						if (k==nullIndex) {
							gamma[i][j][k] = rho*pi[j]*dist.probability(currentPair,mu[j],null,tau[j],nullIndex);
						} else {
							gamma[i][j][k] = (1-rho)*psi[j][k]*dist.probability(currentPair,mu[j],mu[k],tau[j],tau[k]);
						}
						norm += gamma[i][j][k];
					}
				}

				for (int j=0; j<gamma[i].length; j++) {
					for (int k=j; k<gamma[i][j].length; k++) {
						gamma[i][j][k] /= norm;
					}
				}
			}
			//M Step
			float N_pi = 0;
			float N_psi = 0;
			float[][] N = new float[mu.length][mu.length+1];
			//mu (and tau?)
			for (int j=0; j<mu.length; j++) {
				if (componentAlive(j)) {
					Point currentmu = mu[j];
					int bestx = -1;
					float bestlikelihood = Float.NEGATIVE_INFINITY;
					Region currentmuregion = currentmu.expand(dist.radius);
					Pair<Integer,Integer> leftminmaxn = storage.getLeftMinMaxn(currentmuregion);
					int leftminn = leftminmaxn.car();
					int leftmaxn = leftminmaxn.cdr();
					Pair<Integer,Integer> rightminmaxn = storage.getRightMinMaxn(currentmuregion);
					int rightminn = rightminmaxn.car();
					int rightmaxn = rightminmaxn.cdr();
					for (int x = Math.max(0,currentmu.getLocation()-2*dist.radius); x <= currentmu.getLocation()+2*dist.radius; x++) {
						float likelihood = 0;
						for (int n = leftminn; n < leftmaxn; n++) {
							Pair<StrandedPoint,StrandedPoint> r_n = storage.getLeftPair(n);
							likelihood += gamma[n][j][nullIndex]*Math.log(dist.probability(r_n, mu[j], null, tau[j], nullIndex));
							for (int k=j; k<gamma[n].length; k++) {
								float tmpgamma = gamma[n][j][k];
								likelihood += tmpgamma*Math.log(dist.probability(r_n, mu[j], mu[k], tau[j], tau[k]));
							}
						}
						for (int n = rightminn; n < rightmaxn; n++) {
							Pair<StrandedPoint,StrandedPoint> r_n = storage.getRightPair(n);
							for (int k=0; k<j; k++) {
								float tmpgamma = gamma[n][k][j];
								likelihood += tmpgamma*Math.log(dist.probability(r_n, mu[k], mu[j], tau[k], tau[j]));
							}
						}
						// add log prior to likelihood
						if (likelihood>bestlikelihood) {
							bestx = x;
							bestlikelihood = likelihood;
						}
					}
					mu[j] = new Point(currentmu.getGenome(), currentmu.getChrom(), bestx);
					// M step for Ns
					for (int n = leftminn; n < leftmaxn; n++) {
						for (int k=j; k<N[j].length; k++) {
							N[j][k] += gamma[n][j][k];
						}
					}
					for (int n = rightminn; n < rightmaxn; n++) {
						for (int k=0; k<j; k++) {
							N[k][j] += gamma[n][k][j];
						}
					}
				}
			}

			// normalize Ns to get pi, psi, and rho
			for (int j=0; j<N.length; j++) {
				for (int k=j; k<N[j].length-1; k++) {
					N[j][k] = Math.max(N[j][k]-localBeta, 0);
					N_psi += N[j][k];
				}
				N[j][N[j].length-1] = Math.max(N[j][N[j].length-1]-localAlpha, 0);
				N_pi += N[j][N[j].length-1];
			}
			for (int j=0; j<N.length; j++) {
				pi[j] = N[j][N[j].length-1] / N_pi;
				for (int k=j; k<N[j].length-1; k++) {
					psi[j][k] = N[j][k] / N_psi;
				}
			}
			rho = (N_pi+a) / (((float)storage.getLeftN())+a+b);

			// check for convergence
			converged = true;
			if (iter==maxIters) {
				converged = false;
				break;
			}
			if (iter<NO_ALPHA_ITER+ALPHA_ANNEAL_ITER) {
				converged = false;
			} else {
				for (int j=0; j<mu.length; j++) {
					if (!mu[j].equals(oldmu[j])) {
						converged = false;
						break;
					}
				}
			}
			iter++;
			System.err.println(iter+" iterations completed "+dfm.format(new Date()));
			dump();
		}
		if (converged) {
			System.err.println("converged after "+iter+" iterations "+dfm.format(new Date()));
		} else {
			System.err.println("did not converge "+dfm.format(new Date()));
		}
	}
	 */

	private boolean componentAliveStage1(int i) {
		if (deadComponents.contains(i)) {
			return false;
		}
		if (pi[i]>0) {
			return true;
		}
		deadComponents.add(i);
		return false;
	}

	/*
	private boolean componentAlive(int i) {
		if (deadComponents.contains(i)) {
			return false;
		}
		if (pi[i]>0) {
			return true;
		}
		for (int j=0; j<i; j++) {
			if (psi[j][i]>0) {
				return true;
			}
		}
		for (int j=i; j<pi.length; j++) {
			if (psi[i][j]>0) {
				return true;
			}
		}
		deadComponents.add(i);
		return false;
	}
	 */

}
