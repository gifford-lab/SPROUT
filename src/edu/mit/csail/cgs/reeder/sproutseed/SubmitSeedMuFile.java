package edu.mit.csail.cgs.reeder.sproutseed;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.reeder.sprout.SproutUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class SubmitSeedMuFile {

	/**
	 * @param args
	 * @throws NotFoundException 
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws NotFoundException, FileNotFoundException {
		Genome g = SproutUtils.parseGenome(args);
		List<String> chromlist = g.getChromList();
		String genome = Args.parseString(args, "genome", "");
		float rho = Args.parseFloat(args, "rho", 0);
		float alpha = Args.parseFloat(args, "alpha", 0);
		float beta = Args.parseFloat(args, "beta", 0);
		float a = Args.parseFloat(args, "a", 0);
		float b = Args.parseFloat(args, "b", 0);
		String readfile = Args.parseString(args, "readfile", "");
		Region region;
		String dumpfile = Args.parseString(args, "dumpfile", "");
		String outfile = Args.parseString(args, "outfile", "");
		String directory = Args.parseString(args, "directory", "");
		int stage = Args.parseInteger(args, "stage", 0);
		int maxiters = Args.parseInteger(args, "maxiters", 0);
		String mutaubase = Args.parseString(args, "mutaubase", "");
		int mutaunum = Args.parseInteger(args, "mutaunum", 0);
		String wd = Args.parseString(args, "wd", "/");
		String submitfile = Args.parseString(args, "submitfile", "");
		PrintStream out = new PrintStream(submitfile);
		out.println("#!/bin/bash");
		out.println();
		for (int i=0; i<mutaunum; i++) {
			String mutaufile = mutaubase+i+".txt";
			List<String> command = new ArrayList<String>();
			command.add("/usr/bin/qsub");
			command.add("-l");
			command.add("mem_free=2G");
			command.add("-wd");
			command.add("/cluster/ccr/sprout_seed"+wd);
			command.add("-q");
			command.add("batch");
			command.add("-v");
			command.add("'genome="+genome+",rho="+rho+",alpha="+alpha+",beta="+beta+",a="+a+",b="+b+",readfile="+
					readfile+",dumpfile="+dumpfile+i+",outfile="+outfile+i+".txt,directory="+directory+
					",stage="+stage+",maxiters="+maxiters+",mutaufile="+mutaufile+"'");
			command.add("-o");
			command.add("output"+outfile+i+".txt");
			command.add("-e");
			command.add("error"+outfile+i+".txt");
			command.add("../runseedmufile.sh");
			for (int k=0; k<command.size(); k++) {
				out.print(command.get(k)+" ");
			}
			out.println();
		}
		out.flush();
		out.close();
	}

}
