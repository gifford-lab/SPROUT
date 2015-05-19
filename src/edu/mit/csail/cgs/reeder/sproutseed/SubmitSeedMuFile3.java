package edu.mit.csail.cgs.reeder.sproutseed;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.reeder.sprout.SproutUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class SubmitSeedMuFile3 {

	/**
	 * @param args
	 * @throws FileNotFoundException 
	 * @throws NotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException, NotFoundException {
		Genome g = SproutUtils.parseGenome(args);
		List<String> chromlist = g.getChromList();
		String genome = Args.parseString(args, "genome", "");
		float rho = Args.parseFloat(args, "rho", 0);
		float alpha = Args.parseFloat(args, "alpha", 0);
		float beta = Args.parseFloat(args, "beta", 0);
		float a = Args.parseFloat(args, "a", 0);
		float b = Args.parseFloat(args, "b", 0);
		String readfile = Args.parseString(args, "readfile", "");
		String region;
		String dumpfile = Args.parseString(args, "dumpfile", "");
		String outfile = Args.parseString(args, "outfile", "");
		String directory = Args.parseString(args, "directory", "");
		String genericdirectory = Args.parseString(args, "genericdirectory", directory);
		int stage = Args.parseInteger(args, "stage", 0);
		int maxiters = Args.parseInteger(args, "maxiters", 0);
		String stage2file = Args.parseString(args, "stage2file", "");
		String eventout = Args.parseString(args, "eventout", "");
		String interactionout = Args.parseString(args, "interactionout", "");
		String wd = Args.parseString(args, "wd", "/");
		String submitfile = Args.parseString(args, "submitfile", "");
		PrintStream out = new PrintStream(submitfile);
		out.println("#!/bin/bash");
		out.println();
		for (String c : chromlist) {
			region = SproutUtils.chromRegion(c, g).toString();
			List<String> command = new ArrayList<String>();
			command.add("/usr/bin/qsub");
			command.add("-l");
			command.add("mem_free=4G");
			command.add("-wd");
			command.add("/cluster/ccr/sprout_seed"+wd);
			command.add("-q");
			command.add("batch");
			command.add("-v");
			command.add("'genome="+genome+",rho="+rho+",alpha="+alpha+",beta="+beta+",a="+a+",b="+b+",readfile="+
					readfile+",region="+region+",dumpfile="+dumpfile+c+",outfile="+outfile+c+".txt,directory="+directory+
					",genericdirectory="+genericdirectory+",stage="+stage+",maxiters="+maxiters+
					",eventout="+eventout+c+"out.txt,interactionout="+interactionout+c+"out.txt,stage2file="+stage2file+c+".txt'");
			command.add("-o");
			command.add("output"+outfile+c+".txt");
			command.add("-e");
			command.add("error"+outfile+c+".txt");
			command.add("../runseedmufile3.sh");
			for (int k=0; k<command.size(); k++) {
				out.print(command.get(k)+" ");
			}
			out.println();
		}
		out.flush();
		out.close();
	}

}
