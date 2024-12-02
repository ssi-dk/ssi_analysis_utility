#!/usr/bin/env python

import sys;
import os;
import re;

VERSION = "1.0.0";

def openOrDie(filename, mode):
	try:
		fileP = open(filename, mode);
	except:
		sys.stderr.write("No such file or directory:\t%s\n", filename);
		sys.exit(1);
	return fileP;

def helpMessage():	
	sys.stderr.write("LRE-Finder.py calls kma and getGene to find relevant postions in the kma-mapping.\n");
	sys.stderr.write("The options are identical to that of kma.\n");
	sys.exit(0);

path = sys.argv[0][:sys.argv[0].rfind("/") + 1];
if(path == "./"):
	path = "";

kma_args = [];
try:
	infile = open(path + "kma", "rb");
	infile.close();
	kma_args = [path + "./kma"];
except:
	kma_args = ["kma"];

getGene_args = [];
try:
	infile = open(path + "getGene", "rb");
	infile.close();
	getGene_args = [path + "./getGene"];
except:
	getGene_args = ["getGene"];

outputfilename = "";
t_db = "";
argc = len(sys.argv);
args = 1;
while(args < argc):
	if(sys.argv[args] == "-t_db"):
		args += 1;
		if(args < argc):
			t_db = sys.argv[args];
			getGene_args.append("-pos");
			getGene_args.append(t_db + ".genepos");
			kma_args.append("-t_db");
			kma_args.append(t_db);
	elif(sys.argv[args] == "-o"):
		args += 1;
		if(args < argc):
			kma_args.append("-o");
			kma_args.append(sys.argv[args]);
			getGene_args.append("-i");
			getGene_args.append(sys.argv[args]);
			outputfilename = sys.argv[args];
	elif(sys.argv[args] == "-h"):
		helpMessage();
	elif(sys.argv[args] == "-v"):
		sys.stderr.write("%s\n" %(VERSION));
		sys.exit(0);
	else:
		kma_args.append(sys.argv[args]);
	args += 1;

if(len(t_db) == 0):
	sys.stderr.write("Need DB file.\n");
	helpMessage();
elif(len(outputfilename) == 0):
	sys.stderr.write("kma needs an outputfile under option \"-o\"\n");
	helpMessage();
if("-matrix" not in kma_args):
	kma_args.append("-matrix");

#call kma
os.system(" ".join(kma_args));

#get positions
outputfile = open(outputfilename + ".pos", "w");
outputfile.write("#Template\tpos\tref\tA\tC\tG\tT\tN\t-\tA[%]\tC[%]\tG[%]\tT[%]\tN[%]\t-[%]\n");
outputfile.close();
os.system("%s >> %s" %(" ".join(getGene_args), outputfilename + ".pos"));

#Make HTML output
# Get info
potentialSpecies = {"Unknown" : "Unknown"};
TargetMuts = {};
targetMuts = {};
Species = "Unknown";
infile = open(t_db + ".mutdist", "r");
for line in infile:
	line = line.rstrip();
	if(line.startswith('>')):
		TargetMuts[Species] = targetMuts;
		spec = line[1:].split("\t");
		Species = spec[0];
		if(len(spec) == 2):
			potentialSpecies[Species] = spec[1];
		else:
			potentialSpecies[Species] = Species;
		targetMuts = {};
	else:
		spec = line.split("\t");
		
		m = re.search('^\w(\d+)(\w+)$', spec[0]);
		if(m):
			mutNuc = m.group(2);
			targetPos = int(m.group(1));
		else:
			sys.stderr.write("Wrong format of mutdist.\n");
			sys.exit(1);
		#mutNuc = spec[0][-1];
		#targetPos = int(spec[0][1:-1]);
		
		targetMuts[targetPos] = [mutNuc, float(spec[1]) * 100];
TargetMuts[Species] = targetMuts;
infile.close();

bases = "ACGTN-"
cutOff = 10;
Species = "Unknown";
Genes = [["Unknown", 0, 0]];
Muts = [];
info = [];
line = "";
targetPos = 0;
wildRatio = 0.0;
resRatio = 0.0;
pheno = '';

geneFile = open(outputfilename + ".res", "r");
line = geneFile.readline();
for line in geneFile:
	line = line.rstrip();
	info = line.split("\t");
	if(info[0] in potentialSpecies):
		#Check which species prediction is the strongest
		if(float(Genes[0][2]) < float(info[8])):
			Species = info[0];
			Genes[0] = [info[0], float(info[4]), float(info[8])];
	else:
		Genes.append([info[0], float(info[4]), float(info[8])]);
geneFile.close();

mutFile = open(outputfilename + ".pos", "r");
line = mutFile.readline();
targetMuts = TargetMuts[Species];
for line in mutFile:
	line = line.rstrip();
	info = line.split("\t");
	if(info[0] == Species):
		targetPos = int(info[1]);
		wildNuc = info[2];
		wildRatio = float(info[9 + bases.find(wildNuc)]);
		tmp = targetMuts.get(targetPos, []);
		if(len(tmp) == 0):
			sys.stderr.write("Database not updated properly.\n");
			sys.exit(1);
		mutNucs = tmp[0];
		cutOff = tmp[1];
		for mutNuc in mutNucs:
			resRatio = float(info[9 + bases.find(mutNuc)]);
			if(cutOff < resRatio):
				pheno = 'R';
			else:
				pheno = 'S';
			Muts.append(["%s%d%s" %(wildNuc, targetPos, mutNuc), wildRatio, resRatio, pheno]);
mutFile.close();

# Output final HTML
#Species
sys.stdout.write("<b>Species identified:</b>\t<table style=\"width:100%\">\n");
sys.stdout.write("<tr bgcolor=#B22222><th>%s</th></tr><br>\n" %(potentialSpecies[Species]));
sys.stdout.write("</table><br>\n");

#Genes
sys.stdout.write("<b>Genes Identified:</b><br>\n<table style=\"width:100%\">\n");
sys.stdout.write("<tr bgcolor=#B22222><th>Gene</th><th>Template_identity</th><th>Depth</th></tr>\n");
for gene in Genes:
	sys.stdout.write("<tr><td>%s</td><td>%.1f</td><td>%.1f</td><tr>\n" %(gene[0], gene[1], gene[2]));
sys.stdout.write("</table><br>\n");

#Muts
sys.stdout.write("<b>Identified mutations in 23s:</b><br>\n<table style=\"width:100%\">\n");
sys.stdout.write("<tr bgcolor=#B22222><th>Position in reference</th><th>Wild type ratio [%]</th><th>Mutant type ratio [%]</th><th>Predicted phenotype</th></tr>\n");
for mut in Muts:
	sys.stdout.write("<tr><td>%s</td><td>%.1f</td><td>%.1f</td><td>%s</td><tr>\n" %(mut[0], mut[1], mut[2], mut[3]));
sys.stdout.write("</table></tr><br>\n");

#Warning at low depth
if(float(Genes[0][2]) < 100):
	sys.stdout.write("<font color=\"#B22222\" size=\"6\"><b>WARNING</b><br>\n");
	sys.stdout.write("Please note that you have submitted sequencing data with a 23S depth lower than 100.<br>This may influence the interpretation of the data.<br>Please ensure you are uploading raw Illumina sequencing data (not draft genome assemblies) with appropriate depth of coverage.<br>LRE-Finder has only been validated on Illumina raw reads in fastq format.</font><br><br>\n");

