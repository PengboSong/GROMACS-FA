#!/usr/bin/env python3

import os, sys, subprocess, shlex

def main(argv):
	argc = len(argv)
	pwd = os.getcwd()

	gmxbinary = "gmx"
	box = "cubic"
	edged = 1.0
	prof = "pro.gro"
	topolf = "topol.top"
	indexf = "index.ndx"
	genionf = "ions.mdp"
	solvatef = ""
	pion = "NA"
	nion = "CL"
	conc = 0.1
	pn = 0
	nn = 0

	help_flag = False
		
	for i in range(argc):
		if argv[i] == "-dir" and os.path.isdir(argv[i+1]):
			os.chdir(argv[i+1])
		elif argv[i] == "-h":
			help_flag = True
	
	if help_flag is True:
		print("""-h	Print help information.
-dir <path>	Change work directory to <path>.
-f <path>	Specify frame file to <path>.
-top <path>	Specify topology file to <path>.
-ndx <path>	Specify index file to <path>.
-gmx <binary>	Specfiy GROMACS binary name to <binary>.
*** Solvation ***
-box <box-type>	Set solvate box type to <box-type>.
-edged <float>	Set solvate box edge to <float>.
-solvate <path>	Specify solvation template file to <path>.
*** Generate ions ***
-genion <path>	Specify parameters file for generating ions to <path>.
-pion <ion>	Set positive ion name to <ion>.
-nion <ion>	Set negative ion name to <ion>.
-pn <int>	Set positive ion number to <int>.
-nn <int>	Set negative ion number to <int>.
-conc <float>	Set salt concentration to <float>.""")
		return 0
	
	for i in range(argc):
		if argv[i] == "-gmx":
			gmxbinary = argv[i+1]
		elif argv[i] == "-box":
			box = argv[i+1]
		elif argv[i] == "-edged":
			edged = float(argv[i+1])
		elif argv[i] == "-f" and os.path.isfile(argv[i+1]) and os.path.splitext(argv[i+1])[1] == ".gro":
			prof = argv[i+1]
		elif argv[i] == "-top" and os.path.isfile(argv[i+1]) and os.path.splitext(argv[i+1])[1] == ".top":
			topolf = argv[i+1]
		elif argv[i] == "-ndx" and os.path.isfile(argv[i+1]) and os.path.splitext(argv[i+1])[1] == ".ndx":
			indexf = argv[i+1]
		elif argv[i] == "-genion" and os.path.isfile(argv[i+1]) and os.path.splitext(argv[i+1])[1] == ".mdp":
			genionf = argv[i+1]
		elif argv[i] == "-solvate" and os.path.splitext(argv[i+1])[1] == ".gro":
			solvatef = argv[i+1]
		elif argv[i] == "-pion":
			pion = argv[i+1]
		elif argv[i] == "-nion":
			nion = argv[i+1]
		elif argv[i] == "-pn":
			pn = int(argv[i+1])
		elif argv[i] == "-nn":
			nn = int(argv[i+1])
		elif argv[i] == "-conc":
			conc = float(argv[i+1])
	
	prefix = os.path.splitext(prof)[0]
	# Resize box
	inf = prof
	outf = prefix + "_box.gro"
	f = open(prefix + "_pre.log", 'wb')
	ferr = open(prefix + "_pre_error.log", 'wb')
	if os.path.isfile(inf):
		print("Start resizing box.")
		subprocess.check_call(shlex.split("{0} editconf -f {1} -o {2} -bt {3} -c -d {4:.2f}".format(gmxbinary, inf, outf, box, edged)), stdout=f, stderr=ferr)
		if os.path.isfile(outf):
			print("Resizing box succeed.")
		else:
			sys.exit(2)
	else:
		sys.exit(1)
	
	# Solvation
	inf = outf
	outf = prefix + "_solv.gro"
	if os.path.isfile(inf):
		print("Start solvation.")
		solvatecmd = "{0} solvate -cp {1} -o {2} -p {3}".format(gmxbinary, inf, outf, topolf)
		if solvatef:
			solvatecmd += " -cs {0}".format(solvatef)
		subprocess.check_call(shlex.split(solvatecmd), stdout=f, stderr=ferr)
		if os.path.isfile(outf):
			print("Solvation succeed.")
		else:
			sys.exit(4)
	else:
		sys.exit(3)
	
	# Generate ions
	inf = outf
	outf = prefix + "_ions.gro"
	if os.path.isfile(inf):
		print("Start generating ions.")
		genionbinaryf = os.path.splitext(genionf)[0] + ".tpr"
		subprocess.check_call(shlex.split("{0} grompp -c {1} -r {1} -f {2} -p {3} -o {4}".format(gmxbinary, inf, genionf, topolf, genionbinaryf)), stdout=f, stderr=ferr)
		if os.path.isfile(genionbinaryf):
			with open("in.txt", 'w') as comm:
				comm.write("SOL\n")
			incommand = open("in.txt", 'r')
			if pn == 0 and nn == 0:
				subprocess.check_call(shlex.split("{0} genion -s {1} -o {2} -p {3} -pname {4} -nname {5} -neutral -conc {6}".format(gmxbinary, genionbinaryf, outf, topolf, pion, nion, conc)), stdin=incommand, stdout=f, stderr=ferr)
			else:
				subprocess.check_call(shlex.split("{0} genion -s {1} -o {2} -p {3} -pname {4} -np {5} -nname {6} -nn {7}".format(gmxbinary, genionbinaryf, outf, topolf, pion, pn, nion, nn)), stdin=incommand, stdout=f, stderr=ferr)
			incommand.close()
			os.remove("in.txt")
			if os.path.isfile(outf):
				print("Generating ions succeed.")
			else:
				sys.exit(7)
		else:
			sys.exit(6)
	else:
		sys.exit(5)
	
	# Make index
	inf = outf
	outf = indexf
	if os.path.isfile(inf):
		print("Start making system index.")
		with open("in.txt", 'w') as comm:
			comm.write("q\n")
		incommand = open("in.txt")
		subprocess.check_call(shlex.split("{0} make_ndx -f {1} -o {2}".format(gmxbinary, inf, outf)), stdin=incommand, stdout=f, stderr=ferr)
		incommand.close()
		os.remove("in.txt")
		if os.path.isfile(outf):
			print("Making index succeed.")
		else:
			sys.exit(9)
	else:
		sys.exit(8)
	
	f.close()
	ferr.close()
	return 0

if __name__ == '__main__':
	sys.exit(main(sys.argv))
