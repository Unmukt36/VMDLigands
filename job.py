import argparse, os
#---------------------------------------------------------------------------------------#
if __name__ == "__main__":
# Add all arguments
#---------------------------------------------------------------------------------------#
	parser = argparse.ArgumentParser(	description="Generate movie from XYZ file",
																		)

	parser.add_argument(	"xyzName", action="store", default=None,
											help="Relative/Absolute path of the property file"	)

	parser.add_argument(	"-c", "--color", nargs=2, action="store", default=None,
											help="Relative/Absolute path, column number of the property file"	)

	parser.add_argument(	"-n", "--frames", action="store", type=int, default=None,
											help="Reduce the trajectory file to nFrames"	)

	parser.add_argument(	"-s", "--submit", action="store_true", default=False,
											help="Submit the job"	)

	parser.add_argument(	"-f", "--form", action="store", default="mp4",
											help="Specify the format of the output movie"	)

	parser.add_argument(	"-p", "--plot", action="store_true", default=False,
											help="Simultaneously show evolution of (a) propert(y)ies with the movie"	)

	args = parser.parse_args()
#---------------------------------------------------------------------------------------#
# Submit
#---------------------------------------------------------------------------------------#
	# Pass Namespace to workhorse in the form of a dict so as to be safely evaluated
	subCMD = "/usr/bin/python workMovie.py \"%s\"\n" % (str(vars(args)))
	if (args.submit):
		subName = 'run.sh'
#---------------------------------------------------------------------------------------#
# SLURM scheduling batch file
#---------------------------------------------------------------------------------------#
		with open(subName, 'w') as shFile:
			shFile.write('#!/bin/bash\n')
			shFile.write('#SBATCH -J "%s"\n' % "Render")
			shFile.write('#SBATCH -p fe13\n')
			shFile.write('#SBATCH --ntasks=%d\n' % (1))
			#shFile.write("#SBATCH --exclusive\n")
			# Ensures no time limit is set on the run
			shFile.write('#SBATCH --time=1-00:00:00')
			shFile.write('\n')
			## Copy data to a local tmp space on the compute node to reduce I/O
			shFile.write('export MYTMP=/tmp/$USER/png_vmd_$SLURM_JOB_ID\n')
			shFile.write('/usr/bin/mkdir -p $MYTMP || exit $?\n')
			shFile.write('echo \"Copy data to /tmp on the compute node.\"\n')
			shFile.write('module load vmd\n')	
			#shFile.write('cp -rp $SLURM_SUBMIT_DIR/mydatadir $MYTMP || exit $?\n')
			shFile.write('\n')
			shFile.write('#SBATCH -e $MYTMP/err.job\n')
			shFile.write('#SBATCH -o $MYTMP/out.job\n')
			shFile.write('\n')

			shFile.write(subCMD)
			shFile.write('\n')

			shFile.write('echo "Ended at `date` on `hostname`."\n')
			shFile.write('echo "Copy data back to $HOME."\n') 
			#shFile.write('/usr/bin/mkdir -p $SLURM_SUBMIT_DIR/newdatadir || exit $?\n')
			shFile.write('cp -rp $MYTMP $SLURM_SUBMIT_DIR || exit $?\n')
			## remove your data from the compute node /tmp space
			shFile.write('rm -rf $MYTMP\n')
			shFile.write('\n')
			shFile.write('exit 0')
		#---------------------------------------------------------------------------------------#
		subCMD = "sbatch %s" % subName

	print subCMD
	os.system(subCMD)
	if (args.submit == True):
		if subName in vars(__builtins__):
			if os.path.exists(subName):
				os.system("rm -rf %s" % subName)
#---------------------------------------------------------------------------------------#

