import argparse, os, sys, ast, shutil, tempfile
import parser as ps
import numpy as np
#---------------------------------------------------------------------------------------#
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rc('text', usetex=True)
from matplotlib.ticker import MaxNLocator
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 22
#font = {'family':'sans-serif', 'serif':'Helvetica Neue', 'size':22}
#matplotlib.rc('font', **font)
matplotlib.rcParams.update({'figure.autolayout': True})
#----------------------------------------------------------------------------------#
# Generate epsilon ramp map (for temporal only) based on run data in out.* file
#----------------------------------------------------------------------------------#
def _eps(timeList, runData, e_max, Nfreq, loc=-2, e_min=0.0):
	loc = (loc+len(runData))%len(runData)
	rd = np.cumsum(runData)*Nfreq
	# Time list contains actual value of the timestep
	# tLoc contains the interval (for each time value in timeList) in which the data belongs
	tLoc = np.searchsorted(rd, timeList)
	epsArr = np.ones(tLoc.shape)*e_min
	for ti, tVal in enumerate(timeList):
		if tLoc[ti] == loc:
			epsArr[ti] += (e_max-e_min)*(tVal-rd[loc-1])/float(runData[loc]*Nfreq)
		elif tLoc[ti] > loc:
			epsArr[ti] = e_max
	return epsArr
#----------------------------------------------------------------------------------#
def _timeMatch(xyzTime, propTime):
#----------------------------------------------------------------------------------#
	# fMatch is the all of the matching/commom timesteps between XYZ and COM files
	# Ind arrays are the indices where these matches exist in the respective arrays
	fMatch = np.intersect1d(propTime, xyzTime)
	# Using ravel is fine because time values are monotonically increasing i.e. one match per "where
	comInd = np.array( [np.where(x==propTime) for x in fMatch] ).ravel()
	xyzInd = np.array( [np.where(x==xyzTime) for x in fMatch] ).ravel()
	# Make sure that the data for all XYZ snapshots is present in the COM file
	assert len(xyzTime) == len(xyzTime[xyzInd])
	#return common indices
	return comInd
#----------------------------------------------------------------------------------#
def _plot(inFile, xyzTime, seq, tmpDir):
#----------------------------------------------------------------------------------#
	# Turn-on interactive plotting
	plt.ion()
	# Set-up the canvas
	figg, axg = plt.subplots(figsize=(8,5))
	axg.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useMathText=True))
	# Four axes: twin the x-axis twice to make independent y-axes.
	axes = [axg, axg.twinx(), axg.twinx(), axg.twinx()]
	# Make some space on the right side for the extra y-axis.
	figg.subplots_adjust(right=0.75)
	# Move the last y-axis spine over to the right by 20% of the width of the axes
	axes[2].spines['right'].set_position(('axes', 1.25))
	# To make the border of the right-most axis visible, we need to turn the frame
	# on. This hides the other plots, however, so we need to turn its fill off.
	axes[2].set_frame_on(True)
	axes[2].patch.set_visible(False)

	axes[3].spines['right'].set_position(('axes', 1.5))
	#axes[3].tick_params(axis="y",direction="in", left="off",labelleft="on")
	axes[3].set_frame_on(True)
	axes[3].patch.set_visible(False)
	lines = []
	propList = []
#----------------------------------------------------------------------------------#
	markevery = 1
	colors = ['Green', 'Red', 'Blue', 'Black']
	alabel = (r"$\rm\alpha_{m} (deg)$", r"$\rm|H| (\rm\sigma)$", r"$\rmn_{bonds}$", r"$\rm\epsilon_{NP,NP}$")
	mark = ("o", "x", "v", "o")
#----------------------------------------------------------------------------------#
# \alpha_{min}
#----------------------------------------------------------------------------------#
	col = 7
	propName = ps.flipPath(inFile, 'or.', locType="Front", newLoc="Front")
	propVec = ps.MolVec(propName)
	propTime = propVec.timeList()
	# Find the common indices corresponding to time
	comInd = _timeMatch(xyzTime, propTime)
	propVec.headMask = propVec.headLoc[comInd]
	prop = propVec._genProp(col).mean(axis=1)[:,0]
	propList.append(prop)
	# Match with XYZ possibly
	line0, = axes[0].plot(xyzTime, prop, color=colors[0], linestyle="None", marker=mark[0], markevery=markevery, label=alabel[0], alpha=0.8)
	lines.append(line0)
#----------------------------------------------------------------------------------#
# y_{abs}
#----------------------------------------------------------------------------------#
	col = 12
	# Same file
	prop = propVec._genProp(col)
	pabs = (np.absolute(prop)).mean(axis=1)[:,0]
	propList.append(pabs)
	# Match with XYZ possibly
	line0, = axes[1].plot(xyzTime, pabs, color=colors[1], linestyle="None", marker=mark[1], markevery=markevery, label=alabel[1], alpha=0.8)
	lines.append(line0)
#----------------------------------------------------------------------------------#
# n_{bonds}
#----------------------------------------------------------------------------------#
	#propName = ps.flipPath(inFile, 'nb.', locType="Front", newLoc="Front")
	#data = ps.quickRead(propName)
	#propTime = data[:,0]
	# Find the common indices corresponding to time
	#comInd = _timeMatch(xyzTime[:-1], propTime)
	# Multiplied by 2 so as to reflect the actual number of bonds/molecule
	#prop = data[comInd,1]*2.0
	#propList.append(prop)
	#line0, = axes[2].plot(xyzTime[:-1], prop, color=colors[2], linestyle="None", marker=mark[2], markevery=markevery, label=alabel[2], alpha=0.8)
	lines.append(line0)
#----------------------------------------------------------------------------------#
# \epsilon_{NP,NP}
#----------------------------------------------------------------------------------#
	Nfreq = 1000
	eps = 0.6
	propName = ps.flipPath(inFile, 'out.', locType="Front", newLoc="Front")
	scaledTime, bDims, runData = ps.get_box(propName, xyzTime)
	epsArr = _eps(xyzTime, runData, eps, Nfreq, loc=-2, e_min=0)
	print epsArr, len(epsArr)
	#propList.append(epsArr)
	lineEPS = axes[3].fill_between(xyzTime, epsArr, color='k', alpha=0.6, label=alabel[3], zorder=1)
	#lines.append(line0)
	#axes[3].plot(xyzTime, epsArr, color=colors[3], label=alabel[3], alpha=1.0, zorder=1)
	#axes[3].legend(loc='best', frameon=False, prop={'size':18})
	axes[3].set_ylim([0.0, eps+0.001])
#----------------------------------------------------------------------------------#
# Final touches
#----------------------------------------------------------------------------------#
	axes[0].set_xlim([min(xyzTime), max(xyzTime)])
	for ai, ax in enumerate(axes):
		ax.set_ylabel(alabel[ai], color=colors[ai])
		ax.tick_params(axis='y', colors=colors[ai])
		ax.spines['right'].set_color(colors[ai])
	axes[0].yaxis.set_major_locator(MaxNLocator(prune='lower'))
#----------------------------------------------------------------------------------#
	axes[0].set_xlabel(r"Timesteps")
	for ah in np.cumsum(runData)*Nfreq:
		plt.axvline(x=ah, color='k', linestyle='--', alpha=1.0, zorder=0)
#----------------------------------------------------------------------------------#
	# Turn-on interactive plotting
	nd = int(np.log10(len(xyzTime))) + 2
	pSeq = os.path.join(tmpDir, "Plot-%s.png")
	cSeq = os.path.join(tmpDir, "Comb-%s.png")
	vres = 1024
	pres = 1200
	for ii, itime in enumerate(xyzTime):
		for pi, pl in enumerate(propList):
			print xyzTime[:ii], propList[pi][:ii]
			lines[pi].set_data(xyzTime[:ii], propList[pi][:ii])
		lineEPS.remove()
		lineEPS = axes[3].fill_between(xyzTime[:ii+1], epsArr[:ii+1], color='k', alpha=0.1, label=alabel[3], zorder=1)
		figg.canvas.draw()
		ffID = "%%0%dd" % nd
		figg.savefig((pSeq%ffID)%ii, dpi=600, transparent=True)
		# Parameters to remove border
		trim = "-flatten -fuzz 1%% -trim +repage"
		# Call imagemagick to concatenate/stack the pictures together
		os.system("convert \( %s -resize %dx%d \) \( %s -resize %dx%d \) -gravity SouthEast -background white +smush +10 %s %s" % (seq%ii, vres, vres, (pSeq%ffID)%ii, pres, pres, trim, (cSeq%ffID)%ii))
	return pSeq, cSeq
#---------------------------------------------------------------------------------------#
# Helper functions
#---------------------------------------------------------------------------------------#
def _mktemp(path, seq="lig_vmd_"):
	# Temporary directory to save the generated pictures
	tmpDir = tempfile.mkdtemp(prefix=seq, dir=path)
	print "Saving pictures/config files in temporary directory = %s" % tmpDir 
	return tmpDir
#---------------------------------------------------------------------------------------#
def _genColor(outVec, propVec, col, tmpDir, isscaled=True):
#------------------------------------------------------------------------#
# Generate a coloring scheme from a variable in comFile
#------------------------------------------------------------------------#
	colorName = os.path.join(tmpDir, "col.out")
	xyzTime = outVec.timeList
	propTime = propVec.time()
	# fMatch is the all of the matching/commom timesteps between XYZ and COM files
	# Ind arrays are the indices where these matches exist in the respective arrays
	fMatch = np.intersect1d(propTime, xyzTime)
	# Using ravel is fine because time values are monotonically increasing i.e. one match per "where
	comInd = np.array( [np.where(x==propTime) for x in fMatch] ).ravel()
	xyzInd = np.array( [np.where(x==xyzTime) for x in fMatch] ).ravel()

	# Make sure that the data for all XYZ snapshots is present in the COM file
	print "MATCH:\t", len(xyzTime), len(xyzTime[xyzInd])
	assert len(xyzTime) == len(xyzTime[xyzInd])
	propVec.headMask = propVec.headLoc[comInd]

	prop = propVec._genSingle(col)
	if isscaled:
		# For frame-wise scaling
		#pMin = prop.min(axis=1)[:,np.newaxis]
		#pRange = prop.max(axis=1)[:,np.newaxis] - pMin
		# Global scaling
		# Handles divide by zero; returns 0
		prop -= prop.min()
		pRange = float(prop.max() - prop.min())
		if pRange != 0:
			prop /= pRange
		else:
			prop = np.zeros_like(prop)

	with open(colorName, 'w') as outFile:
		outFile.write(ps._toString2D( prop ))

	# Get number of atoms/molecule and number of molecules
	nMol = propVec.nChunk
	nAtom = int( float(outVec.nAtom[0])/nMol )

	return colorName, nMol, nAtom
#---------------------------------------------------------------------------------------#
def _genVMD(tclFile, varList, tmpDir):
#------------------------------------------------------------------------#
# Generate a string of VMD snapshots in $CWD/tmpDir
#------------------------------------------------------------------------#
	if (os.path.realpath(tmpDir) == os.getcwd()):
		raise ValueError("Temporary directory cannot be the CWD to prevent race conditions.")

	varFile = os.path.join(tmpDir, "vars.inp")
	varList.append(tmpDir)
	outStr = "\n".join( map( str, varList ) )
	with open(varFile, 'w') as varTemp:
		varTemp.write( outStr )

	# Absolute path of VMD executable
	#vmdAbs = os.path.join( os.path.expanduser('~'), ".local/bin/vmd")
	vmdAbs = "vmd"
	cmd = "%s -dispdev text -eofexit -args %s < %s > out.log" % (vmdAbs, varFile, tclFile)
	os.system(cmd)
	return tmpDir
#---------------------------------------------------------------------------------------#
def _genMovie(outLoc, seq, form, imDelay=25):
#------------------------------------------------------------------------#
# Generate movie in the CWD
#------------------------------------------------------------------------#
	form = form.strip('.')
	if form == "gif":
		# imDelay: time in 1/100th of a second
		cmd = "convert -delay %d -dispose background -quality 100 -loop 0 %s %s" % ( imDelay, seq, outLoc )

	elif form == "mp4":
		ffmpeg = os.path.join( os.path.expanduser('~'), ".local/bin/ffmpeg-4.0.2-64bit-static/ffmpeg")
		# Convert imDelay to framerate
		framerate = 100.0/imDelay
		cmd = "%s -framerate %f -i %s -c:v libx264 -vf scale=1280:-2 -pix_fmt yuv420p -y %s" % (ffmpeg, framerate, seq, outLoc)

	print cmd
	os.system(cmd)
	return outLoc
#---------------------------------------------------------------------------------------#
def main(args):
#---------------------------------------------------------------------------------------#
	xyzVec = ps.TRAJVec(args.xyzName)
	nXYZ = xyzVec._mask.size
	xyzTime = xyzVec.timeList
	print "Number of frames = %d" % nXYZ
	nf0 = 0
	if (args.color != None):
		propName, col = args.color
		propVec = ps.MolVec(propName)#, start=1200, skip=10)
		propTime = propVec.time()
		# fMatch is the all of the matching/commom timesteps between XYZ and COM files
		# Ind arrays are the indices where these matches exist in the respective arrays
		fMatch = np.intersect1d(propTime, xyzTime)
		# Using ravel is fine because time values are monotonically increasing i.e. one match per "where
		comInd = np.array( [np.where(x==propTime) for x in fMatch] ).ravel()
		xyzInd = np.array( [np.where(x==xyzTime) for x in fMatch] ).ravel()
		nf0 = xyzInd[0]
#---------------------------------------------------------------------------------------#
# If reduce is true
#---------------------------------------------------------------------------------------#
	if (args.frames != None) and (nXYZ > args.frames):
		xyzVec._mask = np.linspace(nf0, nXYZ-1, args.frames).astype(np.int)
		outXYZ = xyzVec.reduce(addStr='_t%d' % args.frames, override=False)
		print outXYZ
		outVec = ps.TRAJVec(outXYZ)
	else:
		if (nXYZ <= args.frames):
			print "Lesser number of frames in the trajectory than requested. Skipping step..."
			args.frames = nXYZ
		outXYZ = args.xyzName
		outVec = xyzVec
#---------------------------------------------------------------------------------------#
# If coloring is required
#---------------------------------------------------------------------------------------#
	tclFile = "baseGen.tcl"
	colorName = None
	nPix = 4096
	varList = [outXYZ, nPix]
#---------------------------------------------------------------------------------------#
# Compatibilty with SLURM scheduler
#---------------------------------------------------------------------------------------#
	if (args.submit):
		try:
			print os.environ["HOME"]
			tmpDir = os.environ["MYTMP"]
		except KeyError:
			print "Please set the environment variable MYTMP"
			sys.exit(1)
		shutil.copy2("run.sh", os.path.join(tmpDir, "run.sh"))
	else:
		cwd = os.getcwd()
		tmpDir = _mktemp(cwd, seq="lig_vmd_")
	# Store the input arguments
	with open(os.path.join(tmpDir, "in.args"), 'w') as argFile:
		argFile.write(str(vars(args)))
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
	if (args.color != None):
		propName, col = args.color
		colorName, nMol, nAtom = _genColor( outVec, propVec, int(col), tmpDir )
		tclFile = "baseGen.tcl"
		varList += [colorName, nAtom, nMol]

	xyzTime = outVec.time()
	tmpDir = _genVMD(tclFile, varList, tmpDir)
	nd = int(np.log10(len(xyzTime))) + 2
	ffID = "%%0%dd" % nd
	#seq = "test-trans-%s.png" % ffID
	seq = "test-%s.tga" % ffID
	seq = os.path.join(tmpDir, seq)
#---------------------------------------------------------------------------------------#
# If simulataneous plot is required
#---------------------------------------------------------------------------------------#
	if (args.plot):
#---------------------------------------------------------------------------------------#
		propName = ps.flipPath(args.xyzName, 'or.', locType="End", newLoc="Front")
		pSeq, cSeq = _plot(propName, xyzTime, seq, tmpDir)
		# Make combined movie
		movName = ps.altPath(outXYZ, args.form, locType="End", addStr="-comb")
		#movName = "comb.%s" % (args.form).strip('.')
		movPath = _genMovie(movName, cSeq%ffID, form=args.form, imDelay=25)
#---------------------------------------------------------------------------------------#
# Generate movie
#---------------------------------------------------------------------------------------#
	movName = ps.altPath(outXYZ, args.form, locType="End")
	#movName = "out.%s" % (args.form).strip('.')
	movPath = _genMovie(movName, seq, form=args.form, imDelay=25)
	#shutil.rmtree(tmpDir)

	# Without the ligands
	#seq = "NL-test-trans-%s.png" % ffID
	seq = "NL-test-%s.tga" % ffID
	seq = os.path.join(tmpDir, seq)
	movName = ps.altPath(outXYZ, args.form, locType="End", addStr="NL")
	movPath = _genMovie(movName, seq, form=args.form, imDelay=25)

	ps._del("patch.pyc")
	ps._del("parser.pyc")
#---------------------------------------------------------------------------------------#
# MAIN
#---------------------------------------------------------------------------------------#
# NOT TO BE RUN BY ITSELF; RUN THROUGH WRAPPER SCRIPT JOB.PY
#---------------------------------------------------------------------------------------#
if __name__ == "__main__":
	if len(sys.argv) == 2: 
		# Get dictionary from string
		iDict = ast.literal_eval(sys.argv[1])
	else:
		raise ValueError("Only accepts one argument in the form of a dictionary")

	# Convert to argparse.Namespace
	args = argparse.Namespace(**iDict)
	main(args)
