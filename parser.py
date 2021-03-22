import numpy as np
import os
import sys
import patch as pt
#---------------------------------------------------------------------------------------#
# Helper functions
#---------------------------------------------------------------------------------------#
# Wrapper function over sed/bash utility
def wsed(cmd,tmp='tmp'):
	print "Generating header list...\n%s" % cmd
	os.system( "%s > %s" % (cmd, tmp) )
	with open(tmp, 'r') as inFile:
		data = np.array( inFile.readlines() )
	_del(tmp)
	return data
#---------------------------------------------------------------------------------------#
def _del(fPath):
	if os.path.exists(fPath):
		os.system('rm -rf %s' % fPath)
#---------------------------------------------------------------------------------------#
def _toString(fList):
	return "\t".join( map(str, fList) )
#---------------------------------------------------------------------------------------#
def _toString2D(fList):
	return "\n".join( [ "\t".join( map(str, x) ) for x in fList ] )
#---------------------------------------------------------------------------------------#
def altPath(path, newExt, locType="Front", addStr=None):
	locDict = { "Front" : 0, "End" : -1}
	dirList = path.split(os.sep)
	filename = dirList[-1].split('.')
	loc = locDict[locType]
	filename[loc] = newExt.strip('.')
	if addStr != None:
		filename[loc-1] += addStr
	dirList[-1] = '.'.join(filename)
	return os.sep.join(dirList)
#---------------------------------------------------------------------------------------#
def flipPath(path, newExt, locType="Front", newLoc="Front", addStr=None):
	locDict = { "Front" : 0, "End" : -1}
	dirPath, fileName = os.path.split(path)
	base = fileName.split('.')
	base.pop(locDict[locType])

	loc = (locDict[newLoc] + len(base)) % len(base)
	base.insert(loc, newExt.strip('.'))
	return os.path.join(dirPath, '.'.join(base))
#---------------------------------------------------------------------------------------#
def dict_prop(prop, header, data):
	data = np.array(data)
	return data[:, header.index(prop)]
#---------------------------------------------------------------------------------------#
def quickRead(filename):
	with open( filename ) as inFile:
		data = np.array( [ [float(x) for x in y.split()] for y in inFile.readlines() ] )
	return data
#---------------------------------------------------------------------------------------#
# Calculate the total number of times "run" command has been invoked in the simulation
#---------------------------------------------------------------------------------------#
def runCount(filename, key):
	cval = 0
	with open(filename) as inFile:
		for line in inFile:
			if key in line.split():
				cval += 1
	return cval
#---------------------------------------------------------------------------------------#
# Read standard LAMMPS output file (log.lammps/thermo output)
#---------------------------------------------------------------------------------------#
def read(infile, start, stop='Loop'):
	lines = []
	count = 0
	for line in infile:
		try:
			header = line.split()
			if start in header:
				break
		except IndexError:
			continue

	for line in infile:
		try:
			if line.split()[0] == stop:
				break
		except IndexError:
			break
		try:
			lines.append([float(x) for x in line.split()])
		except ValueError:
			count += 1
			continue
	if count > 0:
		return lines, header, count, 1
	else:
		return lines, header, count, 0

#---------------------------------------------------------------------------------------#
# Get box dimensions form out.* file
#---------------------------------------------------------------------------------------#
def get_box(datafile, time_series):
	filename = altPath(datafile, 'out')

	# Count the total number of runs in the file
	cval = runCount(filename, "TotEng")
	# Get box dimensions
	dList = []
	runData = []
	gWarn = 0
	gCount = 0
	RC = 0
	with open(filename) as infile:
		for c in xrange(cval):
			RC += 1
			data, header, count, warn = read(infile, "TotEng")
			gWarn = (warn or gWarn)
			gCount += count
			runData.append(len(data)-1)
			if RC == 0:
				dList = dList + data
			else:
				dList = dList + data[1:]

	if gWarn == 1:
		print "WARNING: %d warnings found in file, " % gCount, filename

	df = np.array(dList)
	tList = dict_prop('Time', header, df)
	mask = np.in1d(tList, time_series)
	modList = df[mask]
	#df = pd.DataFrame(data=dList, index=None, columns=header)
	#modList = df[df['Time'].isin(time_series)]
	time = dict_prop("Time", header, modList)
	xlo = dict_prop("Xlo", header, modList)
	xhi = dict_prop("Xhi", header, modList)
	ylo = dict_prop("Ylo", header, modList)
	yhi = dict_prop("Yhi", header, modList)
	zlo = dict_prop("Zlo", header, modList)
	zhi = dict_prop("Zhi", header, modList)
	return time, np.c_[ (xhi-xlo), (yhi-ylo), (zhi-zlo) ], runData
#---------------------------------------------------------------------------------------#
# Get box dimensions form out.* file
#---------------------------------------------------------------------------------------#
def get_dims(datafile, time_series):
	import pandas as pd
	filename = altPath(datafile, 'out')

	# Count the total number of runs in the file
	cval = runCount(filename, "TotEng")
	print cval
	# Get box dimensions
	dList = []
	gWarn = 0
	gCount = 0
	RC = 0
	with open(filename) as infile:
		for c in xrange(cval):
			RC += 1
			data, header, count, warn = read(infile, "TotEng")
			gWarn = (warn or gWarn)
			gCount += count

			if RC == 0:
				dList = dList + data
			else:
				dList = dList + data[1:]

	if gWarn == 1:
		print "WARNING: %d warnings found in file, " % gCount, filename

	df = pd.DataFrame(data=dList, index=None, columns=header)
	modList = df[df['Time'].isin(time_series)]

	time = dict_prop("Time", header, modList)
	lx = dict_prop("Lx", header, modList)
	ly = dict_prop("Ly", header, modList)
	lz = dict_prop("Lz", header, modList)

	return time, np.c_[lx, ly, lz]
#---------------------------------------------------------------------------------------#
def getOP(fFile, Vec, rCols, boxDims):
	# Pre-compile f77 code
	cmd = "ifort -o a.out %s -O2" % fFile
	print cmd
	os.system(cmd)

	lines = []
	for ti, tStep in enumerate(Vec.timeList()):
		if len(boxDims) <= ti:
			bVal = boxDims[-1]
		else:
			bVal = boxDims[ti]
		with open ("coord.inp", 'w') as tempFile:
			tempFile.write( '%s\n' % _toString( [bVal[0], 0.0, 0.0] ) )
 			tempFile.write( '%s\n' % _toString( [0.0, bVal[1], 0.0] ) )
			tempFile.write( '%s\n' % _toString( [0.0, 0.0, bVal[2]] ) )

			tempFile.write( '%s\n' % '\n'.join( [ _toString(x) for x in ( Vec._genStep(ti) ) [:,rCols] ] ) )

		os.system("./a.out")

		with open("NTr.out") as ntrFile:
			lines.append([ float(x) for x in ntrFile.readline().split() ])

	return np.array(lines)
#---------------------------------------------------------------------------------------#
# Molecular vector output file:
# Each block of data (representing a timestep) is demarcated by a header line
#---------------------------------------------------------------------------------------#
class MolVec:
	def __init__(self, filename, start=None, skip=None, end=None):
		self.filename = filename
		headLoc, rawData = self._genLoc()
		self.headLoc = headLoc
		self.rawData = rawData
		self._setmask(start, skip, end)

#---------------------------------------------------------------------------------------#
	# mask = nskip = (# of timesteps to skip + 1)
	def _setmask(self, start=None, skip=None, end=None):
		conv = lambda x, y: y if (x == None) else x
		self._start = conv(start,0)
		self._skip = conv(skip,1)
		self._end = conv(end,len(self.headLoc))
		assert (type(self._start) == int) and (type(self._skip) == int) and (type(self._end) == int)
		self.headMask = self.headLoc[self._start:self._end:self._skip]
		# re-/initialize all variables here

#---------------------------------------------------------------------------------------#
	def time(self):
		return np.array( list(self.rawData[self.headMask]), dtype=np.float ) [:,0]
#---------------------------------------------------------------------------------------#
	# Generate the header and rawData of the file
	def _genLoc(self):
		with open(self.filename) as inFile:
			# Ignore first 3 title lines
			data = np.array( [ np.array( x.split() ).astype(float) for x in inFile.readlines() [3:] ] )
		# Collect location of header rows from the array
		# "header" rows: Timestep, nrows
		headLoc = []
		lNo = 0
		while True:
			if lNo <= len(data)-1:
				headLoc.append( lNo )
				lNo += int(data[lNo][1]) + 1
			else:
				break
		headLoc = np.array(headLoc)
		# Confirm that number of chunks (molecules) does not change during the simulation
		nChunk = list( set(headLoc[1:]-headLoc[:-1]) )
		assert len( nChunk ) == 1
		self.nChunk = nChunk[0]-1
		return headLoc, data
#---------------------------------------------------------------------------------------#
	# Generate a single snapshot (time) as a numpy array slice
	def _genStep(self, step):
		loc = self.headMask[step]
		tStep, nChunks = self.rawData[loc]
		return np.array( list( self.rawData[loc+1:loc+1+int(nChunks)] ), dtype=np.float )
#---------------------------------------------------------------------------------------#
	# Generate a timeseries data of one property for one molecule as a numpy array slice
	def _genMol(self, molID, col):
		assert molID <= self.nChunk
		return np.array( list(self.rawData[self.headMask+molID]), dtype=np.float ) [:,col]
#---------------------------------------------------------------------------------------#
  # Generate a timeseries for all molecules for one property
	def _genProp(self, col):
		# Index array "ind" is a 2D array that stores indices of all molecules for a 
		# given timestep (in headMask) in a row, and each row obviously represents a timestep
		ind = np.array( [ np.arange(x+1, x+1+self.nChunk, 1) for x in self.headMask ] )
		# Extract only one column data
		propList = np.array( map(list,self.rawData[ind.ravel()]), dtype=np.float ) [:,col]

		if type(col) == int:
			col = [col]

		propList.reshape([-1,self.nChunk,len(col)])

		return propList.reshape([-1, self.nChunk, len(col)])
#---------------------------------------------------------------------------------------#
	# Generate a timeseries for all molecules for one property
	def _genSingle(self, col):
		# Index array "ind" is a 2D array that stores indices of all molecules for a 
		# given timestep (in headMask) in a row, and each row obviously represents a timestep
		ind = np.array( [ range(x+1, x+1+self.nChunk, 1) for x in self.headMask ] )
		# Extract only one column data
		propList = np.array( list(self.rawData[ind.ravel()]), dtype=np.float ) [:,col]

		return propList.reshape([-1, self.nChunk])
#---------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------#
# Molecular vector output file:
# Each block of data (representing a timestep) is demarcated by a header line
#---------------------------------------------------------------------------------------#
class XYZVec:
	def __init__(self, filename, start=None, skip=None, end=None, nCols=4):
		self.filename = filename
		self.headLoc, nAtom, self.timeList = self._genLoc()
		self.nAtom = list ( set (nAtom) )
		# Closed (mass-wise) system: no change in the number of atoms as time passes
		assert len(self.nAtom) == 1
		self.nAtom = self.nAtom[0]
		# Pre-allocate zeros array to store 
		self.rawData = np.zeros( [ self.headLoc.size, self.nAtom, nCols ] )
		# Generate a mask (index values) of time points to be considered from global headLoc
		self._mask = self._genMask(start, skip, end)
		self.headMask = self.headLoc[self._mask]

#---------------------------------------------------------------------------------------#
	# Convert input start skip, end to correct integer value
	def _conv(self, val, alt):
		if ( val == None ):
			return alt
		elif ( (type(val) == int) and (val < len(self.headLoc)) ):
			return val
		elif ( ( type(val) == float) and (0.0 <= val <= 1.0) ):
			return val*len(self.headLoc)
		else:
			raise ValueError("Incorrect start/stop/skip value for the XYZ vector.")
#---------------------------------------------------------------------------------------#
	# mask = nskip = (# of timesteps to skip + 1)
	def _genMask(self, start, skip, end):
		lenh = self.headLoc.size
		self._start = self._conv(start,0)
		self._skip = self._conv(skip,1)
		self._end = self._conv(end,lenh)
		assert ( (self._start < self._end) and (self._skip <= len(self.headLoc)) )
		mask = np.arange(self._start, self._end, self._skip)
		# re-/initialize all variables here
		return mask
#---------------------------------------------------------------------------------------#
	# Generate the header and rawData of the file
	def _genLoc(self, strVal="Atoms.", start_l=1):
		outTemp = altPath(self.filename, '.head', locType="End")
		# Print separated by newline in groups of 3:
		# line number\nNumber of beads\nAtoms. Timestep: t\n
		print "Generating header list..."
		cmd = 'sed -n \'%d,$ {/%s/{=;x;p;x;p};h}\' < %s > %s' % (start_l, strVal, self.filename, outTemp)
		print cmd
		os.system(cmd)	
		with open(outTemp, 'r') as inFile:
			data = inFile.readlines()
		_del(outTemp)

		headLoc = np.array( data[::3], dtype=np.int )
		nAtom = np.array( data[1::3], dtype=np.int )
		timeList = np.array( [ x.split()[2] for x in data[2::3] ], dtype=np.int )
		return headLoc, nAtom, timeList
#---------------------------------------------------------------------------------------#
	# Generate a single snapshot (time) as a numpy array slice
	def _genStep(self, step):
		loc = self.headMask[step]
		tStep, nChunks = self.rawData[loc]
		return np.array( list( self.rawData[loc+1:loc+1+int(nChunks)] ), dtype=np.float )
#---------------------------------------------------------------------------------------#
	# Save all masked timepoints from the XYZ file as a numpy array
	def _fillData(self):
		for ti, tVal in enumerate(self._mask):
			self._extract(tVal)
#---------------------------------------------------------------------------------------#
	# Save all masked timepoints from the XYZ file as a numpy array
	def _genCom(self, molList):
		# MolList has the length of nMol and each element is the number of atoms/molecule
		cumList = np.insert(np.cumsum(molList),0,0)
		for ti, tVal in enumerate(self._mask):
			self._extract(tVal)
			comList = np.array([ np.mean(self.rawData[1][cumList[i]:cumList[i+1],1:], axis=0) for i in range(len(cumList)-1) ])
#---------------------------------------------------------------------------------------#
	# Generate a single snapshot (time) from the XYZ file as a numpy array
	def _extract(self, tVal):
		# tVal is the relative position of timestep in [0,len(self.headLoc)]
		self._conv(tVal,-1)
		outTemp = 'tmp'
		# If the data for the required timestep has not been filled
		if not self.rawData[tVal].any():
			# loc is the line number of the Atoms.... line for the required timestep
			loc = self.headLoc[tVal]
			start_l = loc + 1
			end_l = loc + self.nAtom
			# Print lines between start_l and end_l (both inclusive):
			cmd = 'sed -n \'%d,%dp;%dq\' < %s > %s' % (start_l, end_l, end_l, self.filename, outTemp)
			print cmd
			os.system(cmd)	
			with open(outTemp, 'r') as inFile:
				self.rawData[tVal][:,:4] = np.array( [ x.split() for x in inFile.readlines() ], dtype=np.float)
			_del(outTemp)
#---------------------------------------------------------------------------------------#
	# Save all masked timepoints from the XYZ file into a reduced file
	def reduce(self, addStr=None, override=False):
		if (addStr == None):
			addStr = '_t%d-%d-%d' % (self._start, self._end, self._skip)
		outName = altPath(self.filename,'.xyz', locType="End", addStr=addStr)
		# If file already exists and no override flag
		if os.path.exists(outName) and (not override):
			return outName
		#	Otherwise re-run
		_del(outName)
		sedStr = ''
		for ti, tVal in enumerate(self._mask):
			# loc is the line number of the Atoms.... line for the required timestep
			loc = self.headLoc[tVal]
			start_l = loc - 1
			end_l = loc + self.nAtom
			# Print lines between start_l and end_l (both inclusive):
			sedStr += '%d,%dp;' % (start_l, end_l)

		cmd = 'sed -n \'%s%dq\' < %s >> %s' % (sedStr, end_l, self.filename, outName)
		print cmd
		os.system(cmd)	
		return outName
#---------------------------------------------------------------------------------------#
# LAMMPS Trajectory: Molecular vector
# Each block of data (representing a timestep) is demarcated by a header line
#---------------------------------------------------------------------------------------#
class TRAJVec:
	def __init__(self, filename, start=None, skip=None, end=None):
		self.filename = filename
		self._genLoc()
		# Generate a mask (index values) of time points to be considered from global headLoc
		self._mask = self._genMask(start, skip, end)
		self.headMask = self.headLoc[self._mask]
#---------------------------------------------------------------------------------------#
	# Convert input start skip, end to correct integer value
	def _conv(self, val, alt):
		if ( val == None ):
			return alt
		elif ( (type(val) == int) and (val < len(self.headLoc)) ):
			return val
		elif ( ( type(val) == float) and (0.0 <= val <= 1.0) ):
			return val*len(self.headLoc)
		else:
			raise ValueError("Incorrect start/stop/skip value for the XYZ vector.")
#---------------------------------------------------------------------------------------#
	# mask = nskip = (# of timesteps to skip + 1)
	def _genMask(self, start, skip, end):
		lenh = self.headLoc.size
		self._start = self._conv(start,0)
		self._skip = self._conv(skip,1)
		self._end = self._conv(end,lenh)
		assert ( (self._start < self._end) and (self._skip <= len(self.headLoc)) )
		mask = np.arange(self._start, self._end, self._skip)
		# re-/initialize all variables here
		return mask
#---------------------------------------------------------------------------------------#
	def time(self):
		return self.timeList[self._mask]
#---------------------------------------------------------------------------------------#
	# Generate the header and rawData of the file
	def _genLoc(self, start_l=1):
		self.tmp = altPath(self.filename, '.head', locType="End")
		# Print timestep in next line as matched line
		cmd = 'sed -n \'%d,$ {/%s/{n;p};h}\' < %s' % (start_l, "TIMESTEP", self.filename)
		self.timeList = wsed(cmd,self.tmp).astype(np.int)
		# Print number of beads in next line as matched line
		cmd = 'sed -n \'%d,$ {/%s/{n;p};h}\' < %s' % (start_l, "NUMBER", self.filename)
		self.nAtom = wsed(cmd,self.tmp).astype(np.int)
		# Print box bounds in next line as matched line
		cmd = 'sed -n \'%d,$ {/%s/{n;p;n;p;n;p};h}\' < %s' % (start_l, "BOX", self.filename)
		box = np.array( [ np.array( x.split()) for x in wsed(cmd,self.tmp)], dtype=np.float )
		self.box = box.reshape(-1,6)
		# Print line number of per-timestep data
		cmd = 'sed -n \'%d,$ {/%s/{=};h}\' < %s' % (start_l, "ITEM: ATOM", self.filename)
		self.headLoc = wsed(cmd,self.tmp).astype(np.int)
		# Print line number of timestep as matched line
		cmd = 'sed -n \'%d,$ {/%s/{=};h}\' < %s' % (start_l, "TIMESTEP", self.filename)
		f_start = wsed(cmd,self.tmp).astype(np.int)
		# Total number of lines in file
		cmd = 'wc -l < %s' % (self.filename)
		lend = wsed(cmd,self.tmp).astype(np.int)[0]
		f_end = np.r_[f_start[1:]-1, lend]
		# Delimiters for the start and end each frame
		self.delim = np.c_[f_start,f_end]
#---------------------------------------------------------------------------------------#
	# Save all masked timepoints from the XYZ file as a numpy array
	def _fillData(self):
		arr = []
		for ti, tVal in enumerate(self._mask):
			arr.append( self._extract(tVal) )
		return arr
#---------------------------------------------------------------------------------------#
	# Generate a single snapshot (time) from the XYZ file as a numpy array
	def _extract(self, tVal=None):
		# tVal is the relative position of timestep in [0,len(self.headLoc)]
		self._conv(tVal,-1)
		# Print lines between start_l and end_l (both inclusive):
		start_l, end_l = self.delim[tVal]
		cmd = 'sed -n \'%d,%dp;%dq\' < %s' % (start_l, end_l, end_l, self.filename)
		frame = np.array( [ np.array( x.split()) for x in wsed(cmd,self.tmp)], dtype=np.float )
		return frame.reshape(self.nAtom[tVal], -1)
#---------------------------------------------------------------------------------------#
	# Save all masked timepoints from the XYZ file into a reduced file
	def reduce(self, addStr=None, override=False):
		if (addStr == None):
			addStr = '_t%d-%d-%d' % (self._start, self._end, self._skip)
		outName = altPath(self.filename,'.lammpstrj', locType="End", addStr=addStr)
		# If file already exists and no override flag
		if os.path.exists(outName) and (not override):
			return outName
		#	Otherwise re-run
		_del(outName)
		sedStr = "".join( ["%d,%dp;" % tuple(x) for x in self.delim[self._mask] ] )
		end_l = self.delim[self._mask[-1], 1]
		# Delimiters for the start and end each frame
		cmd = 'sed -n \'%s%dq\' < %s >> %s' % (sedStr, end_l, self.filename, outName)
		print cmd
		os.system(cmd)
		return outName
#---------------------------------------------------------------------------------------#
if (__name__ == "__main__"):
#---------------------------------------------------------------------------------------#
	if len(sys.argv) == 2:
		fileName = sys.argv[1]
	else:
		raise ValueError("Incorrect number of input arguments.")
#---------------------------------------------------------------------------------------#
	rVec = MolVec(fileName)

