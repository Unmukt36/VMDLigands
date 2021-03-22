#-----------------------------------------------#
proc lexists {name} {
    expr {![catch {file lstat $name finfo}]}
}
#-----------------------------------------------#
proc repRCO {res}		{
	# Deletes top representation
	#mol delrep 0 top
	mol representation CPK 2.10000 0 res 0
	mol color Element
	mol selection {element H}
	mol addrep top

	mol representation CPK 1.45000 0 res 0
	mol color Element
	mol selection {element He}
	mol addrep top

	mol representation CPK 1.15000 0 res 0
	mol color Element
	mol selection {element Li}
	mol addrep top

	mol representation CPK 1.00000 0 res 0
	mol color Element
	mol selection {element Be}
	mol addrep top
}
#-----------------------------------------------#
proc removeAxes {}	{
	axes Location off
}
#-----------------------------------------------#
proc delall {molid} {
	# delete all representations
	set numrep [molinfo $molid get numreps]
	for {set i 0} {$i < $numrep} {incr i} {
		mol delrep 0 $molid
	}
}
#-----------------------------------------------#
proc comLim {aSelect nMol} {
	upvar $aSelect aLocal
	set xcom {}
	set ycom {}
	set zcom {}
	set lim {}
	# Get the limits of the center of masses
	for {set i 0} {$i < $nMol} {incr i} {
		set xi [measure center $aLocal($i)]
		lappend xcom [lindex $xi 0]
		lappend ycom [lindex $xi 1]
		lappend zcom [lindex $xi 2]
	}
	lappend lim [lindex [lsort -real $xcom] 0] [lindex [lsort -real $xcom] end]
	lappend lim [lindex [lsort -real $ycom] 0] [lindex [lsort -real $ycom] end]
	lappend lim [lindex [lsort -real $zcom] 0] [lindex [lsort -real $zcom] end]
	set rX [expr [lindex $lim 1] - [lindex $lim 0]]
	set rY [expr [lindex $lim 3] - [lindex $lim 2]]
	set rZ [expr [lindex $lim 5] - [lindex $lim 4]]
	set rmax [expr max($rX,$rY)]
	return rmax
}
#-----------------------------------------------#
proc colorby {args}		{
	global aSelect fUser molid nMol
	# Gloablly defined variable names
	set f [molinfo $molid get frame]
	for {set i 0}	{$i < $nMol}	{incr i}	{
		$aSelect($i) set user [lindex $fUser($f) $i]
	}
}
#-----------------------------------------------#
proc readFile {fname fUser nLines}		{
	if {[lexists $fname] == 1}	{
		upvar $fUser fLocal
		puts [format "Reading %s" $fname]
		set fp [open $fname r]
		# n is a counter for number of lines in a file (wc -l fails)
		set n 0
		foreach line [split [read $fp] \n]	{
			set fLocal($n) $line
			incr n
		}
		close $fp
		# Ensure that number of lines match the given value of nLines
		if {$nLines != $n}	{
			error "ERROR: Incorrect # of lines to be read from file."
			exit
		}
	}	else	{
		error [format "File %s does not exist. Exiting..." $fname]
		exit
	}
}
#-----------------------------------------------#
proc selectMol {molid aSelect nAtom nMol}		{
	upvar $aSelect aLocal
	# Atom selection objects array
	for {set i 0}	{$i < $nMol}	{incr i}	{
		set start [expr $i*$nAtom]
		set end [expr ($i+1)*$nAtom-1]
		set aLocal($i) [atomselect $molid [format "index %d to %d" $start $end]]
		$aLocal($i) uplevel 1
	}
}
#-----------------------------------------------#

#-----------------------------------------------#
#	MAIN
#-----------------------------------------------#
package require pbctools
set reso 5000.0000
set zoom 2.5
set rcore 8.5
set rlig 6.5
#source tools.tcl

# Parse command line arguments
# Get command line arguments from an input file
set varFile [lindex $argv 0]

readFile $varFile arg_v 6
set xyzFile $arg_v(0)
set nPix $arg_v(1)
set colFile $arg_v(2)
set nAtom $arg_v(3)
set nMol $arg_v(4)
set tmpDir $arg_v(5)

# Load the XYZ file
mol new $xyzFile waitfor all
set molid [molinfo top]
delall $molid
#-----------------------------------------------#
display projection orthographic
# Sets white background and depth cue
color Display {Background} white
display depthcue on
color Name 1 gray
color Name 2 gray
color Name 3 gray
color change rgb red2 1.0 0.35 0.35
color Name 4 red2
color change rgb green2 0.35 1.0 0.35
color Name 5 green2
#-----------------------------------------------#
# Core beads
for {set i 0} {$i < 4} {incr i} {
	mol representation CPK $rcore 0 $reso 0
	mol color Name
	mol selection "name $i and sqr(x)+sqr(y)+sqr(z) > 0.0001"
	mol material AOShiny
	mol addrep $molid
	mol selupdate $i $molid 1
	}
# Ligand beads
for {set i 4} {$i < 6} {incr i} {
	mol representation CPK $rlig 0 $reso 0
	mol color Name
	mol selection "name $i and sqr(x)+sqr(y)+sqr(z) > 0.0001"
	mol material AOShiny
	mol addrep $molid
	mol selupdate $i $molid 1
	}
# Clip the y-axis planes
#set dim [lindex [pbc get] 0]
#lset dim 1 200.0
#set dim [lmap x [pbc get -all] {lset x 1 200.0}]
set dim [pbc get -all]
set rlim 0.0
for { set i 0 } { $i < [llength $dim] } { incr i } {
	lset dim $i 1 200.0
	set rlim [expr max( [::tcl::mathfunc::max {*}[lrange [lindex $dim 1] 0 2]], $rlim )]
}
# Set VMD variables to clipped $dim
pbc set $dim -all
# Draw the pbc
pbc box -center origin -color black -style dashed -material Transparent -resolution $reso -width 0.5
# Wrap across the PBCs
pbc wrap -all -center origin
# Zoom to an appropriate amout
#set rlim [lindex $dim 0]
# Find largest box size in 
scale to [expr $zoom/$rlim]

# Create new rep based on coloring wrt a user-defined field
# fUser is a per-molecule array that stores eg. efrac, active etc
selectMol $molid aSelect $nAtom $nMol
# File format: NMOL columns, numFrames rows
set n [molinfo $molid get numframes]
readFile $colFile fUser $n

# Scales from Red-White-Blue
mol scaleminmax 0 $molid 0.0 1.0
color scale method BWR
color scale min 0.0
color scale max 1.0
color scale midpoint 0.5

axes Location off
rotate y to 180
rotate x by 90
set xScale [comLim aSelect $nMol]

animate goto start
colorby
# $molid aSelect fUser $nMol

# Total number of frames in the trajectory
set nFrames [molinfo $molid get numframes]
# Pad with appropriate number of zeros, so that order is maintained in convert
set ns [expr int(log10($nFrames))+2]

for {set cFrame 0} {$cFrame < $nFrames}	{incr cFrame}	{
	puts $cFrame
	animate goto $cFrame
	display update
	# Set output name
	set datName [format "%s/scene.dat" $tmpDir]

	# Render with ligands
	mol showrep $molid 4 1
	mol showrep $molid 5 1
	set outName [format "%s/test-%0${ns}d.tga" $tmpDir $cFrame]
	set transName [format "%s/test-trans-%0${ns}d.png" $tmpDir $cFrame]
	set cmd [format "/opt/ohpc/pub/software/vmd-1.9.3/lib/tachyon_LINUXAMD64 -aasamples 12 %s -format TGA -res %d %d -o %s" $datName $nPix $nPix $outName]
	render Tachyon $datName $cmd
	set cmd [format "convert %s -flatten -fuzz 1%% -trim +repage -transparent white %s" $outName $transName]
	puts $cmd
	exec {*}$cmd

	# Render without ligands/ turn-off respective reps
	mol showrep $molid 4 0
	mol showrep $molid 5 0
	set outName [format "%s/NL-test-%0${ns}d.tga" $tmpDir $cFrame]
	set transName [format "%s/NL-test-trans-%0${ns}d.png" $tmpDir $cFrame]
	set cmd [format "/opt/ohpc/pub/software/vmd-1.9.3/lib/tachyon_LINUXAMD64 -aasamples 12 %s -format TGA -res %d %d -o %s" $datName $nPix $nPix $outName]
	render Tachyon $datName $cmd
	set cmd [format "convert %s -flatten -fuzz 1%% -trim +repage -transparent white %s" $outName $transName]
	puts $cmd
	exec {*}$cmd
}
