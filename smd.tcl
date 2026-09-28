# $Id: smd.tcl,v 1.2 2005/02/18 18:07:11 mbach Exp $

# Atoms selected for force application 

global runid

set id1 [atomid PROA 1 N]
set grp1 {}
lappend grp1 $id1
set a1 [addgroup $grp1]

set id2 [atomid PROA 10 N]
set grp2 {}
lappend grp2 $id2
set a2 [addgroup $grp2]

# set the output frequency, initialize the time counter
set Tclfreq 50
set t 0

# contraint points

set c1x -0.04
set c1y  0.31
set c1z -0.40
      
set c2x -6.05
set c2y 1.25
set c2z 14.58

# force constant (kcal/mol/A^2)
set k 7.2

# pulling velocity (A/timestep)
set v 0.00002

#set outfilename da_smd_tcl_1.out
set outfilename "da_smd_tcl_$runid.out"
open $outfilename w

proc calcforces {} {

  global Tclfreq t k v a1 a2 c1x c1y c1z c2x c2y c2z outfilename

  # get coordinates

  loadcoords coordinate

  set r1 $coordinate($a1)
  set r1x [lindex $r1 0]
  set r1y [lindex $r1 1]
  set r1z [lindex $r1 2]

  set r2 $coordinate($a2)
  set r2x [lindex $r2 0]
  set r2y [lindex $r2 1]
  set r2z [lindex $r2 2]

  # calculate forces

  set f1x [expr $k*($c1x-$r1x)]
  set f1y [expr $k*($c1y-$r1y)]
  set f1z [expr $k*($c1z-$v*$t-$r1z)]
  set f1 [list $f1x $f1y $f1z]

  set f2x [expr $k*($c2x-$r2x)]
  set f2y [expr $k*($c2y-$r2y)]
  set f2z [expr $k*($c2z-$r2z)]
  set f2 [list $f2x $f2y $f2z]

  # apply forces

  addforce $a1 $f1
  addforce $a2 $f2

  # output

  set foo [expr $t % $Tclfreq]
  if { $foo == 0 } {
      set outfile [open $outfilename a]
      set time [expr $t*2/1000.0]
      puts $outfile "$time $r1z $f1z"
      close $outfile
  }
  incr t
  return
}
