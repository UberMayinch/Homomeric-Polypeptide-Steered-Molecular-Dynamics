lappend auto_path /home/chinmay/Desktop/vmd/scripts/la1.0
lappend auto_path /home/chinmay/Desktop/vmd/scripts/orient

package require Orient
namespace import Orient::orient

# rotating such that the long axis aligns with z-axis
set sel [atomselect top "protein"]
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 1] {0 1 0}]
$sel move $A
set I [draw principalaxes $sel]

# translating such that N of the first residue is at the origin
set atom [atomselect top "index 0"]
set pos [lindex [$atom get {x y z}] 0]
set all [atomselect top "all"]
$all moveby [vecscale -1.0 $pos]

set molid [molinfo top]
set filename [molinfo $molid get filename]

# Extract just the base name (remove path)
set basename [file tail $filename]

# Remove extension (optional, if you want custom suffix)
set rootname [file rootname $basename]

# Create new filename
set outname "${rootname}_cleaned.pdb"

$all writepdb $outname
