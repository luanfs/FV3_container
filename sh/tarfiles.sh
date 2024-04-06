#!/bin/bash
# Creates a tarball of the bacukp files
directory="FV3_container"

# Source code files
shield_src_dir="../SHiELD_SRC"
shield_build_dir="../SHiELD_build"

srcdir="../SHiELD_SRC/GFDL_atmos_cubed_sphere"
sourcefiles="$srcdir/tools $srcdir/driver $srcdir/model $srcdir/GFDL_tools"
driverfiles="../SHiELD_SRC/atmos_drivers/solo/*.F90"

# Plot files
plotfiles="../plot/*.py"

# Bash files
shfiles="../sh/*sh"

# Test files
testdir="../SHiELD_SRC/test"
testfiles="$testdir/*sh"

#files="$shield_src_dir $shield_build_dir $plotfiles $shfiles"
#files="$sourcefiles $plotfiles $shfiles $testfiles $driverfiles"
files="$plotfiles $shfiles $testfiles $driverfiles $sourcefiles"

# Output name
output=$directory".tar.bz2"

# Creates the tarball
tar cjfv $output $files

echo "File " $output " ready!"
echo
