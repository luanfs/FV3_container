#!/bin/bash
directory="FV3_container" 
output=$directory".tar.bz2"


#-------------------------------------------------------------------------------------------------------
#First we send the code to dropbox
dropdir="/home/luanfs/Dropbox/doc/code/"

if [ -d "$dropdir$directory" ]; then
    echo "Directory $dropdir$directory already exists."
else
    cd $dropdir
    mkdir "$directory"
    echo "Directory $dropdir$directory created successfully."
    cd -
    pwd
fi

echo "Sync with Dropbox:"
rsync -v -t -u $output  "$dropdir$directory/."
echo "Synchronized with Dropbox"
echo

#-------------------------------------------------------------------------------------------------------
# remote host 1 - ime.usp.br
user_remote_host1="luansantos"
remote_host1="ime.usp.br"
remote_host1_dir="/var/tmp/lfs"
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# remote host 2 - ybytu
user_remote_host2="luansantos"
remote_host2="ybytu.ime.usp.br"
remote_host2_dir="/home/luansantos/$directory"
#-------------------------------------------------------------------------------------------------------

 
#-------------------------------------------------------------------------------------------------------
#remote server ime.usp.br backup sync
echo "Sending to $remote_host1:"
rsync -t -v -z -a --progress $output $user_remote_host1@$remote_host1:$remote_host1_dir
echo "Sent to $user_remote_host1@$remote_host1"
echo
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#remote server ybytu backup sync
echo "Sending to $remote_host2:"
ssh -t $user_remote_host1@$remote_host1 "rsync -t -v -z -a --progress $remote_host1_dir/$output $user_remote_host2@$remote_host2:$remote_host2_dir; rm -rf $output"
echo "Sent to $user_remote_host2@$remote_host2"
echo
#-------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------
# untar and compile
echo "Untar at $remote_host2"
ssh -t $user_remote_host1@$remote_host1 "ssh -t $user_remote_host2@$remote_host2 <<EOF
	cd $remote_host2_dir;
	tar -xvf $output;
	rm -rf $output;
EOF"
echo "Untar at $remote_host2 done."
#-------------------------------------------------------------------------------------------------------

# remove tar file
rm -rf $output
