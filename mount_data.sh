bw login
bw unlock
export USER=$(bw get username 8e4a8664-b879-4395-9850-ad19007e871a)
export PASSWD=$(bw get password 8e4a8664-b879-4395-9850-ad19007e871a)

sudo mount -t cifs //uw.lu.se/research /home/yakir/mnt/ -o user=$USER,password=$PASSWD,uid=$(id -u),gid=$(id -g)
