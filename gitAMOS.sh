#/bin/sh
git clone git://amos.git.sourceforge.net/gitroot/amos/amos 
./bootstrap
./configure --prefix=/usr
make install