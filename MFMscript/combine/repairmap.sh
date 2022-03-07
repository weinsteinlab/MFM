#! /bin/bash

dxfiles=`ls proteinfilename/cmb*.dx`
IFS=" " read var1 var2 var3 var4 var5 <<< "$(echo $dxfiles)"

#echo $var1
sed -i -e "1,11 s/\"double\"/double/" -e "1,3 s/#/#\n#/" $var1
sed -i -e "1,11 s/\"double\"/double/" -e "1,3 s/#/#\n#/" $var2
sed -i -e "1,11 s/\"double\"/double/" -e "1,3 s/#/#\n#/" $var3
sed -i -e "1,11 s/\"double\"/double/" -e "1,3 s/#/#\n#/" $var4
sed -i -e "1,11 s/\"double\"/double/" -e "1,3 s/#/#\n#/" $var5


exit
