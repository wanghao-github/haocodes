#!/bin/bash
for i in 0.1 0.3 0.5 0.7 0.9
do
	directoryname="U_"$i""
	mkdir $directoryname	
	cp INCAR POSCAR KPOINTS POTCAR subjob_2022 $directoryname
	cd $directoryname
	sed -i "s/LDAUU =  0 4 0/LDAUU =  0 $i 0/" INCAR
	sbatch subjob_2022	
	cd ../
done
