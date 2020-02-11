# This script runs test of Lagrangian and Winslow linearizations for shapeopt with IGA

OUT="output.txt";
BS="/";
BASE_FOLDER="../results/OptParam3D_2"

#R=0; # numRefine

A=3; # quA
B=3; # quB

#e=100;

# Save the ipopt preferences
# cp ../gismo/filedata/options/ipopt.opt $BASE_FOLDER/ipopt.opt

for R in {0,1,2}
do
	for e in {1000,100,10,1,0.1}
	do
		FOLDER="$BASE_FOLDER/refine${R}_reg${e}/";
		
		mkdir -p "${FOLDER}"
		
		(time ./main -r $R -A $A -B $B -e $e -o "$BS${FOLDER}") |& tee "${FOLDER}$BS$OUT"; 
	done

done
