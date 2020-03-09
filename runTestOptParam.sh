# This script runs test of Lagrangian and Winslow linearizations for shapeopt with IGA

OUT="output.txt";
BS="/";
BASE_FOLDER="../results/OptParamTests/jig3d_deg4/"

A=3; # quA
B=1; # quB

e=0.05;
jig=1;
dim=3;


# Save the ipopt preferences
mkdir -p "${BASE_FOLDER}"
cp ../gismo/filedata/options/ipopt.opt $BASE_FOLDER/ipopt.opt

for R in 1 
do
	FOLDER="$BASE_FOLDER/refine${R}_reg${e}_2nd/";
	
	mkdir -p "${FOLDER}"
	
	(time ./main --optParam -r ${R} -A $A -B $B -e $e -j $jig -d $dim -o "$BS${FOLDER}") |& tee "${FOLDER}$BS$OUT"; 
	
done
