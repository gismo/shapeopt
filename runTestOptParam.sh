# This script runs test of Lagrangian and Winslow linearizations for shapeopt with IGA

OUT="output.txt";
BS="/";
BASE_FOLDER="../results/OptParam2D/"

A=3; # quA
B=2; # quB

dim=2;


# Save the ipopt preferences
mkdir -p "${BASE_FOLDER}"
cp ../gismo/filedata/options/ipopt.opt $BASE_FOLDER/ipopt_decrTau.opt

for R in 0 1  
do
    for e in 64 
    do
        for jig in 1 2
        do
        
        	#FOLDER="$BASE_FOLDER/Jigsaw${jig}/refine${R}_reg${e}/";
        	#
        	#mkdir -p "${FOLDER}";
        	#
        	#(time ./main -r $R -e $e --optParamXML --startFile "parametrizations/XML/2D/Jigsaw${jig}.xml" -d $dim -o "${FOLDER}"
            #) |& tee "${FOLDER}$BS$OUT"; 

        	FOLDER="$BASE_FOLDER/Jigsaw${jig}/refine${R}_decreasingTau${e}/";
        	
        	mkdir -p "${FOLDER}";
        	
        	(time ./main -r $R -e $e --decreasingTau --decrTauFactor 0.25 --optParamXML --startFile "parametrizations/XML/2D/Jigsaw${jig}.xml" -d $dim -o "${FOLDER}"
            ) |& tee "${FOLDER}$BS$OUT"; 
        done
    done
done
