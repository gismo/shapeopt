# This script runs test of Lagrangian and Winslow linearizations for shapeopt with IGA

OUT="output.txt";
BS="/";
BASE_FOLDER="../results/OptParam3D_affine/"

A=3; # quA
B=2; # quB

dim=3;


# Save the ipopt preferences
mkdir -p "${BASE_FOLDER}"

for R in 0  
do
    cp ipopt_optParam3D.opt ../gismo/filedata/options/ipopt.opt

    for e in 1 0.5 0.25 0.125 0.0625 0.03125
    do
        
       FOLDER="$BASE_FOLDER/WaterPassage/refine${R}_reg${e}/";
       
       mkdir -p "${FOLDER}";
       cp ../gismo/filedata/options/ipopt.opt $FOLDER/ipopt.opt
       
       (time ./main -r $R -e $e --glueInterfaces --usePow --optParamXML --startFile "parametrizations/XML/water_passage_fixed.xml" -d $dim -o "${FOLDER}"
       ) |& tee "${FOLDER}$BS$OUT"; 

    done


    cp ipopt_optParam3D_decrTau.opt ../gismo/filedata/options/ipopt.opt

    for e in 0.25 1 8 
    do
        
       FOLDER="$BASE_FOLDER/WaterPassage/refine${R}_decreasingTau${e}/";
       
       mkdir -p "${FOLDER}";
       cp ../gismo/filedata/options/ipopt.opt $FOLDER/ipopt.opt
       
       (time ./main -r $R -e $e --glueInterfaces --usePow--decreasingTau --decrTauFactor 0.25 --optParamXML --startFile "parametrizations/XML/water_passage_fixed.xml" -d $dim -o "${FOLDER}"
       ) |& tee "${FOLDER}$BS$OUT"; 
    done
done
