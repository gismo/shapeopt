# This script runs test of Lagrangian and Winslow linearizations for shapeopt with IGA

# Function to run code
runCode () {
        mkdir -p $3;
        cp ../gismo/filedata/options/ipopt.opt $3/ipopt.opt

        echo $1
        echo $2
        echo $3
        echo $4

        (time ./main -r $1 -a $2 --paramTestXML --startFile $4 -o $3) |& tee "${FOLDER}$BS$OUT"; 
}

OUT="output.txt";
BS="/";
BASE_FOLDER="../results/ParamTests/"

A=3; # quA
B=2; # quB

# Save the ipopt preferences
mkdir -p "${BASE_FOLDER}"

for R in 1 
do
    for jig in 1 
    do
        #FOLDER="$BASE_FOLDER/Jigsaw${jig}_refine${R}/Spring/";
        #runCode $R 0 $FOLDER "parametrizations/XML/2D/Jigsaw${jig}.xml"

        #FOLDER="$BASE_FOLDER/Jigsaw${jig}_refine${R}/ModLiao/";
        #runCode $R 1 $FOLDER "parametrizations/XML/2D/Jigsaw${jig}.xml"

        #FOLDER="$BASE_FOLDER/Jigsaw${jig}_refine${R}/Winslow/";
        #runCode $R 2 $FOLDER "parametrizations/XML/2D/Jigsaw${jig}.xml"

        #FOLDER="$BASE_FOLDER/Jigsaw${jig}_refine${R}/Liao/";
        #runCode $R 3 $FOLDER "parametrizations/XML/2D/Jigsaw${jig}.xml"

        #FOLDER="$BASE_FOLDER/Jigsaw${jig}_refine${R}/Harmonic/";
        #runCode $R 4 $FOLDER "parametrizations/XML/2D/Jigsaw${jig}.xml"
        
        FOLDER="$BASE_FOLDER/Jigsaw${jig}_refine${R}/WinslowNoConst/";
        runCode $R 6 $FOLDER "parametrizations/XML/2D/Jigsaw${jig}.xml"

        FOLDER="$BASE_FOLDER/Jigsaw${jig}_refine${R}/MaxDetJac/";
        runCode $R 7 $FOLDER "parametrizations/XML/2D/Jigsaw${jig}.xml"

        FOLDER="$BASE_FOLDER/Jigsaw${jig}_refine${R}/HarmonicNoConst/";
        runCode $R 8 $FOLDER "parametrizations/XML/2D/Jigsaw${jig}.xml"
        
    done
done

#R=0
#FOLDER="$BASE_FOLDER/Seastar_refine${R}/Spring/";
#runCode $R 0 $FOLDER "parametrizations/XML/2D/Seastar.xml"
#
#FOLDER="$BASE_FOLDER/Seastar_refine${R}/ModLiao/";
#runCode $R 1 $FOLDER "parametrizations/XML/2D/Seastar.xml"
#
#FOLDER="$BASE_FOLDER/Seastar_refine${R}/Winslow/";
#runCode $R 2 $FOLDER "parametrizations/XML/2D/Seastar.xml"
#
#FOLDER="$BASE_FOLDER/Seastar_refine${R}/Liao/";
#runCode $R 3 $FOLDER "parametrizations/XML/2D/Seastar.xml"
#
#FOLDER="$BASE_FOLDER/Seastar_refine${R}/Harmonic/";
#runCode $R 4 $FOLDER "parametrizations/XML/2D/Seastar.xml"
#
#FOLDER="$BASE_FOLDER/Seastar_refine${R}/WinslowNoConst/";
#runCode $R 6 $FOLDER "parametrizations/XML/2D/Seastar.xml"
#
#FOLDER="$BASE_FOLDER/Seastar_refine${R}/MaxDetJac/";
#runCode $R 7 $FOLDER "parametrizations/XML/2D/Seastar.xml"
#
#FOLDER="$BASE_FOLDER/Seastar_refine${R}/HarmonicNoConst/";
#runCode $R 8 $FOLDER "parametrizations/XML/2D/Seastar.xml"
