# This script runs test of Lagrangian and Winslow linearizations for shapeopt with IGA

OUT="output.txt";
BS="/";
BASE_FOLDER="../results/OptAntenna_startGuesses"

R=4; # numRefine
I=9; # maxiter
A=3; # quA
B=1; # quB

N=5;
M=4;

e=0.125;

# Save the ipopt preferences
#cp ../gismo/filedata/options/ipopt.opt $BASE_FOLDER/ipopt.opt
#
#FOLDER="$BASE_FOLDER/refine${R}_reg${e}_fromSquare/";
#
#mkdir -p "${FOLDER}"
#
#(time ./main -n $N -m $M -r $R -A $A -B $B -i $I -a 5 -e $e -o "$BS${FOLDER}" --startFromFile -f "startGuess/square.txt") |& tee "${FOLDER}$BS$OUT"; 

FOLDER="$BASE_FOLDER/refine${R}_winslow44_fromSquare/";

mkdir -p "${FOLDER}"

(time ./main -n $N -m $M -r $R -A $A -B $B -i $I -a 6 -o "$BS${FOLDER}" --startFromFile -f "startGuess/square.txt") |& tee "${FOLDER}$BS$OUT"; 

#FOLDER="$BASE_FOLDER/refine${R}_reg${e}_fromLin/";
#
#mkdir -p "${FOLDER}"
#
#(time ./main -n $N -m $M -r $R -A $A -B $B -i $I -a 5 -e $e -o "$BS${FOLDER}" --startFromFile -f "startGuess/ref4_winslow44_5.txt") |& tee "${FOLDER}$BS$OUT"; 
#
#FOLDER="$BASE_FOLDER/refine${R}_winslow44_fromReg/";
#
#mkdir -p "${FOLDER}"
#
#(time ./main -n $N -m $M -r $R -A $A -B $B -i $I -a 6 -s 20 -o "$BS${FOLDER}" --startFromFile -f "startGuess/ref4_reg0125.txt") |& tee "${FOLDER}$BS$OUT"; 
