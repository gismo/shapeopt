# This script runs test of Lagrangian and Winslow linearizations for shapeopt with IGA

OUT="output.txt";
BS="/";
BASE_FOLDER="../results/OptAntenna_quad_2nd/"

mkdir -p "${BASE_FOLDER}"

R=4; # numRefine
I=9; # maxiter
A=3; # quA
B=1; # quB

N=5;
M=4;

e=$1;

qA=$2
qB=4;

# Save the ipopt preferences
cp ../gismo/filedata/options/ipopt.opt $BASE_FOLDER/ipopt.opt

FOLDER="$BASE_FOLDER/refine${R}_reg${e}_quadWinslow${qA}${qB}/";

mkdir -p "${FOLDER}"

(time ./main -n $N -m $M -r $R -A $A -B $B -i $I -a 5 -e $e --quA_optParam $qA --quB_optParam $qB -o "$BS${FOLDER}" --startFromFile -f "startGuess/circle.txt") |& tee "${FOLDER}$BS$OUT"  ;

