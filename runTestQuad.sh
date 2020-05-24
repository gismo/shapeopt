# This script

OUT="output.txt";
BS="/";
BASE_FOLDER="../results/OptAntenna_quad"

R=4; # numRefine
I=9; # maxiter
A=3; # quA
B=1; # quB

N=5;
M=4;

e=0.125;

# Save the ipopt preferences
cp ../gismo/filedata/options/ipopt.opt $BASE_FOLDER/ipopt.opt

for Q in 0 2 6
do    

    FOLDER="$BASE_FOLDER/refine${R}_reg${e}_quadWinslow${Q}4/";
    
    mkdir -p "${FOLDER}"
    
    (time ./main -n $N -m $M -r $R -A $A -B $B -i $I -a 5 -e $e --quA_optParam $Q --quB_optParam 4 "$BS${FOLDER}" --startFromFile -f "startGuess/circle.txt") |& tee "${FOLDER}$BS$OUT"; 
    
    FOLDER="$BASE_FOLDER/refine${R}_winslow${Q}4/";
    
    mkdir -p "${FOLDER}"
    
    (time ./main -n $N -m $M -r $R -A $A -B $B -i $I -a 6 --quA_optParam $Q --quB_optParam 4 -o "$BS${FOLDER}" --startFromFile -f "startGuess/circle.txt") |& tee "${FOLDER}$BS$OUT"; 

done

