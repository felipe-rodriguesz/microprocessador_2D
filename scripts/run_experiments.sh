#!/bin/bash
# run_experiments.sh

SRC="../src"
DATA="../data/results"

echo "========================================"
echo "  COMPILANDO OS TRÊS ESQUEMAS"
echo "========================================"
gfortran -Wall -O2 $SRC/ftcs.f90          -o sim_ftcs
gfortran -Wall -O2 $SRC/btcs.f90          -o sim_btcs
gfortran -Wall -O2 $SRC/crank_nicolson.f90 -o sim_cn

echo ""
echo "========================================"
echo "  RODANDO SIMULAÇÕES PRINCIPAIS"
echo "========================================"
./sim_ftcs
./sim_btcs
./sim_cn

# Move os .dat para a pasta correta
mv *.dat $DATA/

echo ""
echo "========================================"
echo "  ESTUDO DE CONVERGÊNCIA (CN)"
echo "========================================"
# Limpa dados anteriores
rm -f $DATA/dados_convergencia.dat

for malha in 5 10 20 50 100 200; do
    echo "  Rodando malha ${malha}x${malha}..."
    echo "$malha $malha" | ./sim_cn
    mv crank_nicolson.dat $DATA/cn_${malha}x${malha}.dat
done

mv dados_convergencia.dat $DATA/

echo ""
echo "Todos os experimentos concluídos!"
echo "Resultados em: $DATA/"