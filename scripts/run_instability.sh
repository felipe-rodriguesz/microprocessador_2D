#!/bin/bash
cd "$(dirname "$0")/.."

echo "==================================================="
echo "         INICIANDO ESTUDO DE ESTABILIDADE          "
echo "==================================================="

echo "[..] Compilando os programas..."
gfortran -Wall -O2 src/ftcs.f90           -o sim_ftcs
gfortran -Wall -O2 src/btcs.f90           -o sim_btcs
gfortran -Wall -O2 src/crank_nicolson.f90 -o sim_cn

if [ $? -ne 0 ]; then
    echo "[ERRO] Falha na compilacao."
    exit 1
fi

MALHA=50
FATORES=(1.0 1.02)

for F in "${FATORES[@]}"; do
    echo "==================================================="
    echo " Testando Fator de Passo de Tempo (dt) = $F"

    echo "$MALHA $MALHA $F 1.0e8 0.010 0.010 0.5" | ./sim_ftcs > /dev/null
    mv serie_tempo_ftcs.dat data/results/stab_ftcs_f$F.dat 2>/dev/null

    echo "$MALHA $MALHA $F 1.0e8 0.010 0.010 0.5" | ./sim_btcs > /dev/null
    mv serie_tempo_btcs.dat data/results/stab_btcs_f$F.dat 2>/dev/null

    echo "$MALHA $MALHA $F 1.0e8 0.010 0.010 0.5" | ./sim_cn > /dev/null
    mv serie_tempo_cn.dat data/results/stab_cn_f$F.dat 2>/dev/null
done

rm -f sim_ftcs sim_btcs sim_cn *.dat
echo "==================================================="
echo " TESTE FINALIZADO! Resultados em data/results/     "
echo "==================================================="