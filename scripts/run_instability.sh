#!/bin/bash
cd "$(dirname "$0")/.."

echo "==================================================="
echo "         INICIANDO ESTUDO DE ESTABILIDADE          "
echo "==================================================="

echo "[..] Compilando os programas..."
gfortran -Wall -O2 src/ftcs.f90           -o sim_ftcs || { echo "[ERRO] Falha ao compilar FTCS"; exit 1; }
gfortran -Wall -O2 src/btcs.f90           -o sim_btcs || { echo "[ERRO] Falha ao compilar BTCS"; exit 1; }
gfortran -Wall -O2 src/crank_nicolson.f90 -o sim_cn   || { echo "[ERRO] Falha ao compilar CN"; exit 1; }

echo "-> Executando simulacoes (Lendo de parametros.txt)..."
cp src/parametros.txt . || { echo "[ERRO] parametros.txt nao encontrado na pasta src/"; exit 1; }

MALHA=50
FATORES=(1.0 1.02)

for F in "${FATORES[@]}"; do
    echo "==================================================="
    echo " Testando Fator de Passo de Tempo (dt) = $F"

    # O SED substitui a linha 5 do parametros.txt com a malha fixa (50) e o novo Fator (F)
    sed -i "5s/.*/$MALHA $MALHA $F             ! 5. Malha e Tempo: Nx, Ny, fator_dt/" parametros.txt

    echo "  -> Rodando FTCS..."
    ./sim_ftcs > /dev/null
    mv serie_tempo_ftcs.dat data/results/stab_ftcs_f$F.dat 2>/dev/null

    echo "  -> Rodando Backward Euler..."
    ./sim_btcs > /dev/null
    mv serie_tempo_btcs.dat data/results/stab_btcs_f$F.dat 2>/dev/null

    echo "  -> Rodando Crank-Nicolson..."
    ./sim_cn > /dev/null
    mv serie_tempo_cn.dat data/results/stab_cn_f$F.dat 2>/dev/null
done

rm -f sim_ftcs sim_btcs sim_cn parametros.txt *.dat
echo "==================================================="
echo " TESTE FINALIZADO! Resultados em data/results/     "
echo "==================================================="