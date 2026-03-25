#!/bin/bash
cd "$(dirname "$0")/.."

echo "==================================================="
echo " GERANDO DADOS PARA A TABELA DO RELATÓRIO          "
echo "==================================================="

gfortran -Wall -O2 src/ftcs.f90           -o sim_ftcs || { echo "[ERRO] Falha ao compilar FTCS"; exit 1; }
gfortran -Wall -O2 src/btcs.f90           -o sim_btcs || { echo "[ERRO] Falha ao compilar BTCS"; exit 1; }
gfortran -Wall -O2 src/crank_nicolson.f90 -o sim_cn   || { echo "[ERRO] Falha ao compilar CN"; exit 1; }

echo "100 100 1.0 1.0e8 0.010 0.010 0.5" | ./sim_ftcs > /dev/null
echo "100 100 1.0 1.0e8 0.010 0.010 0.5" | ./sim_btcs > /dev/null
echo "100 100 1.0 1.0e8 0.010 0.010 0.5" | ./sim_cn   > /dev/null

# Move resultados para data/results/
mv log_*.txt metricas_*.dat serie_tempo_*.dat snap_*.dat *.dat data/results/ 2>/dev/null
rm -f sim_ftcs sim_btcs sim_cn

echo "PRONTO! Pode rodar o script no MATLAB."