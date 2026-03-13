#!/bin/bash
cd "$(dirname "$0")/.."

echo "==================================================="
echo " GERANDO DADOS PARA A TABELA DO RELATÓRIO          "
echo "==================================================="

gfortran -Wall -O2 src/ftcs.f90 -o sim_ftcs
gfortran -Wall -O2 src/btcs.f90 -o sim_btcs
gfortran -Wall -O2 src/crank_nicolson.f90 -o sim_cn

# Roda os 3 métodos na malha 100x100 COM o fator_dt = 1.0
echo "100 100 1.0" | ./sim_ftcs > /dev/null
echo "100 100 1.0" | ./sim_btcs > /dev/null
echo "100 100 1.0" | ./sim_cn > /dev/null

# Move as métricas e logs para data/results/
mv log_*.txt metricas_*.dat serie_tempo_*.dat snap_*.dat *.dat data/results/ 2>/dev/null
rm -f sim_ftcs sim_btcs sim_cn

echo "PRONTO! Pode rodar o script no MATLAB."