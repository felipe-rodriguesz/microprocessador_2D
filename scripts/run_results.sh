#!/bin/bash
cd "$(dirname "$0")/.."

echo "==================================================="
echo " GERANDO DADOS PARA A TABELA DO RELATÓRIO          "
echo "==================================================="

echo "-> Compilando os simuladores..."
gfortran -Wall -O2 src/ftcs.f90           -o sim_ftcs || { echo "[ERRO] Falha ao compilar FTCS"; exit 1; }
gfortran -Wall -O2 src/btcs.f90           -o sim_btcs || { echo "[ERRO] Falha ao compilar BTCS"; exit 1; }
gfortran -Wall -O2 src/crank_nicolson.f90 -o sim_cn   || { echo "[ERRO] Falha ao compilar CN"; exit 1; }

echo "-> Executando simulacoes (Lendo de parametros.txt)..."
cp src/parametros.txt . || { echo "[ERRO] parametros.txt nao encontrado na pasta src/"; exit 1; }
./sim_ftcs > /dev/null
./sim_btcs > /dev/null
./sim_cn   > /dev/null

echo "-> Movendo resultados para data/results/..."
# Move resultados para data/results/
mv log_*.txt metricas_*.dat serie_tempo_*.dat snap_*.dat *.dat data/results/ 2>/dev/null
rm -f sim_ftcs sim_btcs sim_cn parametros.txt

echo "==================================================="
echo " PRONTO! Pode rodar o script no MATLAB."
echo "==================================================="