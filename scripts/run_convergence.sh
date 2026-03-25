#!/bin/bash
cd "$(dirname "$0")/.."

echo "==================================================="
echo "         INICIANDO ESTUDO DE CONVERGÊNCIA          "
echo "==================================================="

rm -f data/results/dados_convergencia.dat data/results/conv_ftcs.dat data/results/conv_btcs.dat data/results/conv_cn.dat
echo "[OK] Arquivos antigos removidos de data/results/."

echo "[..] Compilando os programas..."
gfortran -Wall -O2 src/ftcs.f90           -o sim_ftcs || { echo "[ERRO] Falha ao compilar FTCS"; exit 1; }
gfortran -Wall -O2 src/btcs.f90           -o sim_btcs || { echo "[ERRO] Falha ao compilar BTCS"; exit 1; }
gfortran -Wall -O2 src/crank_nicolson.f90 -o sim_cn   || { echo "[ERRO] Falha ao compilar CN"; exit 1; }

# Malhas testadas
MALHAS=(5 10 20 50 100)

# Nx Ny fator_dt Q0    x0     y0     duty  (parâmetros base)
echo "==================================================="
echo " 1/3: Rodando FTCS (Explícito)..."
for N in "${MALHAS[@]}"; do
    echo "  -> Malha $N x $N"
    echo "$N $N 1.0 1.0e8 0.010 0.010 0.5" | ./sim_ftcs > /dev/null
done
mv dados_convergencia.dat data/results/conv_ftcs.dat

echo "==================================================="
echo " 2/3: Rodando Backward Euler (Implícito)..."
for N in "${MALHAS[@]}"; do
    echo "  -> Malha $N x $N"
    echo "$N $N 1.0 1.0e8 0.010 0.010 0.5" | ./sim_btcs > /dev/null
done
mv dados_convergencia.dat data/results/conv_btcs.dat

echo "==================================================="
echo " 3/3: Rodando Crank-Nicolson (Implícito)..."
for N in "${MALHAS[@]}"; do
    echo "  -> Malha $N x $N"
    echo "$N $N 1.0 1.0e8 0.010 0.010 0.5" | ./sim_cn > /dev/null
done
mv dados_convergencia.dat data/results/conv_cn.dat

mv *.dat data/results/ 2>/dev/null
rm -f sim_ftcs sim_btcs sim_cn

echo "==================================================="
echo " AUTOMACAO FINALIZADA COM SUCESSO!                 "
echo " Todos os resultados estao em: data/results/       "
echo "==================================================="