#!/bin/bash

# Descobre a pasta do script e sobe um nível para a raiz do projeto
cd "$(dirname "$0")/.."

echo "==================================================="
echo "         INICIANDO ESTUDO DE CONVERGÊNCIA          "
echo "==================================================="

# Limpa arquivos velhos para evitar lixo de execuções passadas
rm -f data/results/dados_convergencia.dat data/results/conv_ftcs.dat data/results/conv_btcs.dat data/results/conv_cn.dat
echo "[OK] Arquivos antigos removidos de data/results/."

# Compila os três códigos Fortran
echo "[..] Compilando os programas..."
gfortran -Wall -O2 src/ftcs.f90 -o sim_ftcs
gfortran -Wall -O2 src/btcs.f90 -o sim_btcs
gfortran -Wall -O2 src/crank_nicolson.f90 -o sim_cn

# Verifica se a compilação deu certo antes de continuar
if [ $? -ne 0 ]; then
    echo "[ERRO] Falha na compilacao. Verifique os caminhos."
    exit 1
fi

# Define as malhas que você quer testar
MALHAS=(5 10 20 50 100)

echo "==================================================="
echo " 1/3: Rodando FTCS (Explícito)..."
for N in "${MALHAS[@]}"; do
    echo "  -> Malha $N x $N"
    echo "$N $N 1.0" | ./sim_ftcs > /dev/null
done
# Renomeia o arquivo gerado para não ser sobrescrito
mv dados_convergencia.dat data/results/conv_ftcs.dat

echo "==================================================="
echo " 2/3: Rodando Backward Euler (Implícito)..."
for N in "${MALHAS[@]}"; do
    echo "  -> Malha $N x $N"
    echo "$N $N 1.0" | ./sim_btcs > /dev/null
done
mv dados_convergencia.dat data/results/conv_btcs.dat

echo "==================================================="
echo " 3/3: Rodando Crank-Nicolson (Implícito)..."
for N in "${MALHAS[@]}"; do
    echo "  -> Malha $N x $N"
    echo "$N $N 1.0" | ./sim_cn > /dev/null
done
mv dados_convergencia.dat data/results/conv_cn.dat

# Move os outros .dat gerados e apaga os executáveis
mv *.dat data/results/ 2>/dev/null
rm -f sim_ftcs sim_btcs sim_cn

echo "==================================================="
echo " AUTOMACAO FINALIZADA COM SUCESSO!                 "
echo " Todos os resultados estao em: data/results/       "
echo "==================================================="