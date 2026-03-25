#!/bin/bash
cd "$(dirname "$0")/.."

echo "==================================================="
echo "         ESTUDO DE SENSIBILIDADE                   "
echo "==================================================="
echo " Parametros variados:"
echo "   - Posicao do hotspot (x0)"
echo "   - Intensidade Q0"
echo "   - Duty-cycle"
echo "==================================================="

# Compilação
echo "[..] Compilando os programas..."
gfortran -Wall -O2 src/ftcs.f90           -o sim_ftcs || { echo "[ERRO] Falha ao compilar FTCS"; exit 1; }
gfortran -Wall -O2 src/btcs.f90           -o sim_btcs || { echo "[ERRO] Falha ao compilar BTCS"; exit 1; }
gfortran -Wall -O2 src/crank_nicolson.f90 -o sim_cn   || { echo "[ERRO] Falha ao compilar CN"; exit 1; }
echo "[OK] Compilacao concluida."
echo ""

# Funcao auxiliar para rodar os 3 metodos e salvar
rodar_caso() {
    local CASO=$1
    local PARAMS=$2

    echo "  -> Rodando caso: $CASO  (params: $PARAMS)"
    echo "$PARAMS" | ./sim_ftcs > /dev/null
    mv serie_tempo_ftcs.dat data/results/sens_ftcs_${CASO}.dat 2>/dev/null
    mv ftcs.dat             data/results/sens_campo_ftcs_${CASO}.dat 2>/dev/null

    echo "$PARAMS" | ./sim_btcs > /dev/null
    mv serie_tempo_btcs.dat data/results/sens_btcs_${CASO}.dat 2>/dev/null
    mv backward_euler.dat   data/results/sens_campo_btcs_${CASO}.dat 2>/dev/null

    echo "$PARAMS" | ./sim_cn > /dev/null
    mv serie_tempo_cn.dat   data/results/sens_cn_${CASO}.dat 2>/dev/null
    mv crank_nicolson.dat   data/results/sens_campo_cn_${CASO}.dat 2>/dev/null
}

# ===========================================================
# CASO BASE (referência)
# Nx Ny fator_dt Q0    x0     y0     duty
# ===========================================================
echo "==================================================="
echo " 1/4: Caso Base (referência)"
rodar_caso "base"     "100 100 1.0 1.0e8 0.010 0.010 0.5"

# ===========================================================
# VARIAÇÃO 1 — Posição do Hotspot
# ===========================================================
echo "==================================================="
echo " 2/4: Variacao da Posicao do Hotspot (x0)"
rodar_caso "x0_025"   "100 100 1.0 1.0e8 0.005 0.010 0.5"   # x0 = Lx/4
rodar_caso "x0_075"   "100 100 1.0 1.0e8 0.015 0.010 0.5"   # x0 = 3*Lx/4

# ===========================================================
# VARIAÇÃO 2 — Intensidade Q0
# ===========================================================
echo "==================================================="
echo " 3/4: Variacao da Intensidade Q0"
rodar_caso "Q0_05x"   "100 100 1.0 0.5e8 0.010 0.010 0.5"   # Q0 / 2
rodar_caso "Q0_2x"    "100 100 1.0 2.0e8 0.010 0.010 0.5"   # Q0 * 2

# ===========================================================
# VARIAÇÃO 3 — Duty-Cycle
# ===========================================================
echo "==================================================="
echo " 4/4: Variacao do Duty-Cycle"
rodar_caso "duty_025" "100 100 1.0 1.0e8 0.010 0.010 0.25"  # 25%
rodar_caso "duty_075" "100 100 1.0 1.0e8 0.010 0.010 0.75"  # 75%

# Limpeza
mv *.dat data/results/ 2>/dev/null
rm -f sim_ftcs sim_btcs sim_cn

echo ""
echo "==================================================="
echo " CONCLUIDO! Resultados em data/results/            "
echo " Arquivos gerados:                                 "
echo "   sens_ftcs_*.dat  sens_btcs_*.dat  sens_cn_*.dat "
echo "   sens_campo_ftcs_*.dat  ...                      "
echo "==================================================="