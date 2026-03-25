# Simulação Térmica 2D — Microprocessador com Hotspot

Projeto de Iniciação Científica: simulação numérica da equação do calor 2D não-estacionária aplicada a um microprocessador com hotspot gaussiano de atividade intermitente (duty-cycle). Três esquemas numéricos são implementados em Fortran e os resultados são analisados e visualizados em MATLAB.

---

## Estrutura do projeto

```
microprocessador_2D/
├── src/
│   ├── ftcs.f90              # Esquema explícito FTCS
│   ├── btcs.f90              # Backward Euler (implícito, Gauss-Seidel)
│   ├── crank_nicolson.f90    # Crank-Nicolson (implícito, Gauss-Seidel)
│   └── ftcs_instavel.f90     # FTCS instável (demonstração didática)
├── scripts/
│   ├── run_results.sh        # Roda os 3 métodos na malha 100x100
│   ├── run_convergence.sh    # Estudo de convergência espacial
│   ├── run_instability.sh    # Demonstração de instabilidade numérica
│   └── run_sensitivity.sh    # Estudo de sensibilidade paramétrica
├── plots/matlab/
│   ├── plot_comparacoes_analises.mlx   # Comparação dos 3 métodos
│   ├── teste_sensibilidade.mlx         # Análise de sensibilidade
├── data/
│   └── results/              # Arquivos .dat gerados pelas simulações
└── report/                   # Relatório técnico (em desenvolvimento)
```

---

## Formulação matemática

Equação do calor 2D com termo fonte:

$$\rho c_p \frac{\partial T}{\partial t} = k \left(\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2}\right) + Q(x, y, t)$$

**Hotspot gaussiano com duty-cycle:**

$$Q(x, y, t) = Q_0 \exp\left(-\frac{(x-x_0)^2 + (y-y_0)^2}{2\sigma^2}\right) \cdot S(t)$$

onde $S(t)$ é uma função que liga/desliga o hotspot periodicamente.

**Condições de contorno:**
- Laterais: adiabáticas (Neumann) — $\partial T/\partial n = 0$
- Base: temperatura fixa (Dirichlet) — $T = T_\infty$
- Topo: convecção (Robin) — $-k\,\partial T/\partial n = h(T - T_\infty)$

---

## Parâmetros físicos (silício)

| Parâmetro | Valor | Unidade |
|---|---|---|
| $L_x = L_y$ | 0.02 | m |
| $\rho$ | 2330 | kg/m³ |
| $c_p$ | 700 | J/(kg·K) |
| $k$ | 150 | W/(m·K) |
| $\alpha = k/(\rho c_p)$ | 9.20×10⁻⁵ | m²/s |
| $h$ | 10 | W/(m²·K) |
| $T_\infty$ | 27 | °C |
| $Q_0$ | 1×10⁸ | W/m³ |
| $\sigma$ | 0.002 | m |
| Período | 0.5 | s |
| $t_{max}$ | 4.75 | s |

---

## Esquemas numéricos implementados

| Esquema | Tipo | Ordem temporal | Estabilidade | Solver |
|---|---|---|---|---|
| FTCS | Explícito | 1ª | Condicional ($r_x + r_y \leq 0.5$) | Direto |
| Backward Euler | Implícito | 1ª | Incondicional | Gauss-Seidel |
| Crank-Nicolson | Implícito | 2ª | Incondicional | Gauss-Seidel |

---

## Como rodar

### Pré-requisitos
- `gfortran` instalado
- MATLAB (para visualização)
- WSL ou Linux

### 1. Resultados principais (malha 100×100)
```bash
cd scripts/
bash run_results.sh
```
Gera: `metricas_*.dat`, `log_*.txt`, `serie_tempo_*.dat`, `*_campo.dat`

### 2. Estudo de convergência espacial
```bash
bash run_convergence.sh
```
Roda malhas 5×5, 10×10, 20×20, 50×50, 100×100 para os 3 métodos.
Gera: `conv_ftcs.dat`, `conv_btcs.dat`, `conv_cn.dat`

### 3. Demonstração de instabilidade numérica
```bash
bash run_instability.sh
```
Testa FTCS com $\Delta t \times 1.02$ (violação da condição de Von Neumann).
Gera: `stab_ftcs_f1.02.dat`, `stab_btcs_f1.02.dat`, `stab_cn_f1.02.dat`

### 4. Estudo de sensibilidade
```bash
bash run_sensitivity.sh
```
Varia $Q_0$, posição do hotspot $x_0$ e duty-cycle para os 3 métodos.
Gera: `sens_*_<caso>.dat` e `sens_campo_*_<caso>.dat`

### Entrada dos programas Fortran
Todos os executáveis recebem 7 parâmetros via stdin:
```
Nx  Ny  fator_dt  Q0  x0  y0  duty
```
Exemplo:
```bash
echo "100 100 1.0 1e8 0.01 0.01 0.5" | ./sim_cn
```

---

## Visualização (MATLAB)

Abra o MATLAB na pasta `data/results/` e rode:

| Script | Conteúdo |
|---|---|
| `plot_comparacoes_analises.mlx` | Mapas de calor, superfície 3D, T_max(t), perfis, convergência, estabilidade, snapshots, erro vs Δx, CPU, tabela |
| `teste_sensibilidade.mlx` | Análise de sensibilidade: Q₀, x₀, duty-cycle |   

---

## Casos do estudo de sensibilidade

| Caso | Q₀ (W/m³) | x₀ (m) | Duty |
|---|---|---|---|
| `base` | 1.0×10⁸ | 0.010 | 0.50 |
| `Q0_05x` | 0.5×10⁸ | 0.010 | 0.50 |
| `Q0_2x` | 2.0×10⁸ | 0.010 | 0.50 |
| `x0_025` | 1.0×10⁸ | 0.005 | 0.50 |
| `x0_075` | 1.0×10⁸ | 0.015 | 0.50 |
| `duty_025` | 1.0×10⁸ | 0.010 | 0.25 |
| `duty_075` | 1.0×10⁸ | 0.010 | 0.75 |

---

## Resultados obtidos (malha 100×100)

| Esquema | T_max final (°C) | CPU (s) | Limiar 29 °C (s) | Iter/passo |
|---|---|---|---|---|
| FTCS | ~33.0 | ~1.0 | ~0.076 | — |
| Backward Euler | ~33.0 | ~3.3 | ~0.110 | ~1.03 |
| Crank-Nicolson | ~33.0 | ~4.5 | ~0.091 | ~1.06 |