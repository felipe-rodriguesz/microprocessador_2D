# Simulacao Termica 2D - Microprocessador (Chip)

Projeto de simulacao numerica em Fortran para problemas 2D de transferencia de calor: chip transiente com hotspot gaussiano, chip em regime estacionario, Laplace e Poisson com solucao manufaturada. Os resultados sao salvos em `data/results/` e podem ser visualizados no MATLAB.

---

## Estrutura do projeto

```
microprocessador_2D/
├── src/
│   ├── simulador_chip.f90         # Transiente com hotspot gaussiano + duty-cycle
│   ├── chip_estacionario.f90      # Regime estacionario (chip + conveccao)
│   ├── simulador_laplace.f90      # Laplace 2D com CC de Dirichlet
│   └── simulador_poisson.f90      # Poisson 2D (MMS sin/cos)
├── parametros.txt                 # Entrada do simulador_chip
├── parametros_estacionario.txt    # Entrada do chip_estacionario
├── parametros_laplace.txt         # Entrada do simulador_laplace
├── parametros_poisson.txt         # Entrada do simulador_poisson
├── data/
│   └── results/                   # Arquivos .dat gerados pelas simulacoes
└── plots/matlab/
    ├── plot_chip.mlx
    ├── plot_laplace.mlx
    └── plot_poisson.mlx
```

---

## Modelos e condicoes de contorno

- **Chip transiente e estacionario**: equacao do calor 2D com fonte gaussiana. Contornos: base Dirichlet (T0), laterais adiabaticas (Neumann) e topo convectivo (Robin).
- **Laplace**: Dirichlet com topo a 100 e demais bordas a 0.
- **Poisson**: solucao manufaturada $T(x,y)=sin(pi x) cos(pi y)$ para validacao do erro.

---

## Como rodar

### Pre-requisitos
- `gfortran` instalado
- MATLAB (para visualizacao)
- Linux ou WSL

### 1. Chip transiente (hotspot + duty-cycle)
```bash
cd src/
gfortran -O2 -std=f2008 simulador_chip.f90 -o simulador_chip
./simulador_chip
```
Entradas em `parametros.txt` (namelist `/parametros/`).

Saidas em `data/results/`:
- `norma_iter_chip.dat`
- `convergencia_passos.dat`
- `tmax_tempo.dat`
- `sensores_tempo.dat`
- `mapa_calor_chip.dat`

### 2. Chip estacionario
```bash
cd src/
gfortran -O2 -std=f2008 chip_estacionario.f90 -o chip_estacionario
./chip_estacionario
```
Entradas em `parametros_estacionario.txt`.

Saidas em `data/results/`:
- `norma_iter_chip.dat`
- `mapa_calor_chip.dat`

### 3. Laplace 2D
```bash
cd src/
gfortran -O2 -std=f2008 simulador_laplace.f90 -o simulador_laplace
./simulador_laplace
```
Entradas em `parametros_laplace.txt`.

Saidas em `data/results/`:
- `norma_iter_laplace.dat`
- `mapa_calor_laplace.dat`

### 4. Poisson 2D (MMS)
```bash
cd src/
gfortran -O2 -std=f2008 simulador_poisson.f90 -o simulador_poisson
./simulador_poisson
```
Entradas em `parametros_poisson.txt`.

Saidas em `data/results/`:
- `norma_iter_poisson.dat`
- `mapa_calor_poisson.dat`

---

## Arquivos de parametros

Os programas leem um namelist `/parametros/` com os campos de cada arquivo `parametros_*.txt`. Exemplo (transiente):

```
&parametros
 Lx = 0.02d0,
 Ly = 0.02d0,
 Nx = 100,
 Ny = 100,
 k = 150.0d0,
 rho = 2330.0d0,
 cp = 700.0d0,
 Q0 = 1.0d8,
 x0 = 0.01d0,
 y0 = 0.01d0,
 sigma = 0.002d0,
 T0 = 27.0d0,
 Tinf = 27.0d0,
 h = 10.0d0,
 dt = 0.005d0,
 tmax = 5.0d0,
 t_ciclo_ini = 0.0d0,
 t_ciclo_fim = 5.0d0,
 period = 0.050d0,
 duty = 0.50d0,
 tol = 1.0d-8,
 max_iter = 100000,
/
```

Notas:
- O arquivo deve terminar com `/` **e** nova linha ao final (o `gfortran` pode falhar sem a nova linha).
- O simulador transiente usa `t_ciclo_ini`, `t_ciclo_fim`, `period` e `duty` para ligar/desligar a fonte.

---

## Visualizacao (MATLAB)

Abra o MATLAB na pasta `data/results/` e rode:

| Script | Conteudo |
|---|---|
| `plot_chip.mlx` | Mapas de calor, evolucao temporal e sensores do chip |
| `plot_laplace.mlx` | Campo de temperatura de Laplace |
| `plot_poisson.mlx` | Campo e erro maximo do teste de Poisson |