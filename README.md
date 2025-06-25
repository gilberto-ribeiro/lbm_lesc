# LBM-LESC

Simulador de Lattice Boltzmann Method (LBM) com suporte a 2D (D2Q9) e 3D (D3Q27), pós-processamento e paralelização via Rayon.

## Estrutura do Projeto

```
.
├── Cargo.toml
├── src/
│   ├── main.rs                # Ponto de entrada do programa, faz o parsing dos argumentos e executa as rotinas principais
│   ├── lib.rs                 # Biblioteca principal, organiza e expõe os módulos do projeto
│   ├── d2q9.rs                # Implementação do modelo LBM D2Q9 (2D, 9 velocidades)
│   ├── d3q27.rs               # Implementação do modelo LBM D3Q27 (3D, 27 velocidades)
│   ├── global_variables.rs    # Variáveis e constantes globais usadas em toda a simulação
│   ├── io.rs                  # Funções gerais de entrada e saída de dados
│   ├── post.rs                # Funções gerais de pós-processamento
│   ├── d2q9/
│   │   ├── bc.rs              # Condições de contorno específicas para D2Q9
│   │   ├── io.rs              # Entrada e saída de dados para D2Q9
│   │   ├── post.rs            # Pós-processamento para D2Q9
│   │   └── post/vtk.rs        # Exportação de resultados em VTK para D2Q9
│   └── d3q27/
│       ├── bc.rs              # Condições de contorno específicas para D3Q27
│       ├── io.rs              # Entrada e saída de dados para D3Q27
│       ├── post.rs            # Pós-processamento para D3Q27
│       └── post/vtk.rs        # Exportação de resultados em VTK para D3Q27
└── ...
```

- **main.rs**: Responsável por interpretar os argumentos da linha de comando e iniciar as simulações ou o pós-processamento.
- **lib.rs**: Organiza os módulos e funções compartilhadas do projeto.
- **d2q9.rs / d3q27.rs**: Implementam os métodos LBM para 2D e 3D, respectivamente.
- **global_variables.rs**: Define parâmetros globais, como dimensões, viscosidade, etc.
- **io.rs**: Funções para leitura e escrita de arquivos comuns a ambos os modelos.
- **post.rs**: Rotinas de pós-processamento comuns.
- **d2q9/** e **d3q27/**: Contêm módulos específicos para cada modelo, incluindo condições de contorno, IO e pós-processamento.
- **post/vtk.rs**: Exporta os resultados das simulações para o formato VTK, facilitando a visualização em softwares científicos.

Essas divisões facilitam a manutenção, extensão e reutilização do código para diferentes tipos de simulação e pós-processamento.

## Instalação

Requer [Rust](https://www.rust-lang.org/tools/install). Certifique-se de que o Rust está instalado no seu sistema antes de prosseguir. Para instalar o Rust, siga os passos abaixo:

1. Baixe e execute o instalador do Rust com o comando:
    ```sh
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
    ```

2. Siga as instruções exibidas no terminal para concluir a instalação.

3. Após a instalação, adicione o diretório do Rust ao seu PATH (se necessário) e verifique se ele foi instalado corretamente:
    ```sh
    rustc --version
    ```

    Este comando deve exibir a versão do compilador Rust instalada.

Depois de instalar o Rust, clone o repositório do projeto e entre no diretório `./lbm_lesc`:

```sh
mkdir lbm_lesc
git clone https://github.com/gilberto-ribeiro/lbm_lesc.git lbm_lesc
cd lbm_lesc
```

## Tutorial de Utilização

### 1. Compilando o Projeto

Compile o executável com o comando:

```sh
cargo build --release
```

O executável será gerado em `target/release/lbm_lesc`.

---

### 2. Preparando o Diretório de Simulação

Crie (ou entre) em um diretório onde estão os arquivos de configuração da sua simulação (por exemplo, `./cases/case_001`):

```sh
cd ./cases/case_001
```

Copie o executável para este diretório:

```sh
cp ../../target/release/lbm_lesc .
```

---

### 3. Executando a Simulação

Para rodar a simulação, utilize o executável `lbm_lesc` com os parâmetros desejados:

```sh
./lbm_lesc -d <DIMENSIONS> -n <NUMBER_OF_THREADS> run
```

- `-d, --dimensions`: Define a dimensão da simulação (`2` para 2D ou `3` para 3D).
- `-n, --number_of_threads`: Define o número de threads para paralelização (`1`, `2`, `4`, `8`, `16` ou `32`).
- `run`: Subcomando para executar a simulação.

**Exemplo:**
```sh
./lbm_lesc -d 3 -n 8 run
```

#### Executando em modo benchmark

Para rodar a simulação em modo benchmark, adicione a flag `-b` ou `--benchmark`:

```sh
./lbm_lesc -d 3 -n 8 run --benchmark
```

---

### 4. Pós-processamento (Geração de arquivos VTK)

Após a simulação, para gerar arquivos VTK para visualização e análise, utilize o subcomando `post`:

```sh
./lbm_lesc -d <DIMENSIONS> -n <NUMBER_OF_THREADS> post
```

**Exemplo:**
```sh
./lbm_lesc -d 2 -n 4 post
```

---

### 5. Saídas

- Os resultados e arquivos de pós-processamento são salvos no diretório `./post_processing` e subdiretórios.
- Arquivos `.dat` e `.vtk` são gerados para análise e visualização.

---

### Observações

- Certifique-se de que os arquivos de configuração necessários estejam presentes no diretório de simulação antes de executar o programa.
- O número de threads pode ser ajustado conforme o hardware disponível para otimizar o desempenho.
- Para mais detalhes sobre os parâmetros e opções, execute:
  ```sh
  ./lbm_lesc --help
  ```

## Casos Tutoriais

O diretório `./cases` contém exemplos práticos para diferentes tipos de escoamento, servindo como tutoriais para uso do simulador.  
Cada caso possui um diretório `pre_processing` com dois arquivos principais:
- `case_conditions.jou`: Define os parâmetros de simulação, condições iniciais e de contorno.
- `case_setup.jou`: Define configurações da simulação, como tolerâncias e número de iterações.

### Lista de Casos

- **case_001:** Escoamento de Couette usando condições de contorno do tipo bounce-back.
- **case_002:** Escoamento de Couette usando condições de contorno do tipo Zou e He (1997).
- **case_003:** Escoamento de Poiseuille usando condições de contorno do tipo Zou e He (1997).
- **case_004:** Escoamento de Poiseuille usando condições de contorno do tipo bounce-back com velocidade de entrada prescrita.
- **case_005:** Escoamento de Poiseuille usando condições de contorno do tipo Zou e He (1997) com velocidade de entrada prescrita.
- **case_006:** Escoamento de esteira de Von Karman usando condições de contorno do tipo Zou e He (1997).
- **case_007:** Escoamento de cavidade com tampa móvel (lid driven cavity) usando condições de contorno do tipo bounce-back.
- **case_008:** Escoamento de cavidade com tampa móvel (lid driven cavity) usando condições de contorno do tipo Zou e He (1997).
- **case_009:** Escoamento de cavidade com tampa móvel (lid driven cavity) tridimensional usando condições de contorno do tipo bounce-back.
- **case_010:** Escoamento em estrutura TPMS usando condições de contorno do tipo bounce-back (não deslizamento nas faces laterais).
- **case_011:** Escoamento em estrutura TPMS usando condições de contorno do tipo bounce-back (periodicidade nas faces laterais).

Esses casos cobrem diferentes geometrias, condições de contorno e abordagens clássicas de simulação, sendo úteis para validação, aprendizado e comparação

## Licença

MIT

---

> Para detalhes de implementação, consulte os arquivos-fonte em [src/](src).

quero que descreva brevemente os casos tutoriais que estão dentro de ./cases. Os detalhes dos parametros de simulação, condições iniciais, de contorno estão em case_conditions.jou e detalhes da simulação estão em case_setup.jou, todos no diretório pre_processing de cada caso

case_001: Escoamento de couette usando condições de contorno do tipo bounce-back
case_002: Escoamento de couette usando condições de contorno do tipo Zou e He (1997)
case_003: Escoamento de poiseuille usando condições de contorno do tipo Zou e He (1997)
case_004: Escoamento de poiseuille usando condições de contorno do tipo bounce-back tendo a velocidade de entrada prescrita
case_005: Escoamento de poiseuille usando condições de contorno do tipo Zou e He (1997) tendo a velocidade de entrada prescrita
case_006: Escoamento de esteira de Von-karman usando condições de contorno do tipo Zou e He (1997)
case_007: Escoamento de cavidade com tampa móvel (lid driven cavity) usando condições de contorno do tipo bounce-back
case_008: Escoamento de cavidade com tampa móvel (lid driven cavity) usando condições de contorno do tipo Zou e He (1997)
case_009: Escoamento de cavidade com tampa móvel (lid driven cavity) tridimensional usando condições de contorno do tipo bounce-back
case_010: Escoamento em estrutura TPMS usando condições de contorno do tipo bounce-back (Não deslizamento nas faces laterais)
case_011: Escoamento em estrutura TPMS usando condições de contorno do tipo bounce-back (Periodicidade nas faces laterais)
