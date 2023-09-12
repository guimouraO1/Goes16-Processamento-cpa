#!/bin/bash

# Ativar o ambiente Conda "goes"
source /root/miniconda3/bin/activate goes

# Navegar até a pasta "/Scripts/goes16/processamento/"
cd /Scripts/goes16/processamento/

# Inicia o download das imagens e depois processa
python download_amazon.py && python main.py

# Desativar o ambiente Conda quando você terminar, se desejar
conda deactivate