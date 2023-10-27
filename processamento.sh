#!/bin/bash

# Ativar o ambiente Conda "goes"
source /root/miniconda3/bin/activate goes16

# Navegar até a pasta "/Scripts/goes16/processamento/"
cd /Scripts/goes16/Processamento/

# Inicia o download das imagens e depois processa
python download_amazon.py && python main.py

# Desativar o ambiente Conda quando você terminar, se desejar
conda deactivate