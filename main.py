#!/opt/miniconda3/envs/goes/bin/python3
# -*- coding: utf-8 -*-
# Packages = conda create --name goes -c conda-forge matplotlib netcdf4 cartopy boto3 gdal scipy pandas scp
# Packages = apt install ffmpeg
# ================================================================================================== #
# Manipulando imagens GOES-16 NetCDF's
# ===================================# Bibliotecas necessarias #==================================== #
import time
import logging 
from modules.logs import conf_log, finalize_log_time # Cria os arquivos de logs
from modules.quantity_products import quantity_products
from modules.send_products import send_products
from modules.dirs import get_dirs
from modules.processamento import processamento_das_imagens
from modules.check_new_images import checar_imagens
from modules.remove import remover_imagens
from modules.process_gif import process_gif
# ===================================# Bibliotecas necessarias #==================================== #


#==================================#           Dicionário        #==================================#
dirs = get_dirs()
# Importando dirs do modulo dirs.py

dir_in = dirs['dir_in']
dir_out = dirs['dir_out']
arq_log = dirs['arq_log']
#==================================#           Dicionário        #==================================#


#==================================#     Dicionário das bandas   #==================================#
# Dicionarios das bandas key : value
bands = {}
# Todas as bandas da 01 a 21 recebem False      bands = {"01": False, "02": False......
for num in range(1, 22):
    b = str(num).zfill(2)
    bands[f'{b}'] = False
br = True
sp = True
#==================================#     Dicionário das bandas   #==================================#

# Script Antigo - 06/09/2023 11:39:06 Tempo gasto 423.6059 segundos |
# Script Atual  - 06/09/2023 11:31:15 Tempo gasto 291.632 segundos  | - Sem o envio de imagens para o site

# ========================================#     Main     #========================================== #
# configura o log
conf_log(arq_log)
# Log start time
start = time.time()

try:
    # Chama a função checarImagens para verificar a existência de novas imagens.
    bands = checar_imagens(bands, dir_in)
    # Se tiver novas imagens para processamento:
    if any(bands[key] for key in bands):
        # Processa as imagens
        processamento_das_imagens(bands, br, sp, dir_in)
        # Remove os arquivos .nc que já foram processados
        remover_imagens(bands, dir_in)
        # A função é chamada para controlar a quantidade de produtos (imagens) a serem mantidos para a produção de um GIF animado.
        quantity_products(dir_out)
        # Processa o gif
        process_gif(bands, br, sp, dir_out)
    else:
        logging.info('Sem arquivos para processamento. \n')
except Exception as error:
    logging.info('ERROR : ' + str(error))

finalize_log_time(start)
# ========================================#     Main     #========================================== #
