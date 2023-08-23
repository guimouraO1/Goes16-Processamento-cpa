#!/opt/miniconda3/envs/goes/bin/python3
# -*- coding: utf-8 -*-
# Packages = conda create --name goes -c conda-forge matplotlib netcdf4 cartopy boto3 gdal scipy pandas scp
# Packages = apt install ffmpeg
# ================================================================================================= #
# Manipulando imagens GOES-16 NetCDF's
# ===================================# Bibliotecas necessarias ==================================== #
import time    
from modules.logs import conf_log, finalize_log_time # Cria os arquivos de logs
from modules.quantity_products import quantity_products
from modules.send_products import send_products
from modules.dirs import get_dirs
from modules.process import processing
from modules.check_images import checarImagens
from modules.remove import removeImagens
from modules.process_gif import process_gif
# ===================================# Bibliotecas necessarias ==================================== #


# ============================================# Diretórios ========================================= #
dirs = get_dirs()
# Importando dirs do modulo dirs.py
arq_log = dirs['arq_log']
dir_in = dirs['dir_in']
dir_temp = dirs['dir_temp']
dir_out = dirs['dir_out']
# ============================================# Diretórios ========================================= #


# ============================================# Bands Dicionario ================================== #
# Dicionarios das bandas key : value
bands = {}
# Todas as bandas da 01 a 21 recebem False      bands = {"01": False, "02": False......
for num in range(1, 17):
    b = str(num).zfill(2)
    bands[f'{b}'] = False
br = True
sp = True
# ============================================# Bands Dicionario ================================== #


# ============================================# Main ============================================== #
# configura o log
conf_log(arq_log)
# Log start time
start = time.time()

try:
    print('\nChecando arquivos.nc...\n')
    # Chama a função checarImagens para verificar a existência de novas imagens.
    bands = checarImagens(bands, dir_in)

    if any(bands[key] for key in bands):
        print('Processando imagens....\n')
        processing(bands, br, sp, dir_in)
        print('Removendo arquivos.nc....\n')
        removeImagens(bands, dir_in)
        print('quantity prod....\n')
        quantity_products(dir_out)
        print('Processando gif...\n')
        process_gif(bands, br, sp, dir_out)
    else:
        print('Sem arquivos para processamento. \n')
except:
    print('Sem imagens novas \n')
    
    
finalize_log_time(start)
# ============================================# Main ============================================== #