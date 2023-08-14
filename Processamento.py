#!/opt/miniconda3/envs/goes/bin/python3
# -*- coding: utf-8 -*-
# Packages = conda create --name goes -c conda-forge matplotlib netcdf4 cartopy boto3 gdal scipy pandas scp
# Packages = apt install ffmpeg
# ====================================================================================================== #
# Manipulando imagens GOES-16 NetCDF's
# ===================================# Bibliotecas necessarias ========================================= #
import datetime
import logging  # Utilitario para criar os logs
import time     
from modulos.check_new_images import check_images # Checa as se há novas imagens para processamento
from modulos.process import process_gif, processing # Processa imagens
from modulos.logs import conf_log, finalize_log_time # Cria os arquivos de logs
from modulos.remove_images import remove_images
from modulos.quantity_products import quantity_products
from modulos.send_products import send_products

# ===================================# Bibliotecas necessarias ========================================= #


# ============================================# Variaveis ============================================== #
dir_in = "/home/guimoura/download_amazon/goes/"
dir_main = "/home/guimoura/download_amazon/"
dir_out = dir_main + "output/"
dir_libs = dir_main + "libs/"
dir_shapefiles = dir_main + "shapefiles/"
dir_colortables = dir_main + "colortables/"
dir_logos = dir_main + "logos/"
dir_temp = dir_main + "temp/"
arq_log = "/home/guimoura/download_amazon/logs/Processamento-GOES_" + str(datetime.date.today()) + ".log"
# ============================================# Variaveis ============================================== #


# ============================================# Bands Dicionario ============================================== #

# Dicionarios das bandas key : value
bands = {}

# Todas as bandas da 01 a 21 recebem False      bands = {"01": False, "02": False......
for num in range(1, 22):
    b = str(num).zfill(2)
    bands[f'{b}'] = False
br = True
sp = True

# ============================================# Bands Dicionario ============================================== #

# ============================================# Main ============================================== #

# configura o log
conf_log(arq_log)
# Log start time
start = time.time()
# Retorna o dicionario de bandas recebendo key : value (True/False) para saber quais imagens são novas
bands = check_images(bands, dir_in, dir_temp)

# Realiza etapas de processamento se houver alguma nova imagem para todas as bandas
if any(bands[key] for key in bands): # Se qualquer banda 'key': True   execute o processamento da banda ->
    # Realiza processamento da imagens
    bands = processing(bands, br, sp)
    # Remove imagens processadas no caso o .nc
    remove_images(bands, br, sp, dir_temp, dir_in)
    # Controle de quantidade de produtos devemos manter para produção do gif
    quantity_products(dir_out)
    # Realiza processamento do gif
    process_gif(bands, br, sp)
    # Envia as imagens para o site cpa.unicamp
    # send_products(s_br, s_sp, dir_our)
else:
    logging.info("")
    logging.info("SEM NOVAS IMAGENS PARA PROCESSAMENTO")


# Finaliza o script
finalize_log_time(start)

# ============================================# Main ============================================== #