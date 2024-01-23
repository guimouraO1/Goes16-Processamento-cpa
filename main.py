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


if __name__ == "__main__":
    
    # Importando os diretórios do modulo dirs.py
    dirs = get_dirs()
    dir_in = dirs['dir_in']
    dir_out = dirs['dir_out']
    arq_log = dirs['arq_log']
    dir_main = dirs['dir_main']

    # Dicionarios das bandas json   key : value
    bands = {}
    # Todas as bandas da 01 a 24 recebem False
    for num in range(1, 25):
        b = str(num).zfill(2)
        bands[f'{b}'] = False
        
    # Br e Sp representam se terá imagens para SP e BR
    br = True
    sp = True

    # configura o log
    conf_log(arq_log)

    # Log start time
    start = time.time()

    try:
        # Chama a função checarImagens para verificar a existência de novas imagens.
        bands = checar_imagens(bands, dir_in, dir_main)
        # Se tiver novas imagens para processamento:
        if any(bands[key] for key in bands):
            # Processa as imagens
            processamento_das_imagens(bands, br, sp, dir_in, dir_main)
            # Remove os arquivos .nc que já foram processados
            remover_imagens(bands, dir_in)
            # A função é chamada para controlar a quantidade de produtos (imagens) a serem mantidos para a produção de um GIF animado.
            quantity_products(dir_out)
            # Processa o gif
            process_gif(bands, br, sp, dir_out)
            # Envia produtos para o site
            # send_products(br, sp, dir_out)
        else:
            logging.info('Sem arquivos para processamento. \n')
    except Exception as error:
        logging.info('ERROR : ' + str(error))

    finalize_log_time(start)
