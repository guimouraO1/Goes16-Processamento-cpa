#!/opt/miniconda3/envs/goes/bin/python3
# -*- coding: utf-8 -*-
from osgeo import gdal
from modules.logs import conf_log_D, finalize_log_time_D 
from modules.utilities import download_dmw, download_glm, download_prod, download_cmi_joao
from modules.dirs import get_dirs
from datetime import timedelta, datetime
import logging
import os
import time

# Ignore GDAL warnings 
gdal.PushErrorHandler('CPLQuietErrorHandler')

if __name__ == "__main__":
    
    # Acessando os diretórios usando as chaves do dicionário
    dirs = get_dirs()
    dir_in = dirs['dir_in']
    arq_log = dirs['arq_log']

    # Inicia contador para o tempo gasto de download
    start = time.time()

    # Configura o log
    conf_log_D(arq_log)

    # Captura hora atual em UTC para download no site da Amazon
    data_hora_atual = datetime.utcnow()

    # Atrasa 10 min para entrar em conformidade com Amazon
    data_10_min = datetime.strftime(data_hora_atual-timedelta(minutes=10),'%Y%m%d%H%M')

    #Correção para poder fazer download em qualquer horário
    data_hora_download_file = data_10_min[0:11]+ '0'

    # Produtos de 10 em 10 minutos

   # Download das 16 bandas
    for x in range(1, 4):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        # Download band
        logging.info("")
        logging.info(f'Tentando download Band{b}...')
        try:
            download_cmi_joao(data_hora_download_file, x, dir_in + f'band{b}', logging)
        except Exception as e:
            print(f'{e}')
            continue
        
    for x in range(13, 14):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        # Download band
        logging.info("")
        logging.info(f'Tentando download Band{b}...')
        try:
            download_cmi_joao(data_hora_download_file, x, dir_in + f'band{b}', logging)
        except Exception as e:
            print(f'{e}')
            continue



    # download_glm(data_hora_download_file, dir_in + f'glm')

    # Produtos de 1 em 1 hora
    
    # Atrasa 1h para entrar em conformidade com Amazon
    data_1h = datetime.strftime(data_hora_atual-timedelta(hours=1),'%Y%m%d%H00')

    # Donwload do produto LST 2km Full disk
    #logging.info(f'Downloading file UTC lst_{data_1h}.nc')
    #download_prod(data_1h,'ABI-L2-LST2KMF',f'{dir_in}lst/')

    # Donwload do produto DMW
    logging.info(f'Downloading file UTC dmw_{data_1h}.nc')
    download_dmw(data_1h, 14, f'{dir_in}dmw')

    finalize_log_time_D(start)