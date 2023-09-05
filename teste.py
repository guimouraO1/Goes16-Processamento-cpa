import os
import datetime
import re
from modules.processamento import abrir_old_json
import matplotlib
matplotlib.use('Agg')# Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt  # Plotagem de resultados, textos, logos, etc.
from matplotlib import cm  # Utilitario para paletas de cores
import cartopy  # Inserir mapas, shapefiles, paralelos, meridianos, latitudes, longitudes, etc.
import cartopy.crs as ccrs  # Utilitario para sistemas de referência e coordenadas
import cartopy.io.shapereader as shpreader  # Utilitario para leitura de shapefiles
from osgeo import gdal  # Utilitario para a biblioteca GDAL
from osgeo import osr  # Utilitario para a biblioteca GDAL
from netCDF4 import Dataset  # Utilitario para a biblioteca NetCDF4
import numpy as np  # Suporte para arrays e matrizes multidimensionais, com diversas funções matemáticas para trabalhar com estas estruturas
import datetime  # Utilitario para datas e horas
import time  # Utilitario para trabalhar com tempos
import os  # Utilitario para trabalhar com chamadas de sistema
import logging  # Utilitario para criar os logs
from multiprocessing import Process  # Utilitario para multiprocessamento
import re # Utilitario para trabalhar com expressoes regulares
from modules.dirs import get_dirs
from modules.utilities import load_cpt  # Funcao para ler as paletas de cores de arquivos CPT
from modules.utilities import download_prod
from modules.processamento import reproject
import cartopy.feature as cfeature # features
from modules.processamento import area_para_recorte
import json
from modules.check_new_images import modificar_chave_old_bands
import datetime
from shutil import copyfile  # Utilitario para copia de arquivos

oldbands = abrir_old_json()

dirs = get_dirs()

# Acessando os diretórios usando as chaves do dicionário
dir_in = dirs['dir_in']
dir_main = dirs['dir_main']
dir_out = dirs['dir_out']
dir_shapefiles = dirs['dir_shapefiles']
dir_colortables = dirs['dir_colortables']
dir_logos = dirs['dir_logos']
dir_temp = dirs['dir_temp']
arq_log = dirs['arq_log']

oldbands = abrir_old_json()

ch01 = oldbands['01']
ch02 = oldbands['02']
ch03 = oldbands['03']

# Carrega os arquivos de processamento das bandas para composicao do ndvi
file_ch02 = oldbands['02']

# Listar arquivos no diretório da pasta "band03"
file_ch03_dir = f'{dir_in}band03/'
file_ch03_list = os.listdir(file_ch03_dir)

date_now = datetime.datetime.now()
date_ini = datetime.datetime(date_now.year, date_now.month, date_now.day, int(13), int(00))
date_end = datetime.datetime(date_now.year, date_now.month, date_now.day, int(13), int(00)) + datetime.timedelta(hours=5, minutes=1)

date_file = datetime.datetime.strptime(file_ch02[file_ch02.find("M6C02_G16_s") + 11:file_ch02.find("_e") - 1], '%Y%j%H%M%S')

if date_ini <= date_file <= date_end:
   # Converte o nome do arquivo da banda 02 para o formato da banda 03
    file_ch03 = file_ch02[0:43].replace('M6C02', 'M6C03')

    # Verifica se há pelo menos um arquivo correspondente na banda 03
    found_match = any(file_ch03_candidate.startswith(file_ch03) for file_ch03_candidate in file_ch03_list)

    
print(found_match)
    