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

oldbands = abrir_old_json()

dirs = get_dirs()

# Acessando os diretórios usando as chaves do dicionário
dir_in = dirs['dir_in']
dir_main = dirs['dir_main']
dir_out = dirs['dir_out']
dir_libs = dirs['dir_libs']
dir_shapefiles = dirs['dir_shapefiles']
dir_colortables = dirs['dir_colortables']
dir_logos = dirs['dir_logos']
dir_temp = dirs['dir_temp']
arq_log = dirs['arq_log']

