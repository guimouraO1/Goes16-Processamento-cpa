#!/opt/miniconda3/envs/goes/bin/python3
# -*- coding: utf-8 -*-
# Packages = conda create --name goes -c conda-forge matplotlib netcdf4 cartopy boto3 gdal scipy pandas scp
# Packages = apt install ffmpeg
# ======================================================================================================
# Manipulando imagens GOES-16 NetCDF's
# ===================================# Bibliotecas necessarias =========================================
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # Plotagem de resultados, textos, logos, etc.
from matplotlib import cm  # Utilitario para paletas de cores
import cartopy  # Inserir mapas, shapefiles, paralelos, meridianos, latitudes, longitudes, etc.
import cartopy.crs as ccrs  # Utilitario para sistemas de referência e coordenadas
import cartopy.io.shapereader as shpreader  # Utilitario para leitura de shapefiles
from osgeo import gdal  # Utilitario para a biblioteca GDAL
from osgeo import osr  # Utilitario para a biblioteca GDAL
from netCDF4 import Dataset  # Utilitario para a biblioteca NetCDF4
import numpy as np  # Suporte para arrays e matrizes multidimensionais, com diversas funções matemáticas para trabalhar com estas estruturas
from utilities import load_cpt  # Funcao para ler as paletas de cores de arquivos CPT
from utilities import download_prod  # Funcao para download dos produtos do goes disponiveis na amazon
import datetime  # Utilitario para datas e horas
import time  # Utilitario para trabalhar com tempos
import os  # Utilitario para trabalhar com chamadas de sistema
import re  # Utilitario para trabalhar com expressoes regulares
import logging  # Utilitario para criar os logs
from difflib import Differ  # Utilitario para verificar as diferencas entre dois arquivos
from string import ascii_letters, digits  # Utilitario para trabalhar com ascii
from multiprocessing import Process  # Utilitario para multiprocessamento
import paramiko  # Utilitario para gerencia conexao SSH
import scp  # Utilitario para envio de arquivos com SCP
from shutil import copyfile  # Utilitario para copia de arquivos
from shapely.geometry import Point

gdal.PushErrorHandler('CPLQuietErrorHandler')   # Ignore GDAL warnings


# ======================================================================================================


# =============================================# Funcoes ===============================================
def read_process_file(banda):
    # Le o arquivo de processamento e retorna a lista
    with open(f'{dir_temp}{banda}_process.txt', 'r') as fo:
        return fo.readlines()


def check_images(c_bands):
    global dir_in, dir_temp
    logging.info("VERIFICANDO NOVAS IMAGENS")

    def alphanumeric_key(text):
        """Return a key based on letters and digits in `text`."""
        return [c.lower() for c in text if c in ascii_letters + digits]

    def write_new_file(banda, file):
        # Ordena de forma alfabetica a lista
        file.sort(key=alphanumeric_key)
        # Cria o arquivo com as novas imagens que estao na lista
        with open(f'{dir_temp}{banda}_new.txt', 'w') as fo:
            fo.writelines(map(lambda f: f + '\n', file))

    def write_process_file(banda):
        # Cria o arquivo band??_old.txt se nao existe
        if not os.path.isfile(f'{dir_temp}{banda}_old.txt'):
            with open(f'{dir_temp}{banda}_old.txt', 'w') as fo:
                fo.close()

        # Le os arquivos band??_old.txt e band??_new.txt
        with open(f'{dir_temp}{banda}_old.txt', 'r') as old, open(f'{dir_temp}{banda}_new.txt', 'r') as new:
            differ = Differ()
            # Realiza a comparacao entre os arquivos e cria uma lista de imagens que estao unicamente no arquivo band??_new.txt
            process_list = [line.strip()[2::] for line in differ.compare(old.readlines(), new.readlines()) if line.startswith('+')]

        # Cria o arquivo band??_process.txt com as imagens para processamento
        with open(f'{dir_temp}{banda}_process.txt', 'w') as process:
            # Escreve as imagens da lista no arquivo, sendo cada uma em uma linha
            process.writelines(map(lambda f: f + '\n', process_list))

        # Se houveram imagens para processamento, retorna True, caso contrario retorna False
        if process_list:
            return True
        else:
            return False

    # Contado para checagem de novas imagens nas 16 bandas
    for x in range(1, 2):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        # Cria uma lista com os itens presentes no diretorio da banda que sao arquivo e terminam com ".nc"
        imagens = [f for f in os.listdir(f'{dir_in}band{b}') if os.path.isfile(os.path.join(f'{dir_in}band{b}', f)) and re.match('^CG_ABI-L2-CMIPF-M[0-9]C[0-1][0-9]_G16_s.+_e.+_c.+.nc$', f)]
        # Se houver arquivos na lista, realiza o organizacao dos arquivos de controle e processamento, caso contrario aponta False no dicionario de controle das bandas
        if imagens:
            # Cria o arquivo com a nova lista de imagens
            write_new_file(f'band{b}', imagens)
            # Realiza a comparacao da lista nova com a lista antiga, se houver novas imagens, cria o arquivo de processamento e aponta True no dicionario de controle das bandas
            c_bands[b] = write_process_file(f'band{b}')
            if c_bands[b]:
                logging.info(f'Novas imagens Banda {b}')
            else:
                logging.info(f'Sem novas imagens Banda {b}')
                os.remove(f'{dir_temp}band{b}_process.txt')
            # Transforma o arquivo band??_new.txt para band??_old.txt
            os.replace(f'{dir_temp}band{b}_new.txt', f'{dir_temp}band{b}_old.txt')
        else:
            c_bands[b] = False
            logging.info(f'Sem novas imagens Banda {b}')

    # Retorna o dicionario de controle das bandas
    return c_bands

def reproject(reproj_file, reproj_var, reproj_extent, reproj_resolution):
    global dir_in

    def get_geot(ex, nlines, ncols):
        # Compute resolution based on data dimension
        resx = (ex[2] - ex[0]) / ncols
        resy = (ex[3] - ex[1]) / nlines
        return [ex[0], resx, 0, ex[3], 0, -resy]

    if reproj_extent == 'br':
        # Brasil
        r_extent = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat
    elif reproj_extent == 'sp':
        # São Paulo
        r_extent = [-53.25, -26.0, -44.0, -19.5]  # Min lon, Min lat, Max lon, Max lat
    else:
        r_extent = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat

    # GOES-16 Spatial Reference System
    source_prj = osr.SpatialReference()
    source_prj.ImportFromProj4('+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 +lat_0=0.0 +lon_0=-75 +sweep=x +no_defs')
    # Lat/lon WSG84 Spatial Reference System
    target_prj = osr.SpatialReference()
    target_prj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # Abrindo imagem com a biblioteca GDAL
    raw = gdal.Open(f'NETCDF:{reproj_file}:' + reproj_var, gdal.GA_ReadOnly)

    # Lendo os metadados do cabecalho
    metadata = raw.GetMetadata()
    scale = float(metadata.get(reproj_var + '#scale_factor'))
    offset = float(metadata.get(reproj_var + '#add_offset'))
    undef = float(metadata.get(reproj_var + '#_FillValue'))
    file_dtime = metadata.get('NC_GLOBAL#time_coverage_start')
    file_satellite = metadata.get('NC_GLOBAL#platform_ID')[1:3]

    # Setup projection and geo-transformation
    raw.SetProjection(source_prj.ExportToWkt())
    # raw.SetGeoTransform(raw.GetGeoTransform())
    GOES16_EXTENT = [-5434894.885056, -5434894.885056, 5434894.885056, 5434894.885056]
    raw.SetGeoTransform(get_geot(GOES16_EXTENT, raw.RasterYSize, raw.RasterXSize))

    # Compute grid dimension
    KM_PER_DEGREE = 111.32
    sizex = int(((r_extent[2] - r_extent[0]) * KM_PER_DEGREE) / reproj_resolution)
    sizey = int(((r_extent[3] - r_extent[1]) * KM_PER_DEGREE) / reproj_resolution)

    # Get memory driver
    driver = gdal.GetDriverByName('MEM')

    # Create grid
    grid = driver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)

    # Setup projection and geo-transformation
    grid.SetProjection(target_prj.ExportToWkt())
    grid.SetGeoTransform(get_geot(r_extent, grid.RasterYSize, grid.RasterXSize))

    # Perform the projection/resampling
    gdal.ReprojectImage(raw, grid, source_prj.ExportToWkt(), target_prj.ExportToWkt(), gdal.GRA_NearestNeighbour, options=['NUM_THREADS=ALL_CPUS'])

    # Close file
    raw = None
    del raw

    # Read grid data
    array = grid.ReadAsArray()

    # Mask fill values (i.e. invalid values)
    np.ma.masked_where(array, array == -1, False)

    # Aplicando scale, offset
    array = array * scale + offset

    grid.GetRasterBand(1).SetNoDataValue(-1)
    grid.GetRasterBand(1).WriteArray(array)

    # Define the parameters of the output file
    kwargs = {'format': 'netCDF',
              'dstSRS': target_prj,
              'outputBounds': (r_extent[0], r_extent[3], r_extent[2], r_extent[1]),
              'outputBoundsSRS': target_prj,
              'outputType': gdal.GDT_Float32,
              'srcNodata': undef,
              'dstNodata': 'nan',
              'resampleAlg': gdal.GRA_NearestNeighbour}

    reproj_file = reproj_file.split('/')
    reproj_file.reverse()
    r_file = reproj_file[0].replace('.nc', f'_reproj_{reproj_extent}.nc')
    gdal.Warp(f'{dir_in}{reproj_file[1]}/{r_file}', grid, **kwargs)

    return file_dtime, file_satellite, grid

def process_band_cmi(file, ch, v_extent):
    global dir_shapefiles, dir_colortables, dir_logos, dir_out
    file_var = 'CMI'
    # Captura a hora para contagem do tempo de processamento da imagem
    processing_start_time = time.time()
    # Area de interesse para recorte
    if v_extent == 'br':
        # Brasil
        extent = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat
        # Choose the image resolution (the higher the number the faster the processing is)
        resolution = 4.0
    elif v_extent == 'sp':
        # São Paulo
        extent = [-53.25, -26.0, -44.0, -19.5]  # Min lon, Min lat, Max lon, Max lat
        # Choose the image resolution (the higher the number the faster the processing is)
        resolution = 1.0
    else:
        extent = [-115.98, -55.98, -25.01, 34.98]  # Min lon, Min lat, Max lon, Max lat
        resolution = 2.0

    # Reprojetando imagem CMI e recebendo data/hora da imagem, satelite e caminho absoluto do arquivo reprojetado
    dtime, satellite, reproject_band = reproject(file, file_var, v_extent, resolution)

    if 1 <= int(ch) <= 6:
        data = reproject_band.ReadAsArray()
    else:
        data = reproject_band.ReadAsArray() - 273.15
    reproject_band = None
    del reproject_band

    # Comprimentos de ondas do canais ABI
    wavelenghts = ['[]', '[0.47 μm]', '[0.64 μm]', '[0.865 μm]', '[1.378 μm]', '[1.61 μm]', '[2.25 μm]', '[3.90 μm]', '[6.19 μm]',
                   '[6.95 μm]', '[7.34 μm]', '[8.50 μm]', '[9.61 μm]', '[10.35 μm]', '[11.20 μm]', '[12.30 μm]', '[13.30 μm]']

    # Formatando data para plotar na imagem e salvar o arquivo
    date = (datetime.datetime.strptime(dtime, '%Y-%m-%dT%H:%M:%S.%fZ'))
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')
    date_file = date.strftime('%Y%m%d_%H%M%S')

    # Definindo unidade de medida da imagem de acordo com o canal
    if 1 <= int(ch) <= 6:
        unit = "Albedo (%)"
    else:
        unit = "Brightness Temperature [°C]"
    # Formatando a descricao a ser plotada na imagem
    description = f" GOES-{satellite} ABI CMI Band {ch} {wavelenghts[int(ch)]} {unit} {date_img}"
    institution = "CEPAGRI - UNICAMP"

    # Definindo tamanho da imagem de saida
    d_p_i = 150
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

    # Utilizando projecao geoestacionaria no cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())

    if v_extent == 'br':
        # Adicionando o shapefile dos estados brasileiros
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil/BR_UF_2020').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='cyan', facecolor='none', linewidth=0.7)
    elif v_extent == 'sp':
        # Adicionando o shapefile dos estados brasileiros e cidade de campinas
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil/BR_UF_2020').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='cyan', facecolor='none', linewidth=0.7)
        shapefile = list(shpreader.Reader(dir_shapefiles + 'campinas/campinas').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='yellow', facecolor='none', linewidth=1)

    # Adicionando  linhas dos litorais
    ax.coastlines(resolution='10m', color='cyan', linewidth=0.5)
    # Adicionando  linhas das fronteiras
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='cyan', linewidth=0.5)
    # Adicionando  paralelos e meridianos
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=0.7, linestyle='--', linewidth=0.2, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5))
    gl.top_labels = False
    gl.right_labels = False

    # Definindo a paleta de cores da imagem de acordo com o canal
    # Convertendo um arquivo CPT para ser usado em Python atraves da funcao loadCPT
    # CPT archive: http://soliton.vm.bytemark.co.uk/pub/cpt-city/
    if 1 <= int(ch) <= 6:
        cpt = load_cpt(dir_colortables + 'Square Root Visible Enhancement.cpt')
    elif int(ch) == 7:
        cpt = load_cpt(dir_colortables + 'SVGAIR2_TEMP.cpt')
    elif 8 <= int(ch) <= 10:
        cpt = load_cpt(dir_colortables + 'SVGAWVX_TEMP.cpt')
    else:
        cpt = load_cpt(dir_colortables + 'IR4AVHRR6.cpt')
    my_cmap = cm.colors.LinearSegmentedColormap('cpt', cpt)  # Criando uma paleta de cores personalizada

    # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
    img_extent = [extent[0], extent[2], extent[1], extent[3]]  # Min lon, Max lon, Min lat, Max lat

    # Plotando a imagem
    # Definindo os valores maximos e minimos da imagem de acordo com o canal e paleta de cores utilizada
    if 1 <= int(ch) <= 6:
        img = ax.imshow(data, origin='upper', vmin=0, vmax=1, cmap=my_cmap, extent=img_extent)
    elif 7 <= int(ch) <= 10:
        img = ax.imshow(data, origin='upper', vmin=-112.15, vmax=56.85, cmap=my_cmap, extent=img_extent)
    else:
        img = ax.imshow(data, origin='upper', vmin=-103, vmax=84, cmap=my_cmap, extent=img_extent)

    # Adicionando barra da paleta de cores de acordo com o canal
    # Criando novos eixos de acordo com a posicao da imagem
    cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
    if 1 <= int(ch) <= 6:
        cb = plt.colorbar(img, orientation="horizontal", cax=cax0, ticks=[0.2, 0.4, 0.6, 0.8])
        cb.ax.set_xticklabels(['0.2', '0.4', '0.6', '0.8'])
        cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
    elif 7 <= int(ch) <= 10:
        cb = plt.colorbar(img, orientation="horizontal", cax=cax0)
        cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
    else:
        cb = plt.colorbar(img, orientation="horizontal", cax=cax0)
        cb.ax.tick_params(axis='x', colors='gray', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
    cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
    cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
    cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de cores

    # Adicionando descricao da imagem
    # Criando novos eixos de acordo com a posicao da imagem
    cax1 = fig.add_axes([ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.026, ax.get_position().width - 0.003, 0.0125])
    cax1.patch.set_color('black')  # Alterando a cor do novo eixo
    cax1.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
    cax1.text(0.85, 0.13, institution, color='yellow', size=10)  # Adicionando texto
    cax1.xaxis.set_visible(False)  # Removendo rotulos do eixo X
    cax1.yaxis.set_visible(False)  # Removendo rotulos do eixo Y

    # Adicionando os logos
    logo_noaa = plt.imread(dir_logos + 'NOAA_Logo.png')  # Lendo o arquivo do logo
    logo_goes = plt.imread(dir_logos + 'GOES_Logo.png')  # Lendo o arquivo do logo
    logo_cepagri = plt.imread(dir_logos + 'CEPAGRI-Logo.png')  # Lendo o arquivo do logo
    fig.figimage(logo_noaa, 32, 233, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_goes, 10, 150, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_cepagri, 10, 70, zorder=3, alpha=0.8, origin='upper')  # Plotando logo

    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}band{ch}/band{ch}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
    # Fecha a janela para limpar a memoria
    plt.close()
    # Realiza o log do calculo do tempo de processamento da imagem
    logging.info(f'{file} - {v_extent} - {str(round(time.time() - processing_start_time, 4))} segundos')


def processing(p_bands, p_br, p_sp):
    global dir_in, dir_temp
    # Cria lista vazia para controle do processamento paralelo
    process_br = []
    # Cria lista vazia para controle processamento paralelo
    process_sp = []
    # Variavel de processamento diario do NDVI
    ndvi_diario = False
    ultima_diario = []

    # Processando arquivos das bandas do ABI
    # Se a variavel de controle de processamento do brasil for True, realiza o processamento
    if p_br:
        logging.info("")
        logging.info('PROCESSANDO IMAGENS "BR"...')
        # Contador para processamento nas 16 bandas
        for x in range(1, 2):
            # Transforma o inteiro contador em string e com 2 digitos
            b = str(x).zfill(2)
            # Se o dicionario de controle das bandas apontar True para essa banda, realiza o processamento
            if p_bands[b]:
                # Le o arquivo de processamento da banda
                img = read_process_file(f'band{b}')
                # Para cada imagem no arquivo, cria um processo chamando a funcao de processamento
                for i in img:
                    # Remove possiveis espacos vazios no inicio ou final da string
                    i = i.strip()
                    # Tenta realizar o processamento da imagem
                    try:
                        # Cria o processo com a funcao de processamento
                        process = Process(target=process_band_cmi, args=(f'{dir_in}band{b}/{i}', b, "br"))
                        # Adiciona o processo na lista de controle do processamento paralelo
                        process_br.append(process)
                        # Inicia o processo
                        process.start()
                    # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
                    except OSError as ose:
                        # Realiza o log do erro
                        logging.info(f'Erro Arquivo - OSError - {i}')
                        logging.info(str(ose))
                        # Remove a imagem com erro de processamento
                        os.remove(f'{dir_in}band{b}/{i}')
                    except AttributeError as ae:
                        # Realiza o log do erro
                        logging.info(f'Erro Arquivo - AttributeError - {i}')
                        logging.info(str(ae))
                        # Remove a imagem com erro de processamento
                        os.remove(f'{dir_in}band{b}/{i}')
            # Se o dicionario de controle das bandas apontar False para essa banda, nao realiza o processamento e continua a tentativa nas outras bandas
            else:
                continue
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in process_br:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()
        # Limpa lista vazia para controle do processamento paralelo
        process_br = []
    # Se a variavel de controle de processamento do estado de sao paulo for True, realiza o processamento
    if p_sp:
        logging.info("")
        logging.info('PROCESSANDO IMAGENS "SP"...')
        # Contador para processamento nas 16 bandas
        for x in range(1, 17):
            # Transforma o inteiro contador em string e com 2 digitos
            b = str(x).zfill(2)
            # Se o dicionario de controle das bandas apontar True para essa banda, realiza o processamento
            if p_bands[b]:
                # Le o arquivo de processamento da banda
                img = read_process_file(f'band{b}')
                # Para cada imagem no arquivo, cria um processo chamando a funcao de processamento
                for i in img:
                    # Remove possiveis espacos vazios no inicio ou final da string
                    i = i.strip()
                    # Tenta realizar o processamento da imagem
                    try:
                        # Cria o processo com a funcao de processamento
                        process = Process(target=process_band_cmi, args=(f'{dir_in}band{b}/{i}', b, "sp"))
                        # Adiciona o processo na lista de controle do processamento paralelo
                        process_sp.append(process)
                        # Inicia o processo
                        process.start()
                    # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
                    except OSError as ose:
                        # Realiza o log do erro
                        logging.info("Erro Arquivo - OSError")
                        logging.info(str(ose))
                        # Remove a imagem com erro de processamento
                        os.remove(f'{dir_in}band{b}/{i}')
                    except AttributeError as ae:
                        # Realiza o log do erro
                        logging.info(f'Erro Arquivo - AttributeError - {i}')
                        logging.info(str(ae))
                        # Remove a imagem com erro de processamento
                        os.remove(f'{dir_in}band{b}/{i}')
            # Se o dicionario de controle das bandas apontar False para essa banda, nao realiza o processamento e continua a tentativa nas outras bandas
            else:
                continue
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in process_sp:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()
        # Limpa lista vazia para controle do processamento paralelo
        process_sp = []

def remove_images(r_bands, r_br, r_sp):
    global dir_temp, dir_in
    logging.info("")
    logging.info('REMOVENDO IMAGENS PROCESSADAS')

    # Contador para remover imagens nas 16 bandas
    for x in range(1, 2):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        # Se o dicionario de controle das bandas apontar True para essa banda, remove a imagem
        if r_bands[b]:
            # Le o arquivo de processamento da banda
            img = read_process_file(f'band{b}')
            logging.info(f'Removendo imagens banda {b}')
            # Para cada imagem no arquivo, realiza a remocao
            for i in img:
                # Remove possiveis espacos vazios no inicio ou final da string
                i = i.strip()
                try:
                    # Remove a imagem
                    os.remove(f'{dir_in}band{b}/{i}')
                except FileNotFoundError as fnfe:
                    # Realiza o log do erro
                    logging.info(f'Erro Arquivo - FileNotFoundError - {i}')
                    logging.info(str(fnfe))
                try:
                    if r_br:
                        os.remove(f'{dir_in}band{b}/{i.replace(".nc", "_reproj_br.nc")}')
                except FileNotFoundError as fnfe:
                    # Realiza o log do erro
                    logging.info(f'Erro Arquivo - FileNotFoundError - {i.replace(".nc", "_reproj_br.nc")}')
                    logging.info(str(fnfe))
                try:
                    if r_sp:
                        os.remove(f'{dir_in}band{b}/{i.replace(".nc", "_reproj_sp.nc")}')
                except FileNotFoundError as fnfe:
                    # Realiza o log do erro
                    logging.info(f'Erro Arquivo - FileNotFoundError - {i.replace(".nc", "_reproj_sp.nc")}')
                    logging.info(str(fnfe))
            logging.info(f'Removendo arquivo de processamento da banda {b}')
            # Remove o arquivo de processamento
            os.remove(f'{dir_temp}band{b}_process.txt')

def quantity_products():
    global dir_out
    logging.info("")
    logging.info("VERIFICANDO NUMERO DE PRODUTOS...")
    aux = False

    # Contador para verificar produtos nas 16 bandas
    for x in range(1, 17):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        # Verifica a quantidade de itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_br.png$"
        result_br = len([name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_br.png$', name)])
        # Verifica a quantidade de itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_sp.png$"
        result_sp = len([name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_sp.png$', name)])
        # Subtrai o limite de imagens BR
        result_br -= 48
        # Subtrai o limite de imagens SP
        result_sp -= 48
        # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
        if result_br > 0:
            aux = True
            # Cria uma lista com os itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_br.png$"
            prod_br = [name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            prod_br.sort()
            # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
            for y in range(0, result_br):
                # Remove arquivo em excesso
                os.remove(f'{dir_out}band{b}/{prod_br[y]}')
        # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
        if result_sp > 0:
            aux = True
            # Cria uma lista com os itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_sp.png$"
            prod_sp = [name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_sp.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            prod_sp.sort()
            # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
            for y in range(0, result_sp):
                # Remove arquivo em excesso
                os.remove(f'{dir_out}band{b}/{prod_sp[y]}')

    # Verificar produtos truecolor
    # Verifica a quantidade de itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_br.png$"
    result_br = len([name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_br.png$', name)])
    # Verifica a quantidade de itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_sp.png$"
    result_sp = len([name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_sp.png$', name)])
    # Subtrai o limite de imagens BR
    result_br -= 48
    # Subtrai o limite de imagens SP
    result_sp -= 48
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_br > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_br.png$"
        prod_br = [name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_br.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_br.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_br):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}truecolor/{prod_br[y]}')
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_sp > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_sp.png$"
        prod_sp = [name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_sp.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_sp.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_sp):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}truecolor/{prod_sp[y]}')

    # Verificar produtos rrqpef
    # Verifica a quantidade de itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_br.png$"
    result_br = len([name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_br.png$', name)])
    # Verifica a quantidade de itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_sp.png$"
    result_sp = len([name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_sp.png$', name)])
    # Subtrai o limite de imagens BR
    result_br -= 48
    # Subtrai o limite de imagens SP
    result_sp -= 48
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_br > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_br.png$"
        prod_br = [name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_br.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_br.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_br):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}rrqpef/{prod_br[y]}')
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_sp > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_sp.png$"
        prod_sp = [name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_sp.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_sp.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_sp):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}rrqpef/{prod_sp[y]}')

    # Verificar produtos glm
    # Verifica a quantidade de itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_br.png$"
    result_br = len([name for name in os.listdir(f'{dir_out}glm') if os.path.isfile(os.path.join(f'{dir_out}glm', name)) and re.match('^glm_.+_.+_br.png$', name)])
    # Verifica a quantidade de itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_sp.png$"
    result_sp = len([name for name in os.listdir(f'{dir_out}glm') if os.path.isfile(os.path.join(f'{dir_out}glm', name)) and re.match('^glm_.+_.+_sp.png$', name)])
    # Subtrai o limite de imagens BR
    result_br -= 48
    # Subtrai o limite de imagens SP
    result_sp -= 48
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_br > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_br.png$"
        prod_br = [name for name in os.listdir(f'{dir_out}glm') if os.path.isfile(os.path.join(f'{dir_out}glm', name)) and re.match('^glm_.+_.+_br.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_br.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_br):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}glm/{prod_br[y]}')
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_sp > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_sp.png$"
        prod_sp = [name for name in os.listdir(f'{dir_out}glm') if os.path.isfile(os.path.join(f'{dir_out}glm', name)) and re.match('^glm_.+_.+_sp.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_sp.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_sp):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}glm/{prod_sp[y]}')

    # Verificar produtos ndvi
    # Verifica a quantidade de itens no diretorio dos produtos ndvi que sao arquivos e se encaixa na expressao regular "^ndvi_.+_.+_br.png$"
    result_br = len([name for name in os.listdir(f'{dir_out}ndvi') if os.path.isfile(os.path.join(f'{dir_out}ndvi', name)) and re.match('^ndvi_.+_.+_br.png$', name)])
    # Subtrai o limite de imagens BR
    result_br -= 48
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_br > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos ndvi que sao arquivos e se encaixa na expressao regular "^ndvi_.+_.+_br.png$"
        prod_br = [name for name in os.listdir(f'{dir_out}ndvi') if os.path.isfile(os.path.join(f'{dir_out}ndvi', name)) and re.match('^ndvi_.+_.+_br.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_br.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_br):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}ndvi/{prod_br[y]}')

    # Verificar produtos fdcf
    # Verifica a quantidade de itens no diretorio dos produtos fdcf que sao arquivos e se encaixa na expressao regular "^fdcf_.+_.+_br.png$"
    result_br = len([name for name in os.listdir(f'{dir_out}fdcf') if os.path.isfile(os.path.join(f'{dir_out}fdcf', name)) and re.match('^fdcf_.+_.+_br.png$', name)])
    # Subtrai o limite de imagens BR
    result_br -= 48
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_br > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos fdcf que sao arquivos e se encaixa na expressao regular "^fdcf_.+_.+_br.png$"
        prod_br = [name for name in os.listdir(f'{dir_out}fdcf') if os.path.isfile(os.path.join(f'{dir_out}fdcf', name)) and re.match('^fdcf_.+_.+_br.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_br.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_br):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}fdcf/{prod_br[y]}')

    if aux:
        logging.info("Produtos em excesso removidos com sucesso")
    else:
        logging.info("Quantidade de produtos dentro do limite")


def process_gif(g_bands, g_br, g_sp):
    global dir_out
    # Cria lista vazia para controle do processamento paralelo
    gif_br = []
    # Cria lista vazia para controle do processamento paralelo
    gif_sp = []

    # Chamada de sistema para o software ffmpeg realizar a criacao do gif animado
    def create_gif_cmi(banda, roi):
        os.system(f'/usr/bin/ffmpeg -y -v warning -framerate 4 -pattern_type glob -i "{dir_out}band{banda}/band{banda}_*_*_{roi}.png" "{dir_out}band{banda}/band{banda}_{roi}.gif"')

    # Chamada de sistema para o software ffmpeg realizar a criacao do gif animado
    def create_gif(file_type, roi):
        os.system(f'/usr/bin/ffmpeg -y -v warning -framerate 4 -pattern_type glob -i "{dir_out}{file_type}/{file_type}_*_*_{roi}.png" "{dir_out}{file_type}/{file_type}_{roi}.gif"')

    # Se a variavel de controle de processamento do brasil for True, cria o gif
    if g_br:
        logging.info('')
        logging.info('CRIANDO GIF ANIMADO "BR"...')
        # Contador para gif nas 16 bandas
        for x in range(1, 2):
            # Transforma o inteiro contador em string e com 2 digitos
            b = str(x).zfill(2)
            # Se o dicionario de controle das bandas apontar True para essa banda, cria o gif
            if g_bands[b]:
                logging.info('Gif BR banda ' + b)
                # Cria o processo com a funcao gif
                process = Process(target=create_gif_cmi, args=(b, "br"))
                # Adiciona o processo na lista de controle do processamento paralelo
                gif_br.append(process)
                # Inicia o processo
                process.start()
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in gif_br:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()

    # Se a variavel de controle de processamento do estado de sao paulo for True, cria o gif
    if g_sp:
        logging.info('')
        logging.info('CRIANDO GIF ANIMADO "SP"...')
        # Contador para gif nas 16 bandas
        for x in range(1, 2):
            # Transforma o inteiro contador em string e com 2 digitos
            b = str(x).zfill(2)
            # Se o dicionario de controle das bandas apontar True para essa banda, cria o gif
            if g_bands[b]:
                logging.info('Gif SP banda ' + b)
                # Cria o processo com a funcao gif
                process = Process(target=create_gif_cmi, args=(b, "sp"))
                # Adiciona o processo na lista de controle do processamento paralelo
                gif_br.append(process)
                # Inicia o processo
                process.start()
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in gif_sp:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()

def finalize(s):
    # Capturando data/hora final
    fim = datetime.datetime.now().strftime("%d/%m/%Y-%H:%M:%S")
    # Realiza o log do calculo do tempo de execucao
    logging.info("")
    logging.info("Tempo gasto " + str(round(time.time() - s, 4)) + ' segundos')
    logging.info("")
    logging.info("=========================================================================================================")
    logging.info("=                         PROCESSAMENTO ENCERRADO GOES AS " + fim + "                           =")
    logging.info("=========================================================================================================")
    logging.info("")
    logging.info("")


# ======================================================================================================


# ============================================# Variaveis ==============================================
dir_in = "/home/guimoura/download_amazon/goes/"
dir_main = "/home/guimoura/download_amazon/"
dir_out = dir_main + "output/"
dir_libs = dir_main + "libs/"
dir_shapefiles = dir_main + "shapefiles/"
dir_colortables = dir_main + "colortables/"
dir_logos = dir_main + "logos/"
dir_temp = dir_main + "temp/"
arq_log = "/home/guimoura/download_amazon/logs/Processamento-GOES_" + str(datetime.date.today()) + ".log"

bands = {"01": False, "02": False, "03": False, "04": False, "05": False, "06": False, "07": False, "08": False,
         "09": False, "10": False, "11": False, "12": False, "13": False, "14": False, "15": False, "16": False,
         "17": False,  # Band 17 = TrueColor
         "18": False,  # Band 18 = RRQPEF
         "19": False,  # Band 19 = GLM
         "20": False,  # Band 20 = NDVI
         "21": False}  # Band 21 = FDCF
br = True
sp = True

# ======================================================================================================


# ===============================================# Main ================================================
# Configurando log
logging.basicConfig(filename=arq_log, level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S")
# Capturando data/hora inicio
inicio = datetime.datetime.now().strftime("%d/%m/%Y-%H:%M:%S")
logging.info("")
logging.info("")
logging.info("=========================================================================================================")
logging.info("=                             PROCESSANDO IMAGENS GOES AS " + inicio + "                           =")
logging.info("=========================================================================================================")
logging.info("")

# Captura a hora para contagem do tempo de execucao
start = time.time()
# Realizando a checagem de novas imagens para processamento
bands = check_images(bands)

# Realiza etapas de processamento se houver alguma nova imagem
if bands["01"] or bands["02"] or bands["03"] or bands["04"] or bands["05"] or bands["06"] or bands["07"] or bands["08"] or \
        bands["09"] or bands["10"] or bands["11"] or bands["12"] or bands["13"] or bands["14"] or bands["15"] or bands["16"]:
    # Realiza processamento da imagens
    bands = processing(bands, br, sp)
    # Remove imagens processadas
    remove_images(bands, br, sp)
    # Controle de quantidade de produtos devemos manter para produção do gif
    quantity_products()
    # Realiza processamento do gif
    process_gif(bands, br, sp)
    # Envia os produtos para o site
    #send_products(br, sp)
    # Finaliza o script
    finalize(start)
else:
    logging.info("")
    logging.info("SEM NOVAS IMAGENS PARA PROCESSAMENTO")
    # Finaliza o script
    finalize(start)
# ======================================================================================================
