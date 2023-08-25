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
import datetime  # Utilitario para datas e horas
import time  # Utilitario para trabalhar com tempos
import os  # Utilitario para trabalhar com chamadas de sistema
import logging  # Utilitario para criar os logs
from multiprocessing import Process  # Utilitario para multiprocessamento
from modules.dirs import get_dirs
import re # Utilitario para trabalhar com expressoes regulares
import json
from libs.utilities import load_cpt  # Funcao para ler as paletas de cores de arquivos CPT

gdal.PushErrorHandler('CPLQuietErrorHandler')   # Ignore GDAL warnings

# ============================================# Diretórios ========================================= #
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
# ============================================# Diretórios ========================================= #


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
    if reproj_var == 'BCM':  ### O arquivo Clear Sky não possui sacale e offset é um arquivo binário
        metadata = raw.GetMetadata()
        scale = 1
        offset = 0
        undef = float(metadata.get(reproj_var + '#_FillValue'))
        file_dtime = metadata.get('NC_GLOBAL#time_coverage_start')
        file_satellite = metadata.get('NC_GLOBAL#platform_ID')[1:3]

    else:  # Alteração João 26/08/2022
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

# Função para abrir o arquivo.json
def openOld():
    with open('oldBands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
        return oldImages

def process_band_rgb(rgb_type, v_extent, ch01=None, ch02=None, ch03=None):
    global dir_in, dir_shapefiles, dir_colortables, dir_logos, dir_out
    # Captura a hora para contagem do tempo de processamento da imagem
    processing_start_time = time.time()
    # Area de interesse para recorte
    if v_extent == 'br':
        # Brasil
        extent = [-90.0, -40.0, -20.0, 10.0]  # Min lon, Min lat, Max lon, Max lat
    elif v_extent == 'sp':
        # São Paulo
        extent = [-53.25, -26.0, -44.0, -19.5]  # Min lon, Min lat, Max lon, Max lat
    else:
        extent = [-115.98, -55.98, -25.01, 34.98]  # Min lon, Min lat, Max lon, Max lat

    if rgb_type == 'truecolor':
        # https://rammb.cira.colostate.edu/training/visit/quick_guides/
        # http://cimss.ssec.wisc.edu/goes/OCLOFactSheetPDFs/ABIQuickGuide_CIMSSRGB_v2.pdf
        # Lendo imagem CMI reprojetada
        
        reproject_ch01 = Dataset(ch01)
        reproject_ch02 = Dataset(ch02)
        reproject_ch03 = Dataset(ch03)
        
        # Coletando do nome da imagem o satelite e a data/hora
        satellite_ch01 = (ch01[ch01.find("M6C01_G") + 7:ch01.find("_s")])
        dtime_ch01 = ch01[ch01.find("M6C01_G16_s") + 11:ch01.find("_e") - 1]

        data_ch01 = reproject_ch01.variables['Band1'][:]
        data_ch02 = reproject_ch02.variables['Band1'][:]
        data_ch03 = reproject_ch03.variables['Band1'][:]
        reproject_ch01 = None
        del reproject_ch01
        reproject_ch02 = None
        del reproject_ch02
        reproject_ch03 = None
        del reproject_ch03

        # RGB Components
        R = data_ch02
        G = 0.45 * data_ch02 + 0.1 * data_ch03 + 0.45 * data_ch01
        B = data_ch01

        # Minimuns, Maximuns and Gamma
        Rmin = 0.0
        Rmax = 1.0
        Gmin = 0.0
        Gmax = 1.0
        Bmin = 0.0
        Bmax = 1.0

        R[R > Rmax] = Rmax
        R[R < Rmin] = Rmin
        G[G > Gmax] = Gmax
        G[G < Gmin] = Gmin
        B[B > Bmax] = Bmax
        B[B < Bmin] = Bmin

        gamma_R = 1
        gamma_G = 1
        gamma_B = 1

        # Normalize the data
        R = ((R - Rmin) / (Rmax - Rmin)) ** (1 / gamma_R)
        G = ((G - Gmin) / (Gmax - Gmin)) ** (1 / gamma_G)
        B = ((B - Bmin) / (Bmax - Bmin)) ** (1 / gamma_B)

        # Create the RGB
        RGB = np.stack([R, G, B], axis=2)

        # Eliminate values outside the globe
        mask = (RGB == [R[0, 0], G[0, 0], B[0, 0]]).all(axis=2)
        RGB[mask] = np.nan

    else:
        return False

    # Formatando data para plotar na imagem e salvar o arquivo
    date = (datetime.datetime.strptime(dtime_ch01, '%Y%j%H%M%S'))
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')
    date_file = date.strftime('%Y%m%d_%H%M%S')

    # Formatando a descricao a ser plotada na imagem
    description = f' GOES-{satellite_ch01} Natural True Color {date_img}'
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

    # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
    img_extent = [extent[0], extent[2], extent[1], extent[3]]  # Min lon, Max lon, Min lat, Max lat

    # Plotando a imagem
    ax.imshow(RGB, origin='upper', extent=img_extent)

    # Adicionando descricao da imagem
    # Criando novos eixos de acordo com a posicao da imagem
    cax0 = fig.add_axes([ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.0135, ax.get_position().width - 0.003, 0.0125])
    cax0.patch.set_color('black')  # Alterando a cor do novo eixo
    cax0.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
    cax0.text(0.85, 0.13, institution, color='yellow', size=10)  # Adicionando texto
    cax0.xaxis.set_visible(False)  # Removendo rotulos do eixo X
    cax0.yaxis.set_visible(False)  # Removendo rotulos do eixo Y

    # Adicionando os logos
    logo_noaa = plt.imread(dir_logos + 'NOAA_Logo.png')  # Lendo o arquivo do logo
    logo_goes = plt.imread(dir_logos + 'GOES_Logo.png')  # Lendo o arquivo do logo
    logo_cepagri = plt.imread(dir_logos + 'CEPAGRI-Logo.png')  # Lendo o arquivo do logo
    fig.figimage(logo_noaa, 32, 203, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_goes, 10, 120, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_cepagri, 10, 40, zorder=3, alpha=0.8, origin='upper')  # Plotando logo
    
    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}{rgb_type}/{rgb_type}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
    
    # Fecha a janela para limpar a memoria
    plt.close()
    # Realiza o log do calculo do tempo de processamento da imagem
    logging.info(f'{ch01[0:59].replace("M6C01", "M6C0*")} - {v_extent} - {str(round(time.time() - processing_start_time, 4))} segundos')


def processing(bands, p_br, p_sp, dir_in): 
    # Cria lista vazia para controle do processamento paralelo
    process_br = []
    # Cria lista vazia para controle processamento paralelo
    process_sp = []
    
    # Se a variavel de controle de processamento do brasil for True, realiza o processamento das bandas do ABI 1-16
    if p_br:
        logging.info('')
        logging.info('PROCESSANDO IMAGENS "BR"...')
        # Contador para processamento nas 16 bandas
        for x in range(1, 17):
            b = str(x).zfill(2) # Transforma o inteiro contador em string e com 2 digitos
            # Se a banda tiver novas imagens para o dia:
            if bands[b]:
                # Imagens para processamento
                old_bands = openOld()
                # Tentando Processar
                try:
                    # Cria o processo com a funcao de processamento
                    process = Process(target=process_band_cmi, args=(f'{dir_in}band{b}/{old_bands[b]}', b, "br"))
                    # Adiciona o processo na lista de controle do processamento paralelo
                    process_br.append(process)
                    # Inicia o processo
                    process.start()
                # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
                except:
                    # Realiza o log do erro
                    logging.info(f'Erro Arquivo - {old_bands[b]}')
                    # Remove a imagem com erro de processamento
                    os.remove(f'{dir_in}band{b}/{old_bands[b]}')
            else:
                continue
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in process_br:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()
        # Limpa lista vazia para controle do processamento paralelo
        process_br = []

    if p_sp:
        logging.info('')
        logging.info('PROCESSANDO IMAGENS "SP"...')
        # Contador para processamento nas 16 bandas
        for x in range(1, 17):
            b = str(x).zfill(2) # Transforma o inteiro contador em string e com 2 digitos
            # Se a banda tiver novas imagens para o dia:
            if bands[b]:
                # Imagens para processamento
                old_bands = openOld()
                # Tentando Processar
                try:
                    # Cria o processo com a funcao de processamento
                    process = Process(target=process_band_cmi, args=(f'{dir_in}band{b}/{old_bands[b]}', b, 'sp'))
                    # Adiciona o processo na lista de controle do processamento paralelo
                    process_sp.append(process)
                    # Inicia o processo
                    process.start()
                # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
                except:
                    # Realiza o log do erro
                    logging.info(f'Erro Arquivo - {old_bands[b]}')
                    # Remove a imagem com erro de processamento
                    os.remove(f'{dir_in}band{b}/{old_bands[b]}')
            else:
                continue
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in process_sp:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()
        # Limpa lista vazia para controle do processamento paralelo
        process_sp = []

# Checagem se e possivel gerar imagem TrueColor
    if bands['17']:
        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS TRUECOLOR "BR"...')
            
            ch01 = old_bands['01']
            ch02 = old_bands['02']
            ch03 = old_bands['03']
            
            # Montando dicionario de argumentos
            kwargs = {'ch01': f'{dir_in}band01/{ch01.replace(".nc", "_reproj_br.nc")}', 'ch02': f'{dir_in}band02/{ch02.replace(".nc", "_reproj_br.nc")}', 
                      'ch03': f'{dir_in}band03/{ch03.replace(".nc", "_reproj_br.nc")}'}
            
            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process = Process(target=process_band_rgb, args=("truecolor", "br"), kwargs=kwargs)
                # Adiciona o processo na lista de controle do processamento paralelo
                process_br.append(process)
                # Inicia o processo
                process.start()
                
            # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
            except:
                # Realiza o log do erro
                logging.info("Erro Arquivo  truecolor")
                
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in process_br:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()
        # Limpa lista vazia para controle do processamento paralelo
        process_br = []
        
        
        # Se a variavel de controle de processamento sp for True, realiza o processamento
        if p_sp:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS TRUECOLOR "SP"...')
            
            ch01 = old_bands['01']
            ch02 = old_bands['02']
            ch03 = old_bands['03']
            
            # Montando dicionario de argumentos
            kwargs = {'ch01': f'{dir_in}band01/{ch01.replace(".nc", "_reproj_sp.nc")}', 'ch02': f'{dir_in}band02/{ch02.replace(".nc", "_reproj_sp.nc")}', 
                      'ch03': f'{dir_in}band03/{ch03.replace(".nc", "_reproj_sp.nc")}'}
            
            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process = Process(target=process_band_rgb, args=("truecolor", "sp"), kwargs=kwargs)
                # Adiciona o processo na lista de controle do processamento paralelo
                process_sp.append(process)
                # Inicia o processo
                process.start()
                
            # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
            except:
                # Realiza o log do erro
                logging.info("Erro Arquivo  truecolor")
                
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in process_sp:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()
        # Limpa lista vazia para controle do processamento paralelo
        process_sp = []