import matplotlib
matplotlib.use('Agg')# Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt  # Plotagem de resultados, textos, logos, etc.
from matplotlib import cm  # Utilitario para paletas de cores
import cartopy  # Inserir mapas, shapefiles, paralelos, meridianos, latitudes, longitudes, etc.
import cartopy.crs as ccrs  # Utilitario para sistemas de referência e coordenadas
import cartopy.io.shapereader as shpreader  # Utilitario para leitura de shapefiles
import cartopy.feature as cfeature # features
from osgeo import gdal  # Utilitario para a biblioteca GDAL
from osgeo import osr  # Utilitario para a biblioteca GDAL
from netCDF4 import Dataset  # Utilitario para a biblioteca NetCDF4
import numpy as np  # Suporte para arrays e matrizes multidimensionais, com diversas funções matemáticas para trabalhar com estas estruturas
import time  # Utilitario para trabalhar com tempos
import os  # Utilitario para trabalhar com chamadas de sistema
import logging  # Utilitario para criar os logs
from multiprocessing import Process  # Utilitario para multiprocessamento
import re # Utilitario para trabalhar com expressoes regulares
import json
from modules.dirs import get_dirs
from modules.utilities import load_cpt  # Funcao para ler as paletas de cores de arquivos CPT
from modules.utilities import download_prod
from shapely.geometry import Point
import shapely.geometry as sg
from shutil import copyfile  # Utilitario para copia de arquivos
from pyorbital import astronomy
from pyspectral.rayleigh import Rayleigh                     # Correção atmosférica no espectro visível 
from pyorbital.astronomy import get_alt_az
from pyorbital.orbital import get_observer_look
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from datetime import timedelta, datetime
from modules.remap import remap

# Configurar o NumPy para ignorar os avisos
np.seterr(invalid='ignore')

osr.DontUseExceptions()

gdal.PushErrorHandler('CPLQuietErrorHandler')   # Ignore GDAL warnings

# ============================================# Diretórios ========================================= #
dirs = get_dirs()

# Acessando os diretórios usando as chaves do dicionário
dir_in = dirs['dir_in']
dir_out = dirs['dir_out']
dir_shapefiles = dirs['dir_shapefiles']
dir_colortables = dirs['dir_colortables']
dir_logos = dirs['dir_logos']
dir_temp = dirs['dir_temp']
arq_log = dirs['arq_log']
dir_maps = dirs['dir_maps']
# ============================================# Diretórios ========================================= #


def filtrar_imagens_por_intervalo(images, ch13):
    # Extrai a data e hora da string 'ch13' e define um intervalo de 9 minutos e 59 segundos a partir dela.
    glm_list = [] 
    ch13_data = (datetime.strptime(ch13[ch13.find("M6C13_G16_s") + 11:ch13.find("_e") - 1], '%Y%j%H%M%S'))
    date_ini = datetime(ch13_data.year, ch13_data.month, ch13_data.day, ch13_data.hour, ch13_data.minute)
    date_end = datetime(ch13_data.year, ch13_data.month, ch13_data.day, ch13_data.hour, ch13_data.minute) + timedelta(minutes=9, seconds=59)
    # Percorre a lista de nomes de imagens e verifica se a data e hora de cada imagem estão dentro do intervalo.
    for x in images:
        xtime = (datetime.strptime(x[x.find("GLM-L2-LCFA_G16_s") + 17:x.find("_e") - 1], '%Y%j%H%M%S'))
        if date_ini <= xtime <= date_end:
            glm_list.append(x)
        else:
            continue
    
    return glm_list


def area_para_recorte(v_extent):
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
    return extent, resolution


def adicionando_shapefile(v_extent, ax, colors='cyan'):
    if v_extent == 'br':
        # Adicionando o shapefile dos estados brasileiros
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil/BR_UF_2020.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=colors, facecolor='none', linewidth=0.7, zorder=8)
    elif v_extent == 'sp':
        # Adicionando o shapefile dos estados brasileiros e cidade de campinas
        # https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2020/Brasil/BR/BR_UF_2020.zip
        shapefile = list(shpreader.Reader(dir_shapefiles + 'brasil/BR_UF_2020.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=colors, facecolor='none', linewidth=0.7, zorder=8)
        shapefile = list(shpreader.Reader(dir_shapefiles + 'campinas/campinas.shp').geometries())
        ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor=colors, facecolor='none', linewidth=1, zorder=8)


def adicionando_linhas(ax, colors='cyan'):
    # Adicionando  linhas dos litorais
    ax.coastlines(resolution='10m', color=colors, linewidth=0.5, zorder=8)
    # Adicionando  linhas das fronteiras
    ax.add_feature(cartopy.feature.BORDERS, edgecolor=colors, linewidth=0.5, zorder=8)
    # Adicionando  paralelos e meridianos
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=0.7, linestyle='--', linewidth=0.2, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), zorder=8)
    gl.top_labels = False
    gl.right_labels = False


def adicionando_descricao_imagem(description, institution, ax, fig, cruz=False):
    cax1 = fig.add_axes([ax.get_position().x0 + 0.003, ax.get_position().y0 - 0.026, ax.get_position().width - 0.003, 0.0125])
    cax1.patch.set_color('black')  # Alterando a cor do novo eixo
    cax1.text(0, 0.13, description, color='white', size=10)  # Adicionando texto
    if cruz:
        cruzStr = '+'
        cax1.text(0.192, 0.13, cruzStr, color='red', size=12)  # Adicionando símbolo "+"
    cax1.text(0.85, 0.13, institution, color='yellow', size=10)  # Adicionando texto
    cax1.xaxis.set_visible(False)  # Removendo rótulos do eixo X
    cax1.yaxis.set_visible(False)  # Removendo rótulos do eixo Y


def calculating_lons_lats(date, extent, data_ch01, data_ch02, data_ch03):
        
    year = date.strftime('%Y')
    month = date.strftime('%m')
    day = date.strftime('%d')
    hour = date.strftime('%H')
    minutes = date.strftime('%M')
    
    # Criar as latitudes e longitudes com base na extensão
    lat = np.linspace(extent[3], extent[1], data_ch01.shape[0])
    lon = np.linspace(extent[0], extent[2], data_ch01.shape[1])
    xx,yy = np.meshgrid(lon,lat)
    lons = xx.reshape(data_ch01.shape[0], data_ch01.shape[1])
    lats = yy.reshape(data_ch01.shape[0], data_ch01.shape[1])

    # Obter o ano, mês, dia, hora e minuto para aplicar a correção zenital
    utc_time = datetime(int(year), int(month), int(day), int(hour), int(minutes))
    sun_zenith = np.zeros((data_ch01.shape[0], data_ch01.shape[1]))
    sun_zenith = astronomy.sun_zenith_angle(utc_time, lons, lats)

    # Aplicar a correção zenital do sol
    data_ch01 = (data_ch01)/(np.cos(np.deg2rad(sun_zenith)))
    data_ch02 = (data_ch02)/(np.cos(np.deg2rad(sun_zenith)))
    data_ch03 = (data_ch03)/(np.cos(np.deg2rad(sun_zenith)))

    return utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03


def applying_rayleigh_correction(file_ch01, utc_time, lons, lats, sun_zenith, data_ch01, data_ch02, longitude):
    # Altitude do satélite
    sat_h = file_ch01.variables['goes_imager_projection'].perspective_point_height

    sunalt, suna = get_alt_az(utc_time, lons, lats)
    suna = np.rad2deg(suna)
    #sata, satel = get_observer_look(sat_lon, sat_lat, sat_alt, vis.attrs['start_time'], lons, lats, 0)
    sata, satel = get_observer_look(longitude, 0.0, sat_h, utc_time, lons, lats, 0)
    satz = 90 - satel

    # Correção de Rayleigh
    atmosphere = 'us-standard'
    aerosol_type = 'rayleigh_only'
    corrector = Rayleigh('GOES-16', 'abi', atmosphere=atmosphere, aerosol_type=aerosol_type)

    sata = sata % 360.
    suna = suna % 360.
    ssadiff = np.absolute(suna - sata)
    ssadiff = np.minimum(ssadiff, 360 - ssadiff)

    red = data_ch02 * 100

    refl_cor_band_c01 = corrector.get_reflectance(sun_zenith, satz, ssadiff, 'C01', redband=red)
    data_ch01 = data_ch01 - (refl_cor_band_c01 / 100)

    refl_cor_band_c02 = corrector.get_reflectance(sun_zenith, satz, ssadiff, 'C02', redband=red)
    data_ch02 = data_ch02 - (refl_cor_band_c02 / 100)

    return data_ch01, data_ch02 


def apply_cira_stretch(band_data):
    
    log_root = np.log10(0.0223)
    denom = (1.0 - log_root) * 0.75
    band_data *= 0.01
    band_data = band_data.clip(np.finfo(float).eps)
    band_data = np.log10(band_data)
    band_data -= log_root
    band_data /= denom
    return 1 + band_data


def adicionando_logos(fig):
    # Adicionando os logos
    logo_noaa = plt.imread(dir_logos + 'NOAA_Logo.png')  # Lendo o arquivo do logo
    logo_goes = plt.imread(dir_logos + 'GOES_Logo.png')  # Lendo o arquivo do logo
    logo_cepagri = plt.imread(dir_logos + 'CEPAGRI-Logo.png')  # Lendo o arquivo do logo
    fig.figimage(logo_noaa, 32, 240, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_goes, 10, 160, zorder=3, alpha=0.6, origin='upper')  # Plotando logo
    fig.figimage(logo_cepagri, 10, 80, zorder=3, alpha=0.8, origin='upper')  # Plotando logo


def abrir_old_json(dir_main):
    with open(f'{dir_main}old_bands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
        return oldImages


def save_txt(array, nome_arquivo_txt):
    # Checa se a matriz é vazia
    if len(array) == 0:
        logging.info(f'{nome_arquivo_txt} vazia')
        pass
    else:
        # Criando nome do arquivo e diretório -- Mudar as barras para Linux
        with open(f"{dir_out}fdcf/{nome_arquivo_txt}.txt", 'w') as file:  # tomar cuidado se não ele fica criando diretorio
            for valor in array:
                valor = f"{valor[0]},{valor[1]}\n"
                file.write(valor)


def save_log_erro(array_errors, nome_arquivo_txt):
    # Checa se a matriz é vazia
    if len(array_errors) == 0:
        pass
    else:
        # Criando nome do arquivo e diretório 
        with open(f"{dir_out}fdcf/{nome_arquivo_txt}.txt", 'w') as file:
            for valor in array_errors:
                erro = f"{valor}\n"
                file.write(erro)


def fdcf_tabela_hot_spots(date, ax):    
    # Cria uma lista com os itens no diretorio temp que sao arquivos e se encaixa na expressao regular "^fdcf_.+_.+_br.txt$"
    fdcf_list = [name for name in os.listdir(f'{dir_out}fdcf') if os.path.isfile(os.path.join(f'{dir_out}fdcf', name)) and re.match(f'^fdcf_{date.strftime("%Y%m%d")}_.+_br.txt$', name)]

    # Cria matriz diária de pontos e log. de erros e faz o ‘loop’ nos arquivos do diretório.
    matriz_diaria = []

    ## Captura a lista com os "contornos" dos estados para separar o total de pontos individual
    geometry_estados = list(shpreader.Reader(dir_shapefiles + "divisao_estados/gadm40_BRA_1.shp").geometries())

    ## Os 26 estados estão em ordem alfabética no arquivo shapefile e o indice 27 armazena o total diário
    lista_totalpixels_uf = [0 for i in range(28)]
    log_erro = []
    for name in fdcf_list:
        try:
            with open(f'{dir_out}fdcf/{name}', 'r') as file:
                for linha in file.readlines():
                    linha = linha.split(',')
                    p = [float(linha[0]), float(linha[1])]  # Ponto em lat,lon
                    ## Checkagem para evitar contagem de pontos duplicados
                    if not (p in matriz_diaria):   
                        ax.plot(float(linha[1]), float(linha[0]), 'r+', ms=2.5, transform=ccrs.PlateCarree())  # plotando em lon,lat
                        matriz_diaria.append(p)
                        ## Contagem do total por estado
                        for j in range(27):
                            estado = geometry_estados[j]
                            if estado.covers(Point(float(linha[1]),float(linha[0]))):
                                lista_totalpixels_uf[j]+=1
                                break
                    else:
                        continue
        except Exception as erro:
            logging.info(f'Erro no processamento da matriz diária: {erro}')
            log_erro.append(f'{name} - Error: {erro}')
            pass

    # Calcula a soma do total de focos incendios diário
    lista_totalpixels_uf[27]+= np.sum(lista_totalpixels_uf)

    ### Manipulação dos dados para tabela
    lista_siglas_UF = ['AC','AL','AP','AM','BA','CE','DF','ES','GO','MA','MT','MS','MG','PA','PB','PR','PE','PI','RJ','RN',
                'RS','RO','RR','SC','SP','SE','TO','Total']
    columns = ('UF', 'Hot Spot')
    data = []
    for i in range(28):
        table_value = (lista_siglas_UF[i],lista_totalpixels_uf[i])
        data.append(table_value)
    
    ## Cria a tabela com os dados separados anteriormente
    tabela= matplotlib.table.table(ax =ax,cellText=data,colLabels=columns,cellLoc ='center',loc='lower right',rowLoc='center',
                            colLoc='center')
    tabela.auto_set_column_width(col=list(range(len(data))))

    # Cria os nomes dos arquivos diários e salva no "directory".

    logging.info("save_txt: matriz_diaria")
    save_txt(matriz_diaria, f'fdcf_{date.strftime("%Y%m%d")}_br')
    save_log_erro(log_erro, f'fdcf_{date.strftime("%Y%m%d")}_errors')

    return fdcf_list


def reproject(reproj_file, reproj_var, reproj_extent, reproj_resolution, dir_in):

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
    extent, resolution = area_para_recorte(v_extent)
        
    # Reprojetando imagem CMI e recebendo data/hora da imagem, satelite e caminho absoluto do arquivo reprojetado
    dtime, satellite, reproject_band = reproject(file, file_var, v_extent, resolution, dir_in)

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
    date = (datetime.strptime(dtime, '%Y-%m-%dT%H:%M:%S.%fZ'))
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
    
    # Adicionando o shapefile dos estados brasileiros
    adicionando_shapefile(v_extent, ax)

    # Adicionando  linhas dos litorais
    adicionando_linhas(ax)

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
    adicionando_descricao_imagem(description, institution, ax, fig)

    # Adicionando os logos
    adicionando_logos(fig)

    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}band{ch}/band{ch}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
    
    # Fecha a janela para limpar a memoria
    plt.close()
    
    # Realiza o log do calculo do tempo de processamento da imagem
    logging.info(f'{file} - {v_extent} - {str(round(time.time() - processing_start_time, 4))} segundos')


def process_truecolor(rgb_type, v_extent, ch01, ch02, ch03, ch13):
    global dir_maps
    start = time.time()  
    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------
    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)
    # Variable to remap
    variable = "CMI"
    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------

    # Lê a imagem da banda 01
    file_ch01 = Dataset(ch01)
    # Lê a longitude central
    longitude = file_ch01.variables['goes_imager_projection'].longitude_of_projection_origin
    # Lê a data do arquivo
    add_seconds = int(file_ch01.variables['time_bounds'][0])
    date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
    date_file = date.strftime('%Y%m%d_%H%M%S')
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')

    #------------------------------------------------------------------------------------------------------#
    #-------------------------------------------Reprojetando----------------------------------------------#
    #------------------------------------------------------------------------------------------------------#
    # reprojetando band 01
    grid = remap(ch01, variable, extent, resolution)
    # Lê o retorno da função
    data_ch01 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # reprojetando band 02
    grid = remap(ch02, variable, extent, resolution)
    # Lê o retorno da função
    data_ch02 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # reprojetando band 03
    grid = remap(ch03, variable, extent, resolution)
    # Lê o retorno da função 
    data_ch03 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
     # reprojetando band 13
    grid = remap(ch13, variable, extent, resolution)
    # Lê o retorno da função
    data_ch13 = grid.ReadAsArray()
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    
    # Calculando correção zenith
    utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03 = calculating_lons_lats(date, extent, data_ch01, data_ch02, data_ch03)

    # Aplicando a correção de Rayleigh
    data_ch01, data_ch02 = applying_rayleigh_correction(file_ch01, utc_time, lons, lats, sun_zenith, data_ch01, data_ch02, longitude)

    # Calculando as cores verdadeiras (True color)
    R = data_ch02
    G = (data_ch01 + data_ch02) / 2 * 0.93 + 0.07 * data_ch03 
    B = data_ch01
    # Aplicando o estiramento CIRA
    R = apply_cira_stretch(R)
    G = apply_cira_stretch(G)
    B = apply_cira_stretch(B)

    # Create the RGB
    RGB = np.stack([R, G, B], axis=2)		

    # If zenith angle is greater than 85°, the composite pixel is zero
    RGB[sun_zenith > 85] = 0
    # Create the mask for the regions with zero
    mask = (RGB == [0.0,0.0,0.0]).all(axis=2)
    # Apply the mask to overwrite the pixels
    RGB[mask] = [0,0,0]
    # Create the fading transparency between the regions with the sun zenith angle of 75° and 85°
    alphas = sun_zenith / 100
    min_sun_angle = 0.75
    max_sun_angle = 0.85
    # Normalize the transparency mask
    alphas = ((alphas - max_sun_angle) / (min_sun_angle - max_sun_angle))
    RGB = np.dstack((RGB, alphas))
    
    if v_extent == 'sp':
        raster = gdal.Open(f'{dir_maps}BlackMarble_2016_B2_geo.tif')
    else:
        raster = gdal.Open(f'{dir_maps}BlackMarble_2016_01deg_geo.tif')
        
    ulx, xres, xskew, uly, yskew, yres = raster.GetGeoTransform()
    lrx = ulx + (raster.RasterXSize * xres)
    lry = uly + (raster.RasterYSize * yres)
    corners = [ulx, lry, lrx, uly]
    min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]
    raster = gdal.Translate('raster.tif', raster, projWin = [min_lon, max_lat, max_lon, min_lat])
    
    #lendo o RGB 
    array = raster.ReadAsArray()
    R_night = array[0,:,:].astype(float) / 255
    G_night = array[1,:,:].astype(float) / 255
    B_night = array[2,:,:].astype(float) / 255
    
    R_night[R_night==4] = 0
    G_night[G_night==5] = 0
    B_night[B_night==15] = 0
    
    #
    rgb_night = np.stack([R_night,G_night,B_night], axis=2)
    
    # Area de recorte
    img_extent = [extent[0], extent[2], extent[1], extent[3]]  
    
    # Remove o arquivo.tif
    os.remove('raster.tif')
    
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    
    # Formatando a descricao a ser plotada na imagem
    description = f' GOES-16 Natural True Color + Black Marble {date_img}'
    institution = "CEPAGRI - UNICAMP"

    d_p_i = 150
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')
    
    # Utilizando projecao geoestacionaria no cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    #band13
    data1 = data_ch13
    data1 = np.maximum(data1, 90)
    data1 = np.minimum(data1, 313)
    data1 = (data1-90)/(313-90)
    data1 = 1 - data1
    
    # Plotando a imagem night
    ax.imshow(rgb_night, extent=img_extent)

    # Plotando band13
    ax.imshow(data1, cmap='gray', vmin=0.1, vmax=0.25, alpha = 0.3, origin='upper', extent=img_extent)

    # Plotando a imagem  # TrueColor
    ax.imshow(RGB, origin='upper', extent=img_extent)

    # Adicionando o shapefile dos estados brasileiros
    adicionando_shapefile(v_extent, ax)

    # Adicionando  linhas dos litorais
    adicionando_linhas(ax)

    # Adicionando descricao da imagem
    adicionando_descricao_imagem(description, institution, ax, fig)

    # Adicionando os logos
    adicionando_logos(fig)

    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}{rgb_type}/{rgb_type}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)

    # Fecha a janela para limpar a memoria
    plt.close()
    
    logging.info(f'True Color - Tempo de Processamento: {round((time.time() - start), 2)} segundos.')
    

def process_rrqpef(rrqpef, ch13, v_extent):
    global dir_in, dir_shapefiles, dir_colortables, dir_logos, dir_out
    file_var = 'RRQPE'
    # Captura a hora para contagem do tempo de processamento da imagem
    processing_start_time = time.time()
    
    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)

    # Reprojetando imagem CMI e recebendo data/hora da imagem, satelite e caminho absoluto do arquivo reprojetado
    dtime, satellite, reproject_rrqpef = reproject(rrqpef, file_var, v_extent, resolution, dir_in)
    
    data = reproject_rrqpef.ReadAsArray()
    reproject_rrqpef = None
    del reproject_rrqpef
    # Convert from int16 to uint16
    data = data.astype(np.float64)
    data[data == max(data[0])] = np.nan
    data[data == min(data[0])] = np.nan

    reproject_ch13 = Dataset(ch13)
    data_ch13 = reproject_ch13.variables['Band1'][:]

    # Formatando data para plotar na imagem e salvar o arquivo
    date = (datetime.strptime(dtime, '%Y-%m-%dT%H:%M:%S.%fZ'))
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')
    date_file = date.strftime('%Y%m%d_%H%M%S')

    # Formatando a descricao a ser plotada na imagem
    description = f' GOES-{satellite} Rainfall Rate (Quantitative Precipitation Estimate) mm/h {date_img}'
    institution = "CEPAGRI - UNICAMP"

    # Definindo tamanho da imagem de saida
    d_p_i = 150
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

    # Utilizando projecao geoestacionaria no cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

    # Adicionando o shapefile dos estados brasileiros
    adicionando_shapefile(v_extent, ax)

    # Adicionando  linhas dos litorais
    adicionando_linhas(ax)

    # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
    img_extent = [extent[0], extent[2], extent[1], extent[3]]  # Min lon, Max lon, Min lat, Max lat

    # Plotando base
    ax.imshow(data_ch13, origin='upper', cmap='gray_r', extent=img_extent)

    # Plotando a imagem
    img = ax.imshow(data, vmin=0, vmax=150, cmap='jet', origin='upper', extent=img_extent)

    # Criando novos eixos de acordo com a posicao da imagem
    cax0 = fig.add_axes([ax.get_position().x0 + 0.0025, ax.get_position().y0 - 0.01325, ax.get_position().width - 0.0025, 0.0125])
    cb = plt.colorbar(img, orientation="horizontal", cax=cax0)
    cb.ax.tick_params(axis='x', colors='white', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
    cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
    cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
    cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de cores

    # Adicionando descricao da imagem
    adicionando_descricao_imagem(description, institution, ax, fig)

    # Adicionando os logos
    adicionando_logos(fig)

    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}rrqpef/rrqpef_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
    # Fecha a janela para limpar a memoria
    plt.close()
    # Realiza o log do calculo do tempo de processamento da imagem
    logging.info(f'{rrqpef} - {v_extent} - {str(round(time.time() - processing_start_time, 4))} segundos')


def process_glm(ch13, glm_list, v_extent):
    global dir_in, dir_shapefiles, dir_colortables, dir_logos, dir_out
    
    # Captura a hora para contagem do tempo de processamento da imagem
    processing_start_time = time.time()
    
    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)

    # Lendo imagem CMI reprojetada
    reproject_ch13 = Dataset(ch13)
    data_ch13 = reproject_ch13.variables['Band1'][:]

    # Coletando do nome da imagem o satelite e a data/hora
    satellite = (glm_list[0][glm_list[0].find("GLM-L2-LCFA_G") + 13:glm_list[0].find("_s")])
    dtime = (glm_list[0][glm_list[0].find("GLM-L2-LCFA_G16_s") + 17:glm_list[0].find("_e") - 1])

    # Initialize arrays for latitude, longitude, and event energy
    lats = np.array([])
    lons = np.array([])
    energies = np.array([])

    for g in glm_list:
        # Lendo imagem CMI reprojetada
        glm = Dataset(f'{dir_in}glm/{g}')
        # Append lats / longs / event energies
        lats = np.append(lats, glm.variables['event_lat'][:])
        lons = np.append(lons, glm.variables['event_lon'][:])
        energies = np.append(energies, glm.variables['event_energy'][:])

    # Stack and transpose the lat lons
    values = np.vstack((lats, lons)).T
    # Get the counts
    points, counts = np.unique(values, axis=0, return_counts=True)
    # Get the counts indices
    idx = counts.argsort()

    # Formatando data para plotar na imagem e salvar o arquivo
    date = (datetime.strptime(dtime, '%Y%j%H%M%S'))
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')
    date_file = date.strftime('%Y%m%d_%H%M%S')

    # Formatando a descricao a ser plotada na imagem
    description = f' GOES-{satellite} GLM Density {date_img}'
    institution = "CEPAGRI - UNICAMP"

    # Definindo tamanho da imagem de saida
    d_p_i = 150
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

    # Use the Geostationary projection in cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Define the image extent
    ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

    # Adicionando o shapefile dos estados brasileiros
    adicionando_shapefile(v_extent, ax)
    
    # Adicionando  linhas dos litorais
    adicionando_linhas(ax)

    # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
    img_extent = [extent[0], extent[2], extent[1], extent[3]]  # Min lon, Max lon, Min lat, Max lat

    # Plotando base
    ax.imshow(data_ch13, origin='upper', cmap='gray_r', extent=img_extent)

    # Plot the GLM Data
    glm = plt.scatter(points[idx, 1], points[idx, 0], vmin=0, vmax=1000, s=counts[idx] * 0.1, c=counts[idx], cmap="jet", zorder=2)

    cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
    cb = plt.colorbar(glm, orientation="horizontal", cax=cax0, ticks=[200, 400, 600, 800])
    cb.ax.set_xticklabels(['200', '400', '600', '800'])
    cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
    cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
    cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
    cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de cores

    # Adicionando descricao da imagem
    adicionando_descricao_imagem(description, institution, ax, fig)

    # Adicionando os logos
    adicionando_logos(fig)
    
    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}glm/glm_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
    # Fecha a janela para limpar a memoria
    plt.close()
    # Realiza o log do calculo do tempo de processamento da imagem
    logging.info(f'{dir_in}glm/{glm_list[0][0:31]}*.nc - {v_extent} - {str(round(time.time() - processing_start_time, 4))} segundos')


def process_ndvi(ndvi_diario, ch02, ch03, v_extent):
    
    global dir_shapefiles, dir_colortables, dir_logos, dir_in, dir_out
    
    # Captura a hora para contagem do tempo de processamento da imagem
    processing_start_time = time.time()
    
    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)

    # Lendo imagem CMI reprojetada
    reproject_ch02 = Dataset(ch02)
    reproject_ch03 = Dataset(ch03)

    # Coletando do nome da imagem a data/hora
    dtime_ch02 = ch02[ch02.find("M6C02_G16_s") + 11:ch02.find("_e") - 1]

    #//Hora para fazer o download Mask
    dtime_mask = ch02[ch02.find("M6C02_G16_s") + 11:ch02.find("_e") - 3]

    ## Transforma de data Juliana para data normal
    year = datetime.strptime(dtime_mask, '%Y%j%H%M').strftime('%Y')
    day_of_year = datetime.strptime(dtime_mask, '%Y%j%H%M').strftime('%d')
    month = datetime.strptime(dtime_mask, '%Y%j%H%M').strftime('%m')
    hora = datetime.strptime(dtime_mask, '%Y%j%H%M').strftime('%H')
    min = datetime.strptime(dtime_mask, '%Y%j%H%M').strftime('%M')
    yyyymmddhhmn = year+month+day_of_year+hora+min

    # Download arquivo mascara de nuvens
    logging.info("Download File Mask NDVI")
    file_mask = download_prod(yyyymmddhhmn,'ABI-L2-ACMF',f'{dir_in}clsm/')

    # Reprojetando o arquivo de nuvens
    reproject(f'{dir_in}clsm/{file_mask}.nc','BCM',v_extent, 4, dir_in)

    # Capturando o nome do arquivo reprojetado
    r_file = file_mask + f'_reproj_{v_extent}.nc'
    
    #  Abrindo o Dataset do arquivo da mascara
    reproject_mask = Dataset(f'{dir_in}clsm/{r_file}')


    data_ch02 = reproject_ch02.variables['Band1'][:]
    data_ch03 = reproject_ch03.variables['Band1'][:]
    data_mask = reproject_mask.variables['Band1'][:]

    reproject_ch02 = None
    del reproject_ch02
    
    reproject_ch03 = None
    del reproject_ch03

    reproject_mask = None
    del reproject_mask

    # NDVI Components
    NDVI = (data_ch03 - data_ch02) / (data_ch03 + data_ch02)

    # // Capturando somente os valores máximos de NDVI
    NDVI_fmax = np.zeros(np.shape(NDVI))
    NDVI_fmax[NDVI_fmax == 0 ] = np.nan

    NDVI_fmax = np.fmax(NDVI_fmax,NDVI)

    # // Retirando as nuvens do arquivo NDVI
    NDVI_fmax = np.where(data_mask==0,NDVI,data_mask)
    NDVI_fmax[NDVI_fmax == 1]= np.nan

    # Formatando data para plotar na imagem e salvar o arquivo
    date = (datetime.strptime(dtime_ch02, '%Y%j%H%M%S'))
    date_file = date.strftime('%Y%m%d_%H%M%S')

    # Salvando array NDVI
    logging.info('Salvando - NDVI')
    NDVI_fmax.dump(f'{dir_in}ndvi/ndvi_{date_file}_{v_extent}.npy')

    # Removendo os arquivos Clear Sky Mask baixados e reprojetado de nuvens
    logging.info('Removendo arquivo filemask do diretório')
    os.remove(f'{dir_in}clsm/{r_file}')
    os.remove(f'{dir_in}clsm/{file_mask}.nc')


    if ndvi_diario and datetime.now().isoweekday() == 6:
        
        # Captura a data atual e calculando data inicial e final
        date_now = datetime.now()
        date_ini = datetime(date_now.year, date_now.month, date_now.day, int(12), int(00))
        date_end = datetime(date_now.year, date_now.month, date_now.day, int(18), int(00))

        # Cria uma lista com os itens no diretorio temp que sao arquivos e se encaixa na expressao regular "^ndvi_.+_.+_br.npy$"
        ndvi_list = [name for name in os.listdir(f'{dir_in}ndvi/') if os.path.isfile(os.path.join(f'{dir_in}ndvi/', name)) and re.match('^ndvi_.+_.+_br.npy$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        ndvi_list.sort()
        teste = np.load(f'{dir_in}ndvi/{ndvi_list[0]}', allow_pickle=True)
        NDVI_fmax = np.zeros(teste.shape)
        ## // Seta os valores para nan para não haver problema entre fazer o máximo entre um valor negativo e zero 
        ## // max(i,v) = v se i#0 e v >i / v £ R, se i=0 e v<0 max(i,v) = i então se i = nan max(i,v)=v para v £ R
        NDVI_fmax[NDVI_fmax == 0] = np.nan
        del teste

        # Looping para ler e comparar os arquivos NDVI salvos do sabado
        for f in ndvi_list:
            # Captura a data do arquivo
            date = f[f.find("ndvi_") + 5:f.find("ndvi_") + 13]
            # Captura a hora do arquivo
            hour = f[f.find("ndvi_") + 14:f.find("ndvi_") + 20]
            # Monta a data/hora do arquivo
            date_file = (datetime.strptime(date + hour, '%Y%m%d%H%M%S'))
            # Se a data/hora do arquivo estiver dentro do limite de datas
            if date_ini <= date_file <= date_end + timedelta(minutes=1):
                # logging.info(f, date_file)
                # Le o arquivo NDVI
                NDVI_file = np.load(f'{dir_in}ndvi/{f}', allow_pickle=True)
                NDVI_file[NDVI_file > 0.999] = np.nan
                NDVI_file[NDVI_file < -0.999] = np.nan
                # Compara os arrays NDVI e retorna um novo array com os maiores valores de cada elemento
                NDVI_fmax = np.fmax(NDVI_fmax, NDVI_file)
                
            else:
                continue
                
        NDVI_fmax.dump(f'{dir_in}ndvi/ndvi_{date_ini.strftime("%Y%m%d")}_{date_end.strftime("%Y%m%d")}_br_fmax.npy')
        
        # Captura a data atual e calculando data inicial e final do acumulado da semana
        date_ini = datetime(date_now.year, date_now.month, date_now.day, int(23), int(59)) - timedelta(days=6, hours=23, minutes=59)
        date_end = datetime(date_now.year, date_now.month, date_now.day, int(23), int(59))

        # Cria uma lista com os itens no diretorio temp que sao arquivos e se encaixa na expressao regular "^ndvi_.+_.+_br.npy$"
        ndvi_list_fmax = [name for name in os.listdir(f'{dir_in}ndvi') if os.path.isfile(os.path.join(f'{dir_in}ndvi', name)) and re.match('^ndvi_.+_.+_br_fmax.npy$', name)]
        
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        ndvi_list_fmax.sort()
        teste = np.load(f'{dir_in}ndvi/{ndvi_list_fmax[0]}', allow_pickle=True)
        NDVI_fmax = np.zeros(teste.shape)
        NDVI_fmax[NDVI_fmax == 0] = np.nan
        del teste

        # Looping para ler e comparar os arquivos NDVI salvos da semana
        for f in range(0, len(ndvi_list_fmax)):
            # Captura a data do arquivo
            date_i = ndvi_list_fmax[f][ndvi_list_fmax[f].find("ndvi_") + 5:ndvi_list_fmax[f].find("ndvi_") + 13]
            date_e = ndvi_list_fmax[f][ndvi_list_fmax[f].find("ndvi_") + 14:ndvi_list_fmax[f].find("ndvi_") + 22]
            # Monta a data/hora do arquivo
            date_file_i = (datetime.strptime(date_i, '%Y%m%d'))
            date_file_e = (datetime.strptime(date_e, '%Y%m%d'))
            # Se a data/hora do arquivo estiver dentro do limite de datas
            if date_file_i == date_file_e and date_ini <= date_file_i <= date_end:
                # Le o arquivo NDVI
                NDVI_file_fmax = np.load(f'{dir_in}ndvi/{ndvi_list_fmax[f]}', allow_pickle=True)

                # Compara os arrays NDVI e retorna um novo array com os maiores valores de cada elemento
                NDVI_fmax = np.fmax(NDVI_fmax, NDVI_file_fmax)

            else:
                continue

        NDVI_fmax.dump(f'{dir_in}ndvi/ndvi_{date_ini.strftime("%Y%m%d")}_{date_end.strftime("%Y%m%d")}_br.npy')

        # Formatando a descricao a ser plotada na imagem
        description = f' GOES-16 COMPOSIÇÃO MÁXIMA DE VALOR NDVI NO PERÍODO {date_ini.strftime("%d-%b-%Y")}  -  {date_end.strftime("%d-%b-%Y")}'
        institution = "CEPAGRI - UNICAMP"

        # Choose the plot size (width x height, in inches)
        d_p_i = 150
        fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

        # Use the Geostationary projection in cartopy
        ax = plt.axes(projection=ccrs.PlateCarree())

        # Adicionando o shapefile dos estados brasileiros
        adicionando_shapefile(v_extent, ax, colors='dimgray')

        # Adicionando  linhas dos litorais
        adicionando_linhas(ax, colors='dimgray')
    
        # # Add an ocean mask
        ax.add_feature(cfeature.OCEAN, facecolor='cornflowerblue', zorder=3)
        
        # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
        img_extent = [extent[0], extent[2], extent[1], extent[3]]  # Min lon, Max lon, Min lat, Max lat

        colors = ["#0000FF", "#0030FF", "#0060FF", "#0090FF", "#00C0FF", "#003000", "#164A16", "#00C000", "#FFFF00", "#FF8000", "#FF0000"]
        my_cmap = cm.colors.LinearSegmentedColormap.from_list("", colors)  # Create a custom linear colormap

        # Plotando a imagem
        array_NDVI = np.load(f'{dir_in}ndvi/ndvi_{date_ini.strftime("%Y%m%d")}_{date_end.strftime("%Y%m%d")}_br.npy', allow_pickle=True)
        img = ax.imshow(array_NDVI, origin='upper', vmin=-1, vmax=1, cmap=my_cmap, extent=img_extent)

        cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
        cb = plt.colorbar(img, orientation="horizontal", cax=cax0, ticks=[-0.66, -0.33, 0, 0.33, 0.66])
        cb.ax.set_xticklabels(['-0.66', '-0.33', '0', '0.33', '0.66'])
        cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
        cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
        cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
        cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de cores

        # Adicionando descricao da imagem
        adicionando_descricao_imagem(description, institution, ax, fig)

        # Adicionando os logos
        adicionando_logos(fig)

        # Salvando a imagem de saida
        # plt.savefig(f'{dir_out}ndvi/ndvi_{date_end.strftime("%Y%m%d_%H%M%S")}_br.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
        plt.savefig(f'{dir_out}ndvi/ndvi_{date_end.strftime("%Y%m%d_%H%M%S")}_br.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
        # Fecha a janela para limpar a memoria
        plt.close()
        
    #   Cria o arquivo sem ser a imagem 
    elif ndvi_diario:
        # Captura a data atual e calculando data inicial e final
        date_now = datetime.now()
        date_ini = datetime(date_now.year, date_now.month, date_now.day, int(12), int(00))
        date_end = datetime(date_now.year, date_now.month, date_now.day, int(18), int(00))

        # Cria uma lista com os itens no diretorio temp que sao arquivos e se encaixa na expressao regular "^ndvi_.+_.+_br.npy$"
        ndvi_list = [name for name in os.listdir(f'{dir_in}ndvi') if os.path.isfile(os.path.join(f'{dir_in}ndvi', name)) and re.match('^ndvi_.+_.+_br.npy$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        ndvi_list.sort()
        teste = np.load(f'{dir_in}ndvi/{ndvi_list[0]}', allow_pickle=True)
        NDVI_fmax = np.zeros(teste.shape)
        ## // Seta os valores para nan para não haver problema entre fazer o máximo entre um valor negativo e zero 
        ## // max(i,v) = v se i#0 e v >i / v £ R, se i=0 e v<0 max(i,v) = i então se i = nan max(i,v)=v para v £ R
        NDVI_fmax[NDVI_fmax == 0] = np.nan
        del teste

        # Looping para ler e comparar os arquivos NDVI salvos do dia
        for f in ndvi_list:
            # Captura a data do arquivo
            date = f[f.find("ndvi_") + 5:f.find("ndvi_") + 13]
            # Captura a hora do arquivo
            hour = f[f.find("ndvi_") + 14:f.find("ndvi_") + 20]
            # Monta a data/hora do arquivo
            date_file = (datetime.strptime(date + hour, '%Y%m%d%H%M%S'))
            # Se a data/hora do arquivo estiver dentro do limite de datas
            if date_ini <= date_file <= date_end + timedelta(minutes=1):
                # logging.info(f, date_file)
                # Le o arquivo NDVI
                NDVI_file = np.load(f'{dir_in}ndvi/{f}', allow_pickle=True)
                NDVI_file[NDVI_file > 0.999] = np.nan
                NDVI_file[NDVI_file < -0.999] = np.nan
                # Compara os arrays NDVI e retorna um novo array com os maiores valores de cada elemento
                NDVI_fmax = np.fmax(NDVI_fmax, NDVI_file)
                # Remove o arquivo NDVI
                # os.remove(f'{dir_in}ndvi/{f}')
            else:
                continue


        NDVI_fmax.dump(f'{dir_in}ndvi/ndvi_{date_ini.strftime("%Y%m%d")}_{date_end.strftime("%Y%m%d")}_br_fmax.npy')

    # Realiza o log do calculo do tempo de processamento da imagem
    logging.info(f'{ch02[0:59].replace("M6C02", "M6C0*")} - {v_extent} - {str(round(time.time() - processing_start_time, 4))} segundos')


def process_fdcf(fdcf, ch01, ch02, ch03, v_extent, fdcf_diario):
    
    global dir_shapefiles, dir_colortables, dir_logos, dir_in, dir_out

    def degrees(file_id):
        
        proj_info = file_id.variables['goes_imager_projection']
        lon_origin = proj_info.longitude_of_projection_origin
        H = proj_info.perspective_point_height + proj_info.semi_major_axis
        r_eq = proj_info.semi_major_axis
        r_pol = proj_info.semi_minor_axis

        # Data info
        lat_rad_1d = file_id.variables['x'][:]
        lon_rad_1d = file_id.variables['y'][:]

        # Create meshgrid filled with radian angles
        lat_rad, lon_rad = np.meshgrid(lat_rad_1d, lon_rad_1d)

        # lat/lon calculus routine from satellite radian angle vectors
        lambda_0 = (lon_origin * np.pi) / 180.0

        a_var = np.power(np.sin(lat_rad), 2.0) + (np.power(np.cos(lat_rad), 2.0) * (
                np.power(np.cos(lon_rad), 2.0) + (((r_eq * r_eq) / (r_pol * r_pol)) * np.power(np.sin(lon_rad), 2.0))))
        b_var = -2.0 * H * np.cos(lat_rad) * np.cos(lon_rad)
        c_var = (H ** 2.0) - (r_eq ** 2.0)

        r_s = (-1.0 * b_var - np.sqrt((b_var ** 2) - (4.0 * a_var * c_var))) / (2.0 * a_var)

        s_x = r_s * np.cos(lat_rad) * np.cos(lon_rad)
        s_y = - r_s * np.sin(lat_rad)
        s_z = r_s * np.cos(lat_rad) * np.sin(lon_rad)

        lat = (180.0 / np.pi) * (np.arctan(((r_eq * r_eq) / (r_pol * r_pol)) * (s_z / np.sqrt(((H - s_x) * (H - s_x)) + (s_y * s_y)))))
        lon = (lambda_0 - np.arctan(s_y / (H - s_x))) * (180.0 / np.pi)

        return lat, lon

    file_var = 'FDCF'
    
    # Captura a hora para contagem do tempo de processamento da imagem
    processing_start_time = time.time()
    
    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)

    # Lendo imagem FDCF
    fire_mask = Dataset(fdcf)

    # Coletando do nome da imagem a data/hora
    dtime_fdcf = fdcf[fdcf.find("ABI-L2-FDCF-M6_G16_s") + 20:fdcf.find("_e") - 1]
    
    # Formatando data para plotar na imagem e salvar o arquivo
    date = (datetime.strptime(dtime_fdcf, '%Y%j%H%M%S'))
    date_img = date.strftime('%d-%b-%Y')
    date_file = date.strftime('%Y%m%d_%H%M%S')

    matriz_pixels_fogo = []
    fire_mask_values = fire_mask.variables['Mask'][:, :]
    selected_fires = (fire_mask_values == 10) | (fire_mask_values == 11) | (fire_mask_values == 13) | (
            fire_mask_values == 30) | (fire_mask_values == 33)

    lat, lon = degrees(fire_mask)
    
    # separando Latitudes e Longitudes dos pontos
    p_lat = lat[selected_fires]
    p_lon = lon[selected_fires]
    brasil = list(shpreader.Reader(dir_shapefiles + "divisao_estados/gadm36_BRA_0.shp").geometries())
    for i in range(len(p_lat)):
        if brasil[0].covers(Point(p_lon[i], p_lat[i])):
            p = (p_lat[i], p_lon[i])
            matriz_pixels_fogo.append(p)

    save_txt(matriz_pixels_fogo, f'fdcf_{date.strftime("%Y%m%d_%H%M%S")}_br')
    
    # Le o arquivo de controle de quantidade de pontos
    try:
        with open(f'{dir_temp}band21_control.txt', 'r') as fo:
            control = fo.readline()
            logging.info(f'tamanho control withopen: {int(control)}')
    except:
        control = 0
        logging.info(f'tamanho control except: {int(control)}')

    # Verifica se as ocorrencias de pontos é maior que as anteriores, se sim, armazena a quantidade e as imagens para gerar fundo
    logging.info(f'Len matriz_pixels_fogo:{len(matriz_pixels_fogo)}  int control:  {int(control)}')
    date_ini = datetime(date.year, date.month, date.day, int(13), int(00))
    date_end = datetime(date.year, date.month, date.day, int(18), int(00))
    
    # Pega a imagem com mais incidencia e salva
    if len(matriz_pixels_fogo) > int(control) and date_ini <= date <= date_end:
        copyfile(ch01, f'{dir_in}fdcf/ch01.nc')
        copyfile(ch02, f'{dir_in}fdcf/ch02.nc')
        copyfile(ch03, f'{dir_in}fdcf/ch03.nc')
        logging.info('Copiando bands 01 a 03 para fdcf')
        with open(f'{dir_temp}band21_control.txt', 'w') as fo:
            fo.write(str(len(matriz_pixels_fogo)))
    
    if fdcf_diario:
        # Reiniciar contagem para verificar imagem com maior quantidade de pontos no dia
        with open(f'{dir_temp}band21_control.txt', 'w') as fo:
            fo.write(str(0))

        ch01 = f'{dir_in}fdcf/ch01.nc'
        ch02 = f'{dir_in}fdcf/ch02.nc'
        ch03 = f'{dir_in}fdcf/ch03.nc'
        
        file_ch01 = Dataset(ch01)

        # Lê a longitude central
        longitude = file_ch01.variables['goes_imager_projection'].longitude_of_projection_origin

        # Lê a data do arquivo
        add_seconds = int(file_ch01.variables['time_bounds'][0])
        date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
        date_file = date.strftime('%Y%m%d_%H%M%S')
        date_img = date.strftime('%d-%b-%Y UTC')

        # Area de interesse para recorte
        extent, resolution = area_para_recorte(v_extent)
        variable = "CMI"
        
        #------------------------------------------------------------------------------------------------------#
        #-------------------------------------------Reprojetando----------------------------------------------#
        #------------------------------------------------------------------------------------------------------#
        # reprojetando band 01
        grid = remap(ch01, variable, extent, resolution)
        # Lê o retorno da função
        data_ch01 = grid.ReadAsArray()
        #------------------------------------------------------------------------------------------------------
        #------------------------------------------------------------------------------------------------------
        # reprojetando band 02
        grid = remap(ch02, variable, extent, resolution)
        # Lê o retorno da função
        data_ch02 = grid.ReadAsArray()
        #------------------------------------------------------------------------------------------------------
        #------------------------------------------------------------------------------------------------------
        # reprojetando band 03
        grid = remap(ch03, variable, extent, resolution)
        # Lê o retorno da função 
        data_ch03 = grid.ReadAsArray()
        #------------------------------------------------------------------------------------------------------
        #------------------------------------------------------------------------------------------------------

        #------------------------------------------------------------------------------------------------------
        # Calculando correção zenith
        utc_time, lats, lons, sun_zenith, data_ch01, data_ch02, data_ch03 = calculating_lons_lats(date, extent, data_ch01, data_ch02, data_ch03)

        # Aplicando a correção de Rayleigh
        data_ch01, data_ch02 = applying_rayleigh_correction(file_ch01, utc_time, lons, lats, sun_zenith, data_ch01, data_ch02, longitude)

        # Calculando as cores verdadeiras (True color)
        R = data_ch02
        G = (data_ch01 + data_ch02) / 2 * 0.93 + 0.07 * data_ch03 
        B = data_ch01

        # Aplicando o estiramento CIRA
        R = apply_cira_stretch(R)
        G = apply_cira_stretch(G)
        B = apply_cira_stretch(B)

        # Create the RGB
        RGB = np.stack([R, G, B], axis=2)		

        # Adicionando descricao da imagem.
        description = f"GOES-16 Natural True Color,    Fire Hot Spot em {date_img}"  # Esse espaço é necessário para adicionar o caractere na imagem
        institution = f'CEPAGRI - UNICAMP'

        # Definindo tamanho da imagem de saida.
        d_p_i = 150
        fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')
        
        # Utilizando projecao geoestacionaria no cartopy
        img_extent = [extent[0], extent[2], extent[1], extent[3]]
        ax = plt.axes(projection=ccrs.PlateCarree(), extent=img_extent)
        ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

        # Adicionando o shapefile dos estados brasileiros
        adicionando_shapefile(v_extent, ax)
        
        # Plotando a imagem RGB
        ax.imshow(RGB, origin='upper', extent=img_extent)

        # Adicionando descricao da imagem ##Apenas fire hot spot tem a cruz = true para descrição da imagem
        adicionando_descricao_imagem(description, institution, ax, fig, cruz=True)

        # Adicionando os logos
        adicionando_logos(fig)
        
        # Faz tabela de hot spots
        fdcf_list = fdcf_tabela_hot_spots(date, ax)
        
        # Salva imagem
        plt.savefig(f'{dir_out}fdcf/fdcf_{date_file}_br.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)
        # Fecha a janela para limpar a memoria
        plt.close()

        # Remove arquivos separados do dia
        logging.info('Removendo os arquivos fdcf do dia do diretório')
        for name in fdcf_list:
            os.remove(f'{dir_out}fdcf/{name}')

    # Realiza o log do calculo do tempo de processamento da imagem
    logging.info(f'{fdcf} - {v_extent} - {str(round(time.time() - processing_start_time, 4))} segundos')


def process_airmass(rgb_type, v_extent, path_ch08, path_ch10, path_ch12, path_ch13):
    global dir_out
    start = time.time()  

    # Read the file using the NetCDF library
    file_ch08 = Dataset(path_ch08)

    # Lê o identificador do satélite
    satellite = getattr(file_ch08, 'platform_ID')
    variable = "CMI"
    
    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)

    # Read the resolution
    band_resolution_km = getattr(file_ch08, 'spatial_resolution')
    band_resolution_km = float(band_resolution_km[:band_resolution_km.find("km")])

    # Getting the file time and date
    add_seconds = int(file_ch08.variables['time_bounds'][0])
    date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)
    date_file = date.strftime('%Y%m%d_%H%M%S')
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')


    #------------------------------------------------------------------------------------------------------#
    #-------------------------------------------Reprojetando----------------------------------------------#
    #------------------------------------------------------------------------------------------------------#
    # Call the reprojection funcion
    grid = remap(path_ch08, variable, extent, resolution)
    # Read the data returned by the function 
    data_ch08 = grid.ReadAsArray() - 273.15
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # Call the reprojection funcion
    grid = remap(path_ch10, variable, extent, resolution)
    # Read the data returned by the function 
    data_ch10 = grid.ReadAsArray() - 273.15
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # Call the reprojection funcion
    grid = remap(path_ch12, variable, extent, resolution)
    # Read the data returned by the function 
    data_ch12 = grid.ReadAsArray() - 273.15
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # Call the reprojection funcion
    grid = remap(path_ch13, variable, extent, resolution)
    # Read the data returned by the function 
    data_ch13 = grid.ReadAsArray() - 273.15
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------

    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------
    # RGB Components
    R = data_ch08 - data_ch10
    G = data_ch12 - data_ch13
    B = data_ch08

    # Minimuns and Maximuns
    Rmin = -26.2
    Rmax = 0.6
    
    Gmin = -43.2
    Gmax = 6.7

    Bmin = -29.25
    Bmax = -64.65

    R[R<Rmin] = Rmin
    R[R>Rmax] = Rmax

    G[G<Gmin] = Gmin
    G[G>Gmax] = Gmax

    B[B<Bmax] = Bmax
    B[B>Bmin] = Bmin

    # Choose the gamma
    gamma = 1

    # Normalize the data
    R = ((R - Rmin) / (Rmax - Rmin)) ** (1/gamma)
    G = ((G - Gmin) / (Gmax - Gmin)) ** (1/gamma)
    B = ((B - Bmin) / (Bmax - Bmin)) ** (1/gamma) 

    # Create the RGB
    RGB = np.stack([R, G, B], axis=2)
    #------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------

    # Formatando a descricao a ser plotada na imagem
    description = f' GOES-{satellite} Air Mass {date_img}'
    institution = "CEPAGRI - UNICAMP"

    # Definindo tamanho da imagem de saida
    d_p_i = 150
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

    # Utilizando projecao geoestacionaria no cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Adicionando o shapefile dos estados brasileiros
    adicionando_shapefile(v_extent, ax)

    # Adicionando  linhas dos litorais
    adicionando_linhas(ax)    
    
    # Defina as cores da colorbar
    colors = ['#b62007','#6f008b', '#0a0a8e', '#538234','#5C8C3A','#335a25', '#704c02', '#b57350', '#ffffff']

    # Crie uma lista de posições normalizadas para as cores
    color_positions = np.linspace(0, 1, len(colors))

    # Crie um dicionário de cores segmentadas
    cmap_dict = {'red': [], 'green': [], 'blue': [], 'alpha': []}

    for color in colors:
        rgba = to_rgba(color)  # Use a função to_rgba para converter a cor
        cmap_dict['red'].append((color_positions[colors.index(color)], rgba[0], rgba[0]))
        cmap_dict['green'].append((color_positions[colors.index(color)], rgba[1], rgba[1]))
        cmap_dict['blue'].append((color_positions[colors.index(color)], rgba[2], rgba[2]))
        cmap_dict['alpha'].append((color_positions[colors.index(color)], rgba[3], rgba[3]))

    # Crie a paleta de cores personalizada
    custom_cmap = LinearSegmentedColormap('CustomCmap', cmap_dict)
    
    # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
    img_extent = [extent[0], extent[2], extent[1], extent[3]]

    # Plot the image
    img = ax.imshow(RGB, origin='upper', cmap=custom_cmap, extent=img_extent)
        
    # Adicionando barra da paleta de cores de acordo com o canal
    cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
    cb = plt.colorbar(img, orientation="horizontal", cax=cax0, ticks=[0.2, 0.4, 0.6, 0.8])
    cb.ax.set_xticklabels(['0.2', '0.4', '0.6','0.8'])
    cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
    cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
    cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
    cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de coreses
    
    # Adicionando descricao da imagem
    adicionando_descricao_imagem(description, institution, ax, fig)

    # Adicionando os logos
    adicionando_logos(fig)

    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}{rgb_type}/{rgb_type}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)

    # Fecha a janela para limpar a memoria
    plt.close()

    # Tempo de processamento True color
    logging.info(f'Total processing time Airmass: {round((time.time() - start),2)} seconds.') 


def process_lst(file, v_extent):
    global dir_maps, dir_in, dir_out
    # Captura a hora para contagem do tempo de processamento da imagem
    start = time.time()

    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)
    satellite = '16'
    variable = 'LST'
    
    # Reprojetando
    grid = remap(file, variable, extent, resolution)

    # Lê o retorno da função
    data_lst = grid.ReadAsArray()- 273.15

    # Abrindo imagem com a biblioteca GDAL
    raw = gdal.Open(f'NETCDF:{file}:' + 'LST', gdal.GA_ReadOnly)
    metadata = raw.GetMetadata()
    dtime = metadata.get('NC_GLOBAL#time_coverage_start')
    date = (datetime.strptime(dtime, '%Y-%m-%dT%H:%M:%S.%fZ'))
    date_img = date.strftime('%d-%b-%Y %H:%M UTC')
    date_file = date.strftime('%Y%m%d_%H%M%S')

    # Define a temperatura minima
    min_temp = -26
    
    # Mask values less than -20 degrees
    data_lst = np.ma.masked_where((data_lst < min_temp), data_lst)

    # Formatando a descricao a ser plotada na imagem
    description = f' GOES-{satellite} Land Surface Temperature (°C) {date_img}'
    institution = "CEPAGRI - UNICAMP"

    # Definindo tamanho da imagem de saida
    d_p_i = 150
    fig = plt.figure(figsize=(2000 / float(d_p_i), 2000 / float(d_p_i)), frameon=True, dpi=d_p_i, edgecolor='black', facecolor='black')

    # Utilizando projecao geoestacionaria no cartopy
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Formatando a extensao da imagem, modificando ordem de minimo e maximo longitude e latitude
    img_extent = [extent[0], extent[2], extent[1], extent[3]]  # Min lon, Max lon, Min lat, Max lat
    
    # Criando imagem de fundo Natural Earth 1
    raster = gdal.Open(f'{dir_maps}HYP_HR_SR_OB_DR.tif')
    ulx, xres, xskew, uly, yskew, yres = raster.GetGeoTransform()
    lrx = ulx + (raster.RasterXSize * xres)
    lry = uly + (raster.RasterYSize * yres)
    min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]
    raster = gdal.Translate('teste.tif', raster, projWin = [min_lon, max_lat, max_lon, min_lat])
    # Lendo o RGB 
    array = raster.ReadAsArray()
    R = array[0,:,:].astype(float) / 255
    G = array[1,:,:].astype(float) / 255
    B = array[2,:,:].astype(float) / 255
    R[R==4] = 0
    G[G==5] = 0
    B[B==15] = 0
    rgb = np.stack([R, G, B], axis=2)
    
    # PLotando imagem de fundo
    ax.imshow(rgb, extent=img_extent)
    
    # Remove imagem de fundo tif criada
    os.remove('teste.tif')
    
    # Plotando a imagem Spectral_r
    img = ax.imshow(data_lst, origin='upper',vmin=-25, vmax=60, extent=img_extent, zorder=2, cmap='jet')

    # Adicionando barra da paleta de cores de acordo com o canal
    cax0 = fig.add_axes([ax.get_position().x0, ax.get_position().y0 - 0.01325, ax.get_position().width, 0.0125])
    cb = plt.colorbar(img, orientation="horizontal", cax=cax0, ticks=[-10, 10, 30, 50])
    cb.ax.set_xticklabels(['-10', '10', '30','50'])
    cb.ax.tick_params(axis='x', colors='black', labelsize=8)  # Alterando cor e tamanho dos rotulos da barra da paleta de cores
    cb.outline.set_visible(False)  # Removendo contorno da barra da paleta de cores
    cb.ax.tick_params(width=0)  # Removendo ticks da barra da paleta de cores
    cb.ax.xaxis.set_tick_params(pad=-13)  # Colocando os rotulos dentro da barra da paleta de coreses

    # Adicionando o shapefile dos estados brasileiros
    adicionando_shapefile(v_extent, ax, colors='dimgray')
    
    # Adicionando  linhas dos litorais
    adicionando_linhas(ax, colors='dimgray')

    # Adicionando descricao da imagem
    adicionando_descricao_imagem(description, institution, ax, fig)

    # Adicionando os logos
    adicionando_logos(fig)

    type_lst = 'lst'
    
    # Salvando a imagem de saida
    plt.savefig(f'{dir_out}{type_lst}/{type_lst}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)

    # Fecha a janela para limpar a memoria
    plt.close()
    
    logging.info(f'Total processing time: {round((time.time() - start),2)} seconds.')


def iniciar_processo_cmi(p_br, p_sp, bands, process_br, process_sp, new_bands):
    global dir_out
    # Checagem se e possivel gerar imagem bandas 1-16
    if p_br:
        logging.info('')
        logging.info('PROCESSANDO IMAGENS "BR"...')
        # Contador para processamento nas 16 bandas
        for x in range(1, 17):
            # Transforma o inteiro contador em string e com 2 digitos
            b = str(x).zfill(2)
            # Se a banda tiver TRUE processa
            if bands[b]:
                # Tentando Processar
                try:
                    # Cria o processo com a funcao de processamento
                    process = Process(target=process_band_cmi, args=(f'{dir_in}band{b}/{new_bands[b]}', b, 'br'))
                    # Adiciona o processo na lista de controle do processamento paralelo
                    process_br.append(process)
                    # Inicia o processo
                    process.start()
                # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
                except Exception as e:
                    logging.info(f'Erro {e} no Arquivo - {new_bands[b]}')
                    # Remove a imagem com erro de processamento
                    os.remove(f'{dir_in}band{b}/{new_bands[b]}')
            else:
                continue
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in process_br:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()
        # Limpa a lista de processos
        process_br.clear()

    if p_sp:
        logging.info('')
        logging.info('PROCESSANDO IMAGENS "SP"...')
        # Contador para processamento nas 16 bandas
        for x in range(1, 17):
            # Transforma o inteiro contador em string e com 2 digitos
            b = str(x).zfill(2)
            # Se a banda tiver novas glm_list para o dia:
            if bands[b]:
                # Tentando Processar
                try:
                    # Cria o processo com a funcao de processamento
                    process = Process(target=process_band_cmi, args=(f'{dir_in}band{b}/{new_bands[b]}', b, 'sp'))
                    # Adiciona o processo na lista de controle do processamento paralelo
                    process_sp.append(process)
                    # Inicia o processo
                    process.start()
                # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
                except Exception as e:
                    logging.info(f'Erro {e} no Arquivo - {new_bands[b]}')
                    # Remove a imagem com erro de processamento
                    os.remove(f'{dir_in}band{b}/{new_bands[b]}')
            else:
                continue
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in process_sp:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()
        # Limpa a lista de processos
        process_sp.clear()


def iniciar_processo_truecolor(p_br, p_sp, bands, new_bands):
    # Checagem se e possivel gerar imagem TrueColor
    if bands['17']:        
        # Pegando nome dos produtos 01, 02, 03, 13
        ch01 = new_bands['01']
        ch02 = new_bands['02']
        ch03 = new_bands['03']
        ch13 = new_bands['13']
        
        # Pegando nome e local e do produto
        ch01 = f'{dir_in}band01/{ch01}'
        ch02 = f'{dir_in}band02/{ch02}'
        ch03 = f'{dir_in}band03/{ch03}'
        ch13 = f'{dir_in}band13/{ch13}'
        
        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS TRUECOLOR WITH NIGHT "BR"...')
            # Tenta realizar o processamento da imagem
            try:
                # Inicia a funcao de processamento
                process_truecolor("truecolor", "br", ch01, ch02, ch03, ch13)
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")

        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_sp:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS TRUECOLOR WITH NIGHT "BR"...')
            # Tenta realizar o processamento da imagem
            try:
                # Inicia a funcao de processamento
                process_truecolor("truecolor", "sp", ch01, ch02, ch03, ch13)
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")
                
                
def iniciar_processo_rrqpef(p_br, p_sp, bands, new_bands):
    # Checagem se e possivel gerar imagem RRQPEF   
    if bands['18']:
        
        # Pega o nome netCDF da banda 13
        ch13 = new_bands['13']
        rrqpef = new_bands['18']
        rrqpef = f'{dir_in}rrqpef/{rrqpef}'

        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS RRQPEF "BR"...')

            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process_rrqpef(f'{rrqpef}', f'{dir_in}band13/{ch13.replace(".nc", "_reproj_br.nc")}', "br")
            # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
            except Exception as e:
                # Realiza o log do erro
                logging.info(f'Erro {e} Arquivo - {rrqpef}')
                # Remove a imagem com erro de processamento
                os.remove(f'{rrqpef}')
                os.remove(f'{dir_in}band13/{ch13.replace(".nc", "_reproj_br.nc")}')
            
                # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_sp:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS RRQPEF "SP"...')
            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process_rrqpef(f'{rrqpef}', f'{dir_in}band13/{ch13.replace(".nc", "_reproj_sp.nc")}', "sp")
            # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
            except:
                # Realiza o log do erro
                logging.info(f'Erro Arquivo - {rrqpef}')
                # Remove a imagem com erro de processamento
                os.remove(f'{rrqpef}')
                os.remove(f'{dir_in}band13/{ch13.replace(".nc", "_reproj_sp.nc")}')


def iniciar_processo_glm(p_br, bands, process_br, dir_in, new_bands):
    # Pega o nome netCDF da banda 13
    ch13 = new_bands['13']
    # Checagem se e possivel gerar imagem GLM
    if bands['19']:
        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS GLM "BR"...')
            # Cria uma lista com os itens presentes no diretório da banda que são arquivos e terminam com ".nc"
            glm_list = [f for f in os.listdir(f'{dir_in}glm') if os.path.isfile(os.path.join(f'{dir_in}glm', f)) and re.match('^OR_GLM-L2-LCFA_G16_s.+_e.+_c.+.nc$', f)]
            # Ordena a lista
            glm_list.sort()
            # Filtra os arq glm para pegar somente os no intervalo ini < glm < fim
            glm_list = filtrar_imagens_por_intervalo(glm_list, ch13)
            
            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process = Process(target=process_glm, args=(f'{dir_in}band13/{ch13.replace(".nc", "_reproj_br.nc")}', glm_list, "br"))
                # Adiciona o processo na lista de controle do processamento paralelo
                process_br.append(process)
                # Inicia o processo
                process.start()
            # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
            except Exception as e:
                # Realiza o log do erro
                logging.info(f'Erro {e} Arquivo ')
                # Remove a imagem com erro de processamento
                os.remove(f'{dir_in}band13/{ch13.replace(".nc", "_reproj_br.nc")}')
        
            # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
            for process in process_br:
                # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
                process.join()
            # Limpa a lista de processos
            process_br.clear()        
    else:
        logging.info(f'')
        logging.info(f'Sem imagens correspondentes a data para glm ')
        logging.info(f'')


def iniciar_processo_ndvi(p_br, bands, process_br, dir_in, new_bands):
    # Pega o nome netCDF das bandas 02, 03
    ch02 = new_bands['02']
    ch03 = new_bands['03']

    # Checagem se e possivel gerar imagem NDVI
    if bands['20']:
        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:
            logging.info('')
            logging.info('PROCESSANDO IMAGENS NDVI...')
            
            # Captura a data do arquivo
            date_file = (datetime.strptime(ch02[ch02.find("M6C02_G16_s") + 11:ch02.find("_e") - 1], '%Y%j%H%M%S'))
            # Captura a data atual
            date_now = datetime.now()
            # Aponta o horario 18h para a data atual
            date = datetime(date_now.year, date_now.month, date_now.day, int(18), int(00))
            
            # Se a data do arquivo for maior ou igual as 18h da data atual e o dia da semana atual for sabado (6)
            if date_file >= date:
                # Adiciona true para a variavel de processamento semanal
                ndvi_diario = True
            else:
                ndvi_diario = False
            
            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process = Process(target=process_ndvi, 
                                  args=(ndvi_diario, 
                                        f'{dir_in}band02/{ch02.replace(".nc", "_reproj_br.nc")}', 
                                        f'{dir_in}band03/{ch03.replace(".nc", "_reproj_br.nc")}', "br"))
                # Adiciona o processo na lista de controle do processamento paralelo
                process_br.append(process)
                # Inicia o processo
                process.start()
            
            # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
            except Exception as e:
                # Realiza o log do erro
                logging.info("Erro no processamento do Arquivo NDVI")
                logging.info(str(e))
                # Remove a imagem com erro de processamento
                os.remove(f'{dir_in}band02/{ch02.replace(".nc", "_reproj_br.nc")}')
                os.remove(f'{dir_in}band03/{ch03.replace(".nc", "_reproj_br.nc")}')
            
            # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
            for process in process_br:
                # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
                process.join()
            
            # Limpa a lista de processos
            process_br.clear()


def iniciar_processo_fdcf(p_br, bands, dir_in, new_bands):
    
    if bands['21']:
        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br: 
            # Checagem se e possivel gerar imagem FDCF
            fdcf_diario = False
            # Pegando nome das bands 01, 02, 03
            ch01 = new_bands['01']
            ch02 = new_bands['02']
            ch03 = new_bands['03']
            
            # Pega o nome do arquivo band01 para fazer comparação
            fdcf = new_bands['21']
            
            logging.info("")
            logging.info('PROCESSANDO IMAGENS FDCF "BR"...')
            
            # Captura a data do arquivo
            date_file = (datetime.strptime(fdcf[fdcf.find("ABI-L2-FDCF-M6_G16_s") + 20:fdcf.find("_e") - 1], '%Y%j%H%M%S'))
            # Captura a data atual
            date_now = datetime.now()
            # Aponta o horario 23h50 para o dia anterior                           
            date = datetime(date_now.year, date_now.month, date_now.day, int(23), int(40))
                        
            #Checagem para ver se é 23:50 para processamento do acumulado diário
            if date_file.year == date.year and date_file.month == date.month and date_file.day == date.day and date_file >= date:
                # Adiciona true para a variavel de processamento semanal
                fdcf_diario = True
                logging.info(f'fdcf_diario: {fdcf_diario}')
            else:
                # Adiciona false para a variavel de processamento semanal
                fdcf_diario = False
                logging.info(f'fdcf_diario: {fdcf_diario}')
            
            # Tenta realizar o processamento da imagem
            try:
                # Cria o processo com a funcao de processamento
                process_fdcf(f'{dir_in}fdcf/{fdcf}', f'{dir_in}band01/{ch01}', f'{dir_in}band02/{ch02}', f'{dir_in}band03/{ch03}', "br", fdcf_diario)
            # Caso seja retornado algum erro do processamento, realiza o log e remove a imagem com erro de processamento
            except Exception as e:
                # Realiza o log do erro
                logging.info("Erro Arquivo fdcf")
                logging.info(str(e))


def iniciar_processo_airmass(p_br, p_sp, bands, new_bands):
    # Checagem se e possivel gerar imagem Air Mass
    if bands['22']:
        
        logging.info("")
        logging.info('PROCESSANDO IMAGENS AIR MASS "BR"...')
        # Pegando nome das bandas 08, 10, 12, 13
        ch08 = new_bands['08']
        ch10 = new_bands['10']
        ch12 = new_bands['12']
        ch13 = new_bands['13']
        
        ch08 = f'{dir_in}band08/{ch08}'
        ch10 = f'{dir_in}band10/{ch10}'
        ch12 = f'{dir_in}band12/{ch12}'
        ch13 = f'{dir_in}band13/{ch13}'
        
        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:            
            # Tenta realizar o processamento da imagem
            try:
                # Processando imagens 
                process_airmass('airmass', 'br', ch08, ch10, ch12, ch13)
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")

        # Se a variavel de controle de processamento sp for True, realiza o processamento
        if p_sp:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS AIRMASS "SP"...')        
            # Tenta realizar o processamento da imagem
            try:
                # Processando imagens 
                process_airmass('airmass', 'sp', ch08, ch10, ch12, ch13)
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")


def iniciar_processo_lst(p_br, p_sp, bands, new_bands):
    # Checagem se e possivel gerar imagem Land Surface Temperature
    if bands['23']:
        
        # Pega o nome do produto LST2KMF 
        lst = new_bands['23']
        # Pega o local do produto 
        file = f'{dir_in}lst/{lst}'

        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if p_br:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS LAND SURFACE TEMPERATURE "BR"...')
            try:
                # Inicia o Processamento
                process_lst(file, 'br')
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")
                
        if p_sp:
            logging.info("")
            logging.info('PROCESSANDO IMAGENS LAND SURFACE TEMPERATURE "SP"...')
            try:
                # Inicia o Processamento
                process_lst(file, 'sp')
            # Caso seja retornado algum erro do processamento, realiza o log 
            except Exception as e:
                # Registra detalhes da exceção, como mensagem e tipo
                logging.error(f"Erro ao criar processo: {e}")

# ========================================#     Main     #========================================== #

def processamento_das_imagens(bands, p_br, p_sp, dir_in, dir_main): 
   
    # Cria lista vazia para controle do processamento paralelo
    process_br = []
    # Cria lista vazia para controle processamento paralelo
    process_sp = []
    
    new_bands = abrir_old_json(dir_main)
    
    try:
        iniciar_processo_cmi(p_br, p_sp, bands, process_br, process_sp, new_bands)
        
        iniciar_processo_truecolor(p_br, p_sp, bands, new_bands)
        
        iniciar_processo_rrqpef(p_br, p_sp, bands, new_bands)

        iniciar_processo_glm(p_br, bands, process_br, dir_in, new_bands)
        
        iniciar_processo_ndvi(p_br, bands, process_br, dir_in, new_bands)
        
        iniciar_processo_fdcf(p_br, bands, dir_in, new_bands)
        
        iniciar_processo_airmass(p_br, p_sp, bands, new_bands)
        
        iniciar_processo_lst(p_br, p_sp, bands, new_bands)
        
    except Exception as e:
        logging.info(f'Ocorreu um Erro {e} no Processamento')
    
    # Realiza log do encerramento do processamento
    logging.info("")
    logging.info("PROCESSAMENTO ENCERRADO")