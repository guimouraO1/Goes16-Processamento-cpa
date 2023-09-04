import os  # Importa a biblioteca os para operações de sistema.
import re  # Importa a biblioteca re para expressões regulares.
import json  # Importa a biblioteca json para manipulação de arquivos JSON.
import logging  # Importa a biblioteca logging para registrar informações.
import shutil
import datetime
from modules.utilities import download_prod

# Função para remover todos os arquivos de uma pasta, exceto um específico.
def remover_todos_exceto(nome_arquivo, pasta):
    for arquivo in os.listdir(pasta):
        caminho_arquivo = os.path.join(pasta, arquivo)
        # Verifica se o arquivo é um arquivo e se não é o arquivo especificado.
        if os.path.isfile(caminho_arquivo) and arquivo != nome_arquivo and len(arquivo) > 1:
            os.remove(caminho_arquivo)


# Função para abrir o arquivo "oldBands.json" e retornar a lista de "oldImagesName".
def abrir_old_json():
    with open('oldBands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
        return oldImages


# Função para modificar um valor em um arquivo JSON.
def modificar_chave_old_bands(caminho_arquivo, chave, novo_valor):
    with open(caminho_arquivo, 'r') as arquivo_json:
        dados = json.load(arquivo_json)
    dados['oldImagesName'][chave] = novo_valor
    with open(caminho_arquivo, 'w') as arquivo_json:
        json.dump(dados, arquivo_json, indent=4)


# Checa bandas 1-16
def checar_bandas(bands, dir_in):
    # Checagem imagens ABI 1-16
    for x in range(1, 17):
        # Formata um int para um string de dois dígitos (01, 02, ..., 16).
        b = str(x).zfill(2)
        # Obtém uma lista de imagens que correspondem a um padrão específico na pasta.
        imagens = [f for f in os.listdir(f'{dir_in}band{b}') if os.path.isfile(os.path.join(f'{dir_in}band{b}', f)) and re.match('^CG_ABI-L2-CMIPF-M[0-9]C[0-1][0-9]_G16_s.+_e.+_c.+[0-9].nc$', f)]
        # Se houver imagens na pasta:
        if imagens:
            # Encontra a imagem mais recente na lista.
            latestBand = max(imagens)
            # Obtém as imagens antigas.
            old_bands = abrir_old_json()
            
            # Se houver uma imagem mais recente e ela for diferente das antigas:
            if latestBand and latestBand != old_bands[b]: 
                logging.info(f'Novas imagens para o dia band{b}')
                
                # Se houver mais de uma imagem na pasta:
                if len(imagens) > 1:
                    # Remove os arquivos netCDF menos o mais atual
                    remover_todos_exceto(latestBand, f'{dir_in}band{b}/') 
                
                # Modifica o arquivo JSON com a imagem mais recente.               
                modificar_chave_old_bands('oldBands.json', b, latestBand)  
                # Atualiza o dicionário "bands" com true para novas imagens.
                bands[b] = True
            else:
                logging.info(f'Sem imagens para o dia band{b}')
                # Atualiza o dicionário "bands" com false sem novas imagens.
                bands[b] = False
        else:
            logging.info(f'Sem imagens para o dia band{b}') 
            # Atualiza o dicionário "bands" com false sem novas imagens.
            bands[b] = False


# Checa se há bandas 1, 2, 3 para true color
def checar_truecolor(bands):
        # Checagem de novas imagens truecolor (Band 17) se todas as bands 1, 2, 3 forem True
    if all(bands[str(x).zfill(2)] for x in range(1, 4)):
        # Se Todas as três bandas são True
        bands['17'] = True
        logging.info(f'Novas imagens TRUECOLOR')
    else:
        bands["17"] = False
        logging.info(f'Sem novas imagens TRUECOLOR')


# Checa se há bandas 13 para rrqpef e baixa o produto
def checar_rrqpef(bands, dir_in):
    old_bands = abrir_old_json()
    # Checagem de novas imagens rrqpef (Band 18)
    if bands['13']:
        ch13 = old_bands['13']
        # Extrair o data/hora do arquivo da banda 13 para download do arquivo RRQPEF
        ftime = (datetime.datetime.strptime(ch13[ch13.find("M6C13_G16_s") + 11:ch13.find("_e") - 1], '%Y%j%H%M%S'))
        try:
            # Download arquivo rrqpef
            download_prod(datetime.datetime.strftime(ftime, '%Y%m%d%H%M'), 'ABI-L2-RRQPEF', f'{dir_in}rrqpef/')
            bands['18'] = True
            logging.info(f'Novas imagens RRQPEF')
        except:
            logging.info(f'Erro ao baixar novas imagens RRQPEF')
            bands['18'] = False
    else:
        logging.info(f'Sem novas imagens RRQPEF')
        bands['18'] = False


# Checa se há bandas 13 para glm
def checar_glm(bands, dir_in):
    # Checagem de novas imagens GLM (Band 19)
    if bands['13']:
        # Pega lista de glm para verificação
        glm_list = [f for f in os.listdir(f'{dir_in}glm') if os.path.isfile(os.path.join(f'{dir_in}glm', f)) and re.match('^OR_GLM-L2-LCFA_G16_s.+_e.+_c.+.nc$', f)]
        glm_list.sort()
        # Se a lista for maior que 0 True
        if len(glm_list) > 0:
            bands['19'] = True
            logging.info('Novas imagens GLM')
        else:
            bands['19'] = False
            logging.info('Sem novas imagens GLM')
    else:
        bands['19'] = False
        logging.info('Sem novas imagens GLM')


# Checa se há bandas 2,3 para ndvi
def checar_ndvi(bands, dir_in):
    
    old_bands = abrir_old_json()
    
    # Checagem de novas imagens ndvi (Band 20)
    if bands['02'] and bands['03']:
        # Carrega os arquivos de processamento das bandas para composicao do ndvi
        file_ch02 = old_bands['02']
        
        # Listar arquivos no diretório da pasta "band03"
        file_ch03_dir = f'{dir_in}band03/'
        file_ch03_list = os.listdir(file_ch03_dir)

        date_now = datetime.datetime.now()
        date_ini = datetime.datetime(date_now.year, date_now.month, date_now.day, 13, 0)
        date_end = date_ini + datetime.timedelta(hours=5, minutes=1)

        date_file = datetime.datetime.strptime(file_ch02[file_ch02.find("M6C02_G16_s") + 11:file_ch02.find("_e") - 1], '%Y%j%H%M%S')

        if date_ini <= date_file <= date_end:
            # Verifica se há arquivo correspondente na banda 03
            matches_ch03 = [z for z in file_ch03_list if z.startswith(file_ch02[0:43].replace('M6C02', 'M6C03'))]
            
            if file_ch02 and matches_ch03:
                bands['20'] = True
                logging.info(f'Novas imagens NDVI')
            else:
                bands['20'] = False
                logging.info(f'Sem novas imagens NDVI')
        else:
            bands['20'] = False
            logging.info(f'Sem novas imagens NDVI')


def checar_fdcf(bands, dir_in):  
    
    # Checagem de novas imagens fdcf (Band 21)
    if bands['17']:
        # Coleto o nome das novas bandas
        old_bands = abrir_old_json()
        # Extrai a data/hora do arquivo para download do arquivo FDCF
        ch01 = old_bands['01']
        # Extrai a data/hora do arquivo ch01
        ftime = (datetime.datetime.strptime(ch01[ch01.find("M6C01_G16_s") + 11:ch01.find("_e") - 1], '%Y%j%H%M%S'))
        try:
            # Download arquivo fdcf
            name_fdcf = download_prod(datetime.datetime.strftime(ftime, '%Y%m%d%H%M'), "ABI-L2-FDCF", f'{dir_in}fdcf/')
            print(name_fdcf)
            
            # Modifica o arquivo JSON com a imagem mais recente.               
            modificar_chave_old_bands(f'oldBands.json', '21', name_fdcf)
            bands['21'] = True
            logging.info(f'novas imagens FDCF')
        except:
            bands['21'] = False
            logging.info(f'Sem novas imagens FDCF')
    else:
        bands['21'] = False
        logging.info(f'Sem novas imagens FDCF')


# ========================================#     Main     #========================================== #

# Função para verificar a existência de novas imagens
def checar_imagens(bands, dir_in):
    
    logging.info("VERIFICANDO NOVAS IMAGENS")

    checar_bandas(bands, dir_in)

    checar_truecolor(bands)
    
    checar_rrqpef(bands, dir_in)
    
    checar_glm(bands, dir_in)   
    
    checar_ndvi(bands, dir_in)

    checar_fdcf(bands, dir_in)
    
    print(bands)
    
    # Retorna o dicionário "bands".      
    return bands