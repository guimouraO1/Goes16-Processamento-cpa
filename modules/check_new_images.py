import os  # Importa a biblioteca os para operações de sistema.
import re  # Importa a biblioteca re para expressões regulares.
import json  # Importa a biblioteca json para manipulação de arquivos JSON.
import logging  # Importa a biblioteca logging para registrar informações.
import shutil
import datetime
from libs.utilities import download_prod

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


# Função para verificar a existência de novas imagens.
def checar_imagens(bands, dir_in):
    logging.info("VERIFICANDO NOVAS IMAGENS")
    # Chegagem imagens ABI 1-16
    for x in range(1, 17):
        # Formata um int para um string de dois dígitos (01, 02, ..., 16).
        b = str(x).zfill(2)
        # Obtém uma lista de imagens que correspondem a um padrão específico na pasta.
        imagens = [f for f in os.listdir(f'{dir_in}band{b}') if os.path.isfile(os.path.join(f'{dir_in}band{b}', f)) and re.match('^CG_ABI-L2-CMIPF-M[0-9]C[0-1][0-9]_G16_s.+_e.+_c.+.nc$', f)]
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
            

    # Checagem de novas imagens truecolor (Band 17) if bands 1, 2, 3
    if all(bands[str(x).zfill(2)] for x in range(1, 4)):
        # Se Todas as três bandas são True
        bands['17'] = True
        logging.info(f'Novas imagens TRUECOLOR')
    else:
        bands["17"] = False
        logging.info(f'Sem novas imagens TRUECOLOR')


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


    # Checagem de novas imagens GLM (Band 19)
    if bands['13']:
        # Pega lista de glm para verificação
        glm_list = [f for f in os.listdir(f'{dir_in}glm') if os.path.isfile(os.path.join(f'{dir_in}glm', f)) and re.match('^OR_GLM-L2-LCFA_G16_s.+_e.+_c.+.nc$', f)]
        glm_list.sort()
        
        # Se não tem arquivos glm para processar
        if not glm_list:
            bands['19'] = False
            logging.info('Sem novas imagens GLM')
        else:
            bands['19'] = True
            logging.info('Novas imagens GLM')
    else:
        bands['19'] = False
        logging.info('Sem novas imagens GLM')


    # Retorna o dicionário "bands".      
    return bands  