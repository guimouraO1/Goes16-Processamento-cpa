import os  # Importa a biblioteca os para operações de sistema.
import re  # Importa a biblioteca re para expressões regulares.
import json  # Importa a biblioteca json para manipulação de arquivos JSON.
import logging  # Importa a biblioteca logging para registrar informações.
import datetime
from modules.utilities import download_prod


# ============================================# Diretórios ========================================= #
# Função para abrir o arquivo "oldBands.json" e retornar a lista de "oldImagesName".
def abrir_old_json(dir_main):
    with open(f'{dir_main}old_bands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
        return oldImages


# filtra arquivos glm para verificar se correspondem a data
def filtrar_imagens_por_intervalo(images, ch13):
    # Extrai a data e hora da string 'ch13' e define um intervalo de 9 minutos e 59 segundos a partir dela.
    glm_list = [] 
    ch13_data = (datetime.datetime.strptime(ch13[ch13.find("M6C13_G16_s") + 11:ch13.find("_e") - 1], '%Y%j%H%M%S'))
    date_ini = datetime.datetime(ch13_data.year, ch13_data.month, ch13_data.day, ch13_data.hour, ch13_data.minute)
    date_end = datetime.datetime(ch13_data.year, ch13_data.month, ch13_data.day, ch13_data.hour, ch13_data.minute) + datetime.timedelta(minutes=9, seconds=59)
    # Percorre a lista de nomes de imagens e verifica se a data e hora de cada imagem estão dentro do intervalo.
    for x in images:
        xtime = (datetime.datetime.strptime(x[x.find("GLM-L2-LCFA_G16_s") + 17:x.find("_e") - 1], '%Y%j%H%M%S'))
        if date_ini <= xtime <= date_end:
            glm_list.append(x)
        else:
            continue
    
    return glm_list


# Função para remover todos os arquivos de uma pasta, exceto um específico.
def remover_todos_exceto(nome_arquivo, pasta):
    for arquivo in os.listdir(pasta):
        caminho_arquivo = os.path.join(pasta, arquivo)
        # Verifica se o arquivo é o arquivo novo, se não for apaga(Serve para apagar arquivos duplicados ou pasta goes/ cheia de arquivos nc não processados)
        if os.path.isfile(caminho_arquivo) and arquivo != nome_arquivo and len(arquivo) > 1:
            os.remove(caminho_arquivo)


# Função para modificar um valor em um arquivo JSON.
def modificar_chave_old_bands(caminho_arquivo, chave, novo_valor):
    with open(caminho_arquivo, 'r') as arquivo_json:
        dados = json.load(arquivo_json)
    dados['oldImagesName'][chave] = novo_valor
    with open(caminho_arquivo, 'w') as arquivo_json:
        json.dump(dados, arquivo_json, indent=4)


# Checa bandas 1-16
def checar_bandas(bands, dir_in, dir_main):
    
    # Obtém as imagens antigas.
    old_bands = abrir_old_json(dir_main)
    # Checagem imagens ABI 1-16
    for x in range(1, 17):
        # Formata um int para um string de dois dígitos (01, 02, ..., 16).
        b = str(x).zfill(2)
        # Obtém uma lista de imagens que correspondem a um padrão específico na pasta.
        imagens = [f for f in os.listdir(f'{dir_in}band{b}') if os.path.isfile(os.path.join(f'{dir_in}band{b}', f)) and re.match('^CG_ABI-L2-CMIPF-M[0-9]C[0-1][0-9]_G16_s.+_e.+_c.+[0-9].nc$', f)]
        # Se houver imagens na pasta e a imagem mais recente for diferente da ultimo processamento
        if imagens and max(imagens) != old_bands[b]:
            logging.info(f'Novas imagens para o dia band{b}')
            # latestBand recebe a nova imagem
            latestBand = max(imagens)
            # Se houver mais de uma imagem na pasta:
            if len(imagens) > 1:
                # Remove os arquivos netCDF menos o mais atual
                remover_todos_exceto(latestBand, f'{dir_in}band{b}/') 
            # Modifica o arquivo JSON com a imagem mais recente.               
            modificar_chave_old_bands('old_bands.json', b, latestBand)  
            # Atualiza o dicionário "bands" com true para novas imagens.
            bands[b] = True
        else:
            logging.info(f'Sem imagens para o dia band{b}') 
            # Atualiza o dicionário "bands" com false sem novas imagens.
            bands[b] = False


# Checa se há bandas 1, 2, 3 para true color
def checar_truecolor(bands):
    # Checagem de novas imagens truecolor (Band 17) se todas as bands 1, 2, 3 forem True
    if all(bands[str(x).zfill(2)] for x in range(1, 4)) and bands['13'] == True:
        # Se Todas as três bandas são True
        bands['17'] = True
        logging.info(f'Novas imagens TRUECOLOR')
    else:
        bands["17"] = False
        logging.info(f'Sem novas imagens TRUECOLOR')


# Checa se há banda 13 para rrqpef e baixa o produto
def checar_rrqpef(bands, dir_in, old_bands):
    # Checagem de novas imagens rrqpef (Band 18)
    if bands['13']:
        ch13 = old_bands['13']
        # Extrair o data/hora do arquivo da banda 13 para download do arquivo RRQPEF
        ftime = (datetime.datetime.strptime(ch13[ch13.find("M6C13_G16_s") + 11:ch13.find("_e") - 1], '%Y%j%H%M%S'))
        try:
            # Download arquivo rrqpef
            rrqpef = download_prod(datetime.datetime.strftime(ftime, '%Y%m%d%H%M'), 'ABI-L2-RRQPEF', f'{dir_in}rrqpef/')
            modificar_chave_old_bands(f'old_bands.json', '18', f'{rrqpef}.nc')
            bands['18'] = True
            logging.info(f'Novas imagens RRQPEF')
        except:
            logging.info(f'Erro ao baixar novas imagens RRQPEF')
            bands['18'] = False
    else:
        logging.info(f'Sem novas imagens RRQPEF')
        bands['18'] = False


# Checa se há bandas 13 para glm
def checar_glm(bands, dir_in, old_bands):
    # Checagem de novas imagens GLM (Band 19)
    if bands['13']:
            ch13 = old_bands['13'] 
            # Cria uma lista com os itens presentes no diretório da banda que são arquivos e terminam com ".nc"
            glm_list = [f for f in os.listdir(f'{dir_in}glm') if os.path.isfile(os.path.join(f'{dir_in}glm', f)) and re.match('^OR_GLM-L2-LCFA_G16_s.+_e.+_c.+.nc$', f)]
            # Ordena a lista
            glm_list.sort()
            # Filtra os arq glm para pegar somente os no intervalo ini < glm < fim
            glm_list = filtrar_imagens_por_intervalo(glm_list, ch13)
            
            if len(glm_list) > 0:
                # Tenta realizar o processamento da imagem
                bands['19'] = True
                logging.info('Novas imagens GLM')
            else:
                bands['19'] = False
                logging.info('Sem novas imagens GLM')
    else:
        bands['19'] = False
        logging.info('Sem novas imagens GLM')


# Checa se há bandas 2, 3 para ndvi band 20
def checar_ndvi(bands, dir_in, old_bands):
    # Checagem de novas imagens ndvi (Band 20)
    if bands['02'] and bands['03']:
        # Carrega os arquivos de processamento das bandas para composicao do ndvi
        file_ch02 = old_bands['02']
        # Converte o nome do arquivo da banda 02 para o formato da banda 03
        file_ch03 = file_ch02[0:43].replace('M6C02', 'M6C03')
        
        # Listar arquivos no diretório da pasta "band03"
        file_ch03_dir = f'{dir_in}band03/'
        file_ch03_list = os.listdir(file_ch03_dir)
        
        # Verifica date now, nosso dia colocando horario 13h utc, e end 18h utc
        date_now = datetime.datetime.now()
        date_ini = datetime.datetime(date_now.year, date_now.month, date_now.day, int(13), int(00))
        date_end = datetime.datetime(date_now.year, date_now.month, date_now.day, int(18), int(1))
        
        # Data do file band02 para fazer a verificação da hora
        date_file = datetime.datetime.strptime(file_ch02[file_ch02.find("M6C02_G16_s") + 11:file_ch02.find("_e") - 1], '%Y%j%H%M%S')
        
        # Processa apenas se o arquivo estiver entre 13h utc as 18h utc
        if date_ini <= date_file <= date_end:
            
            # Verifica se há pelo menos um arquivo correspondente na banda 03
            found_match = any(file_ch03_candidate.startswith(file_ch03) for file_ch03_candidate in file_ch03_list)
            
            # Se encontrar um arquivo correspondente
            if found_match:
                bands['20'] = True
                logging.info(f'Novas imagens NDVI')
            else:
                bands['20'] = False
                logging.info(f'Sem novas imagens NDVI')
        else:
            bands['20'] = False
            logging.info(f'Sem novas imagens NDVI')


# Checa se há bandas 1, 2, 3 (truecolor) para fdcf band 21
def checar_fdcf(bands, dir_in, old_bands):  
    # Checagem de novas imagens fdcf (Band 21)
    if bands['17']:
        # Extrai a data/hora do arquivo para download do arquivo FDCF
        ch01 = old_bands['01']
        # Extrai a data/hora do arquivo ch01
        ftime = (datetime.datetime.strptime(ch01[ch01.find("M6C01_G16_s") + 11:ch01.find("_e") - 1], '%Y%j%H%M%S'))
        try:
            # Download arquivo fdcf
            fdcf = download_prod(datetime.datetime.strftime(ftime, '%Y%m%d%H%M'), "ABI-L2-FDCF", f'{dir_in}fdcf/')
            modificar_chave_old_bands(f'old_bands.json', '21', f'{fdcf}.nc')
            bands['21'] = True
            logging.info(f'Novas imagens FDCF')
        except:
            bands['21'] = False
            logging.info(f'Falha ao baixar FDCF')
    else:
        bands['21'] = False
        logging.info(f'Sem novas imagens FDCF')
        
        
# Checa se há bandas  bandas 08, 10, 12, 13 (truecolor) para Airmass
def checar_airmass(bands):  
    # Checagem de novas imagens fdcf (Band 22)
    if bands['08'] and bands['10'] and bands['12'] and bands['13']:
        bands['22'] = True
        logging.info(f'Novas imagens AIRMASS')
    else:
        bands['22'] = False
        logging.info(f'Sem novas imagens AIRMASS')
        
        
# Checa se há o Produto Land Sarface Temperature
def checar_lst(bands, dir_in, old_bands):

    # Obtém uma lista de imagens que correspondem a um padrão específico na pasta.
    imagens = [f for f in os.listdir(f'{dir_in}lst/') if os.path.isfile(os.path.join(f'{dir_in}lst/', f)) and re.match('^OR_ABI-L2-LST2KMF-M[0-9]_G16_s.+_e.+_c.+[0-9].nc$', f)]
    # Se houver imagens na pasta e a imagem mais recente for diferente da ultimo processamento                           
    if imagens and max(imagens) != old_bands['23']:
        logging.info(f'Novas imagens para Land Surface Temperature')
        # latestBand recebe a nova imagem
        latestBand = max(imagens)
        # Se houver mais de uma imagem na pasta:
        if len(imagens) > 1:
            # Remove os arquivos netCDF menos o mais atual
            remover_todos_exceto(latestBand, f'{dir_in}lst/') 
        # Modifica o arquivo JSON com a imagem mais recente.               
        modificar_chave_old_bands('old_bands.json', '23', latestBand)  
        # Atualiza o dicionário "bands" com true para novas imagens.
        bands['23'] = True
    else:
        logging.info(f'Sem imagens para Land Surface Temperature') 
        # Atualiza o dicionário "bands" com false sem novas imagens.
        bands['23'] = False
    

# ========================================#     Main     #========================================== #

# Função para verificar a existência de novas imagens
def checar_imagens(bands, dir_in, dir_main):
    
    logging.info("VERIFICANDO NOVAS IMAGENS")

    checar_bandas(bands, dir_in, dir_main)

    new_bands = abrir_old_json(dir_main)    
    
    checar_truecolor(bands)
    
    checar_rrqpef(bands, dir_in, new_bands)
    
    checar_glm(bands, dir_in, new_bands)   
    
    checar_ndvi(bands, dir_in, new_bands)

    checar_fdcf(bands, dir_in, new_bands)
    
    checar_airmass(bands)
    
    checar_lst(bands, dir_in, new_bands)
    
    print(bands)

    return bands

# ========================================#     Main     #========================================== #