import logging
from string import ascii_letters, digits  # Utilitario para trabalhar com ascii
import os
from difflib import Differ  # Utilitario para verificar as diferencas entre dois arquivos
import re  # Utilitario para trabalhar com expressoes regulares
import datetime
from modules.utilities import download_prod
from modules.dirs import get_dirs


# ============================================# Diretórios ========================================= #
dirs = get_dirs()
dir_temp = dirs['dir_temp']
# ============================================# Diretórios ========================================= #


# ============================================#   Funções  ============================================== #

# Le o arquivo de processamento e retorna a lista
def read_process_file(banda):
    with open(f'{dir_temp}{banda}_process.txt', 'r') as fo:
        return fo.readlines()

# Ordena a lista em ordem alfabética
def alphanumeric_key(text):
    """Return a key based on letters and digits in `text`."""
    return [c.lower() for c in text if c in ascii_letters + digits]

# Ordena lista e cria novo arquivo new.txt em temp/
def write_new_file(banda, file):
    # Ordena de forma alfabetica a lista
    file.sort(key=alphanumeric_key)
    # Cria o arquivo _new.txt com as novas imagens que estao na lista
    with open(f'{dir_temp}{banda}_new.txt', 'w') as fo:
        fo.writelines(map(lambda f: f + '\n', file))

def write_process_file(banda):
    # Cria o arquivo band??_old.txt se nao existe
    if not os.path.isfile(f'{dir_temp}{banda}_old.txt'):
        with open(f'{dir_temp}{banda}_old.txt', 'w') as fo:
            fo.close()

    # Compara os arquivos band??_old.txt e band??_new.txt
    with open(f'{dir_temp}{banda}_old.txt', 'r') as old, open(f'{dir_temp}{banda}_new.txt', 'r') as new:
        differ = Differ()
        # Realiza a comparacao entre os arquivos e cria uma lista de imagens que estao unicamente no arquivo band??_new.txt
        process_list = [line.strip()[2::] for line in differ.compare(old.readlines(), new.readlines()) if line.startswith('+')]
        print(process_list)

    # Cria o arquivo band??_process.txt com as imagens para processamento
    with open(f'{dir_temp}{banda}_process.txt', 'w') as process:
        # Escreve as imagens da lista no arquivo, sendo cada uma em uma linha
        process.writelines(map(lambda f: f + '\n', process_list))

    # Se houveram imagens para processamento, retorna True, caso contrario retorna False
    if process_list:
        return True
    else:
        return False
# ============================================# Funções ============================================== #
# Checagem de imagens novas
def check_images(c_bands, dir_in, dir_temp):
    logging.info("VERIFICANDO NOVAS IMAGENS")

# ============================================# bands 1-16 ============================================== #
    # Contado para checagem de novas imagens nas 16 bandas
    for x in range(1, 17):
        # Transforma o inteiro co   ntador em string e com 2 digitos
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
# ============================================# bands 1-16 ============================================== #
    

# ============================================# truecolor ============================================== #
    # Checagem de novas imagens truecolor (Band 17)
    if c_bands["01"] and c_bands["02"] and c_bands["03"]:
        # Carrega os arquivos de processamento das bandas para composicao do truecolor
        file_ch01 = read_process_file('band01')
        file_ch02 = read_process_file('band02')
        file_ch03 = read_process_file('band03')
        # Cria lista vazia para adicionar os produtos que forem baixados
        product_list = []
        # Looping a partir dos arquivos da banda 01 que compoem o truecolor
        for x in file_ch01:
            # Verifica se ha arquivo correspondente na banda 02
            matches_ch02 = [y for y in file_ch02 if y.startswith(x[0:43].replace('M6C01', 'M6C02'))]
            # Verifica se ha arquivo correspondente na banda 03
            matches_ch03 = [z for z in file_ch03 if z.startswith(x[0:43].replace('M6C01', 'M6C03'))]
            # Se houver arquivos de mesma data nas 3 bandas
            if x and matches_ch02 and matches_ch03:
                product_list.append(f'{x.strip()};{matches_ch02[0].strip()};{matches_ch03[0].strip()}')
        if product_list:
            # Cria o arquivo com a nova lista de imagens
            write_new_file("band17", product_list)
            # Realiza a comparacao da lista nova com a lista antiga, se houver novas imagens, cria o arquivo de processamento e aponta True no dicionario de controle das bandas
            c_bands["17"] = write_process_file("band17")
            if c_bands["17"]:
                logging.info(f'Novas imagens TRUECOLOR')
            else:
                logging.info(f'Sem novas imagens TRUECOLOR')
                os.remove(f'{dir_temp}band17_process.txt')
            # Transforma o arquivo band??_new.txt para band??_old.txt
            os.replace(f'{dir_temp}band17_new.txt', f'{dir_temp}band17_old.txt')
        else:
            c_bands["17"] = False
            logging.info(f'Sem novas imagens TRUECOLOR')
    else:
        logging.info(f'Sem novas imagens TRUECOLOR')
# ============================================# truecolor ============================================== #


# ============================================# rrqpef ============================================== #
    # Checagem de novas imagens rrqpef (Band 18)
    if c_bands["13"]:
        # Carrega a banda 13 que sera utilizada para compor o fundo
        base = read_process_file('band13')
        # Cria lista vazia para adicionar os produtos que forem baixados
        product_list = []
        # Looping para fazer o download de cada arquivo RRQPEF correspondente ao arquivo da banda 13 existente
        for f in base:
            # Extrair o data/hora do arquivo da banda 13 para download do arquivo RRQPEF
            ftime = (datetime.datetime.strptime(f[f.find("M6C13_G16_s") + 11:f.find("_e") - 1], '%Y%j%H%M%S'))
            try:
                # Download arquivo rrqpef
                file_name = download_prod(datetime.datetime.strftime(ftime, '%Y%m%d%H%M'), 'ABI-L2-RRQPEF', f'{dir_in}rrqpef/')
            except:
                continue
            # Adicona o nome do arquivo na lista
            product_list.append(f'{file_name}.nc;{f.strip()}')
        if product_list:
            # Cria o arquivo com a nova lista de imagens
            write_new_file("band18", product_list)
            # Realiza a comparacao da lista nova com a lista antiga, se houver novas imagens, cria o arquivo de processamento e aponta True no dicionario de controle das bandas
            c_bands["18"] = write_process_file("band18")
            if c_bands["18"]:
                logging.info(f'Novas imagens RRQPEF')
            else:
                logging.info(f'Sem novas imagens RRQPEF')
                os.remove(f'{dir_temp}band18_process.txt')
            # Transforma o arquivo band??_new.txt para band??_old.txt
            os.replace(f'{dir_temp}band18_new.txt', f'{dir_temp}band18_old.txt')
        else:
            c_bands["18"] = False
            logging.info(f'Sem novas imagens RRQPEF')
    else:
        logging.info(f'Sem novas imagens RRQPEF')
# ============================================# rrqpef ============================================== #


# ============================================# GLM ============================================== #
    # Checagem de novas imagens GLM (Band 19)
    if c_bands["13"]:
        # Carrega a banda 13 que sera utilizada para compor o fundo
        base = read_process_file('band13')
        # Cria uma lista com os itens presentes no diretorio da banda que sao arquivo e terminam com ".nc"
        imagens = [f for f in os.listdir(f'{dir_in}glm') if os.path.isfile(os.path.join(f'{dir_in}glm', f)) and re.match('^OR_GLM-L2-LCFA_G16_s.+_e.+_c.+.nc$', f)]
        imagens.sort()
        # Cria lista vazia para adicionar os produtos
        product_list = []
        aux_list = []
        # Looping a partir dos arquivos da banda 13 que compoem o fundo
        for f in base:
            # Extrair o data/hora do arquivo da banda 13 para download do arquivo GLM
            ftime = (datetime.datetime.strptime(f[f.find("M6C13_G16_s") + 11:f.find("_e") - 1], '%Y%j%H%M%S'))
            date_ini = datetime.datetime(ftime.year, ftime.month, ftime.day, ftime.hour, ftime.minute)
            date_end = datetime.datetime(ftime.year, ftime.month, ftime.day, ftime.hour, ftime.minute) + datetime.timedelta(minutes=9, seconds=59)
            for x in imagens:
                xtime = (datetime.datetime.strptime(x[x.find("GLM-L2-LCFA_G16_s") + 17:x.find("_e") - 1], '%Y%j%H%M%S'))
                if date_ini <= xtime <= date_end:
                    aux_list.append(x.strip())
                else:
                    continue
            if aux_list:
                product_list.append(f'{f.strip()};{aux_list}')
            else:
                continue
        # Se houver arquivos na lista, realiza o organizacao dos arquivos de controle e processamento, caso contrario aponta False no dicionario de controle das bandas
        if product_list:
            # Cria o arquivo com a nova lista de imagens
            write_new_file("band19", product_list)
            # Realiza a comparacao da lista nova com a lista antiga, se houver novas imagens, cria o arquivo de processamento e aponta True no dicionario de controle das bandas
            c_bands["19"] = write_process_file("band19")
            if c_bands["19"]:
                logging.info(f'Novas imagens GLM')
            else:
                logging.info(f'Sem novas imagens GLM')
                os.remove(f'{dir_temp}band19_process.txt')
            # Transforma o arquivo band??_new.txt para band??_old.txt
            os.replace(f'{dir_temp}band19_new.txt', f'{dir_temp}band19_old.txt')
        else:
            c_bands["19"] = False
            logging.info(f'Sem novas imagens GLM')
    else:
        logging.info(f'Sem novas imagens GLM')
# ============================================# GLM ============================================== #


# ============================================# NDVI ============================================== #
    # Checagem de novas imagens ndvi (Band 20)
    if c_bands["02"] and c_bands["03"]:
        # Carrega os arquivos de processamento das bandas para composicao do ndvi
        file_ch02 = read_process_file('band02')
        file_ch03 = read_process_file('band03')
        # Cria lista vazia para adicionar os produtos que forem baixados
        product_list = []
        # Looping a partir dos arquivos da banda 02 que compoem o ndvi
        for x in file_ch02:
            date_now = datetime.datetime.now()
            date_ini = datetime.datetime(date_now.year, date_now.month, date_now.day, int(13), int(00))
            date_end = datetime.datetime(date_now.year, date_now.month, date_now.day, int(13), int(00)) + datetime.timedelta(hours=5, minutes=1)
            date_file = (datetime.datetime.strptime(x[x.find("M6C02_G16_s") + 11:x.find("_e") - 1], '%Y%j%H%M%S'))
            if date_ini <= date_file <= date_end:
                # Verifica se ha arquivo correspondente na banda 03
                matches_ch03 = [z for z in file_ch03 if z.startswith(x[0:43].replace('M6C02', 'M6C03'))]
                # Se houver arquivos de mesma data nas 2 bandas
                if x and matches_ch03:
                    product_list.append(f'{x.strip()};{matches_ch03[0].strip()}')
            else:
                continue
        if product_list:
            # Cria o arquivo com a nova lista de imagens
            write_new_file("band20", product_list)
            # Realiza a comparacao da lista nova com a lista antiga, se houver novas imagens, cria o arquivo de processamento e aponta True no dicionario de controle das bandas
            c_bands["20"] = write_process_file("band20")
            if c_bands["20"]:
                logging.info(f'Novas imagens NDVI')
            else:
                logging.info(f'Sem novas imagens NDVI')
                os.remove(f'{dir_temp}band20_process.txt')
            # Transforma o arquivo band??_new.txt para band??_old.txt
            os.replace(f'{dir_temp}band20_new.txt', f'{dir_temp}band20_old.txt')
        else:
            c_bands["20"] = False
            logging.info(f'Sem novas imagens NDVI')
    else:
        logging.info(f'Sem novas imagens NDVI')

# ============================================# NDVI perguntar para o joao ============================================== #
    # Checagem de novas imagens fdcf (Band 21)
    if c_bands["17"]:
        # Carrega os arquivos de processamento do truecolor
        file_truecolor = read_process_file('band17')
        # Cria lista vazia para adicionar os produtos que forem baixados
        product_list = []
        for i in file_truecolor:
            # Remove possiveis espacos vazios no inicio ou final da string e separa cada termo como um elemento
            i = i.strip().split(';')
            # Montando dicionario de argumentos
            kwargs = {'ch01': f'{dir_in}band01/{i[0].replace(".nc", "_reproj_br.nc")}', 'ch02': f'{dir_in}band02/{i[1].replace(".nc", "_reproj_br.nc")}',
                    'ch03': f'{dir_in}band03/{i[2].replace(".nc", "_reproj_br.nc")}'}
            # Extrair o data/hora do arquivo para download do arquivo FDCF
            f = kwargs["ch01"]
            ftime = (datetime.datetime.strptime(f[f.find("M6C01_G16_s") + 11:f.find("_e") - 1], '%Y%j%H%M%S'))
            try:
                # Download arquivo fdcf
                file_name = download_prod(datetime.datetime.strftime(ftime, '%Y%m%d%H%M'), "ABI-L2-FDCF", f'{dir_in}fdcf/')
            except:
                continue
            # Adicona o nome do arquivo na lista
            product_list.append(f'{dir_in}fdcf/{file_name}.nc;{kwargs["ch01"]};{kwargs["ch02"]};{kwargs["ch03"]}')

        if product_list:
            # Cria o arquivo com a nova lista de imagens
            write_new_file("band21", product_list)
            # Realiza a comparacao da lista nova com a lista antiga, se houver novas imagens, cria o arquivo de processamento e aponta True no dicionario de controle das bandas
            c_bands["21"] = write_process_file("band21")
            if c_bands["21"]:
                logging.info(f'Novas imagens FDCF')
            else:
                logging.info(f'Sem novas imagens FDCF')
                os.remove(f'{dir_temp}band21_process.txt')
            # Transforma o arquivo band??_new.txt para band??_old.txt
            os.replace(f'{dir_temp}band21_new.txt', f'{dir_temp}band21_old.txt')
        else:
            c_bands["21"] = False
            logging.info(f'Sem novas imagens FDCF')
    else:
        logging.info(f'Sem novas imagens FDCF')

    # Retorna o dicionario de controle das bandas
    return c_bands
# ============================================# NDVI perguntar para o joao ============================================== #