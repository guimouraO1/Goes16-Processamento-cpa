import logging
from string import ascii_letters, digits  # Utilitario para trabalhar com ascii
import os
from difflib import Differ  # Utilitario para verificar as diferencas entre dois arquivos
import re  # Utilitario para trabalhar com expressoes regulares
import datetime
from modules.utilities import download_prod
from modules.dirs import get_dirs

# diretórios utilizados
dirs = get_dirs()
dir_temp = dirs['dir_temp']

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

# Compara os arquivos band??_old.txt e band??_new.txt
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

    # Cria o arquivo band??_process.txt com as imagens para processamento
    with open(f'{dir_temp}{banda}_process.txt', 'w') as process:
        # Escreve as imagens da lista no arquivo, sendo cada uma em uma linha
        process.writelines(map(lambda f: f + '\n', process_list))

    # Se houveram imagens para processamento, retorna True, caso contrario retorna False
    if process_list:
        return True
    else:
        return False

# Checagem de imagens novas
def check_images(c_bands, dir_in, dir_temp):
    logging.info("VERIFICANDO NOVAS IMAGENS")

    # Contado para checagem de novas imagens nas 16 bandas
    for x in range(1, 17):
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