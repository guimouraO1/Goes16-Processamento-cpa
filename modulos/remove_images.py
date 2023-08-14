import logging
import os

dir_main = "/home/guimoura/download_amazon/"
dir_temp = dir_main + "temp/"

def read_process_file(banda):
    # Le o arquivo de processamento e retorna a lista
    with open(f'{dir_temp}{banda}_process.txt', 'r') as fo:
        return fo.readlines()

def remove_images(r_bands, r_br, r_sp, dir_temp, dir_in):

    logging.info("")
    logging.info('REMOVENDO IMAGENS PROCESSADAS')

    # Contador para remover imagens nas 16 bandas
    for x in range(1, 17):
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

    if r_bands["17"]:
        logging.info(f'Removendo arquivo de processamento TRUECOLOR')
        # Remove o arquivo de processamento truecolor
        os.remove(f'{dir_temp}band17_process.txt')

    if r_bands["18"]:
        # Le o arquivo de processamento da banda
        img = read_process_file(f'band18')
        logging.info(f'Removendo imagens RRQPEF')
        for i in img:
            # Remove possiveis espacos vazios no inicio ou final da string
            i = i.strip().split(';')
            try:
                # Remove a imagem
                os.remove(f'{dir_in}rrqpef/{i[0]}')
            except FileNotFoundError as fnfe:
                # Realiza o log do erro
                logging.info(f'Erro Arquivo - FileNotFoundError - {i[0]}')
                logging.info(str(fnfe))
            try:
                if r_br:
                    os.remove(f'{dir_in}rrqpef/{i[0].replace(".nc", "_reproj_br.nc")}')
            except FileNotFoundError as fnfe:
                # Realiza o log do erro
                logging.info(f'Erro Arquivo - FileNotFoundError - {i[0].replace(".nc", "_reproj_br.nc")}')
                logging.info(str(fnfe))
            try:
                if r_sp:
                    os.remove(f'{dir_in}rrqpef/{i[0].replace(".nc", "_reproj_sp.nc")}')
            except FileNotFoundError as fnfe:
                # Realiza o log do erro
                logging.info(f'Erro Arquivo - FileNotFoundError - {i[0].replace(".nc", "_reproj_sp.nc")}')
                logging.info(str(fnfe))
        logging.info(f'Removendo arquivo de processamento RRQPEF')
        # Remove o arquivo de processamento rrqpef
        os.remove(f'{dir_temp}band18_process.txt')

    if r_bands["19"]:
        # Le o arquivo de processamento da banda
        img = read_process_file("band19")
        logging.info(f'Removendo imagens GLM')
        # Para cada imagem no arquivo, realiza a remocao
        for i in img:
            # Remove possiveis espacos vazios no inicio ou final da string e separa cada termo como um elemento
            i = i.strip().split(';')
            # Coletando lista de GLM no indice 1
            j = i[1].replace("'", "").strip('][').split(', ')
            for g in j:
                try:
                    # Remove a imagem
                    os.remove(f'{dir_in}glm/{g}')
                except FileNotFoundError as fnfe:
                    # Realiza o log do erro
                    logging.info(f'Erro Arquivo - FileNotFoundError - {g}')
                    logging.info(str(fnfe))
        logging.info(f'Removendo arquivo de processamento do GLM')
        # Remove o arquivo de processamento
        os.remove(f'{dir_temp}band19_process.txt')

    if r_bands["20"]:
        logging.info(f'Removendo arquivo de processamento NDVI')
        # Remove o arquivo de processamento ndvi
        os.remove(f'{dir_temp}band20_process.txt')

    if r_bands["21"]:
        # Le o arquivo de processamento da banda
        img = read_process_file(f'band21')
        logging.info(f'Removendo imagens FDCF')
        for i in img:
            # Remove possiveis espacos vazios no inicio ou final da string
            i = i.strip().split(';')
            try:
                # Remove a imagem
                os.remove(f'{i[0]}')
            except FileNotFoundError as fnfe:
                # Realiza o log do erro
                logging.info(f'Erro Arquivo - FileNotFoundError - {i[0]}')
                logging.info(str(fnfe))
        logging.info(f'Removendo arquivo de processamento FDCF')
        # Remove o arquivo de processamento fdcf
        os.remove(f'{dir_temp}band21_process.txt')
