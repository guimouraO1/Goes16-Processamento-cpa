import logging
import os

def read_process_file(banda, dir_temp):
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
            img = read_process_file(f'band{b}', dir_temp)
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