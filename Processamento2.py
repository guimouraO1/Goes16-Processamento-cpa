#!/opt/miniconda3/envs/goes/bin/python3
# -*- coding: utf-8 -*-
# Packages = conda create --name goes -c conda-forge matplotlib netcdf4 cartopy boto3 gdal scipy pandas scp
# Packages = apt install ffmpeg
# ======================================================================================================
# Manipulando imagens GOES-16 NetCDF's
# ===================================# Bibliotecas necessarias =========================================

from osgeo import gdal  # Utilitario para a biblioteca GDAL
from osgeo import osr  # Utilitario para a biblioteca GDAL
import datetime  # Utilitario para datas e horas
import os
import logging  # Utilitario para criar os logs
import paramiko  # Utilitario para gerencia conexao SSH
import scp  # Utilitario para envio de arquivos com SCP
from shutil import copyfile  # Utilitario para copia de arquivos
from shapely.geometry import Point
import cartopy.feature as cfeature # features
import time
from modulos.check_new_images import check_images # Checa as se há novas imagens para processamento
from modulos.processing import process_gif, processing, read_process_file # Processa imagens
from modulos.logs import conf_log, finalize_log_time, start_log_time # Cria os arquivos de logs

gdal.PushErrorHandler('CPLQuietErrorHandler')   # Ignore GDAL warnings

def remove_images(r_bands, r_br, r_sp):
    global dir_temp, dir_in
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

def send_products(s_br, s_sp):
    global dir_out
    logging.info('')
    logging.info('ENVIANDO PRODUTOS PARA O SITE')
    try:
        # Criar objeto cliente SSH
        ssh = paramiko.SSHClient()
        # Carrega as chaves disponiveis no sistema
        ssh.load_system_host_keys()
        # Cria a conexao com o endereco informando, na porta informada e com o usuario informado
        ssh.connect('143.106.227.94', username="cpaunicamp", port=22000)
        # Cria o objeto cliente SCP com a conexao SSH
        scp_client = scp.SCPClient(ssh.get_transport())

        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if s_br:
            logging.info('')
            logging.info('Enviando produtos "BR"')
            # Contador para envio nas 16 bandas
            for x in range(1, 17):
                # Transforma o inteiro contador em string e com 2 digitos
                b = str(x).zfill(2)
                # Cria uma lista com os itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_br.png$"
                ultima_br = [name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_br.png$', name)]
                # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
                ultima_br.sort()
                # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
                ultima_br.reverse()
                # Envia o arquivo "png" mais recente para o site, renomeando no destino
                scp_client.put(f'{dir_out}band{b}/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}_br.png')
                # Cria um arquivo menor do arquivo "png" mais recente
                os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}band{b}/{ultima_br[0]} -vf scale=448:321 {dir_out}band{b}/band{b}.png')
                # Envia o arquivo menor do "png" mais recente para o site
                scp_client.put(f'{dir_out}band{b}/band{b}.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}.png')
                # Envia o arquivo "gif" para o site
                scp_client.put(f'{dir_out}band{b}/band{b}_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}_br.gif')

            # Envio TrueColor
            # Cria uma lista com os itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_br.png$"
            ultima_br = [name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_br.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_br.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}truecolor/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/truecolor/truecolor_br.png')
            # Cria um arquivo menor do arquivo "png" mais recente
            os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}truecolor/{ultima_br[0]} -vf scale=448:321 {dir_out}truecolor/truecolor.png')
            # Envia o arquivo menor do "png" mais recente para o site
            scp_client.put(f'{dir_out}truecolor/truecolor.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/truecolor/truecolor.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}truecolor/truecolor_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/truecolor/truecolor_br.gif')

            # Envio RRQPEF
            # Cria uma lista com os itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_br.png$"
            ultima_br = [name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_br.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_br.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}rrqpef/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/rrqpef/rrqpef_br.png')
            # Cria um arquivo menor do arquivo "png" mais recente
            os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}rrqpef/{ultima_br[0]} -vf scale=448:321 {dir_out}rrqpef/rrqpef.png')
            # Envia o arquivo menor do "png" mais recente para o site
            scp_client.put(f'{dir_out}rrqpef/rrqpef.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/rrqpef/rrqpef.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}rrqpef/rrqpef_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/rrqpef/rrqpef_br.gif')

            # Envio GLM
            # Cria uma lista com os itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_br.png$"
            ultima_br = [name for name in os.listdir(f'{dir_out}glm') if os.path.isfile(os.path.join(f'{dir_out}glm', name)) and re.match('^glm_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_br.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_br.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}glm/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/glm/glm_br.png')
            # Cria um arquivo menor do arquivo "png" mais recente
            os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}glm/{ultima_br[0]} -vf scale=448:321 {dir_out}glm/glm.png')
            # Envia o arquivo menor do "png" mais recente para o site
            scp_client.put(f'{dir_out}glm/glm.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/glm/glm.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}glm/glm_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/glm/glm_br.gif')

            # Envio NDVI
            # Cria uma lista com os itens no diretorio dos produtos ndvi que sao arquivos e se encaixa na expressao regular "^ndvi_.+_.+_br.png$"
            ultima_br = [name for name in os.listdir(f'{dir_out}ndvi') if os.path.isfile(os.path.join(f'{dir_out}ndvi', name)) and re.match('^ndvi_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_br.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_br.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}ndvi/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/ndvi/ndvi_br.png')
            # Cria um arquivo menor do arquivo "png" mais recente
            os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}ndvi/{ultima_br[0]} -vf scale=448:321 {dir_out}ndvi/ndvi.png')
            # Envia o arquivo menor do "png" mais recente para o site
            scp_client.put(f'{dir_out}ndvi/ndvi.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/ndvi/ndvi.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}ndvi/ndvi_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/ndvi/ndvi_br.gif')

            # Envio FDCF
            # Cria uma lista com os itens no diretorio dos produtos fdcf que sao arquivos e se encaixa na expressao regular "^fdcf_.+_.+_br.png$"
            ultima_br = [name for name in os.listdir(f'{dir_out}fdcf') if os.path.isfile(os.path.join(f'{dir_out}fdcf', name)) and re.match('^fdcf_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_br.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_br.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}fdcf/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/fdcf/fdcf_br.png')
            # Cria um arquivo menor do arquivo "png" mais recente
            os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}fdcf/{ultima_br[0]} -vf scale=448:321 {dir_out}fdcf/fdcf.png')
            # Envia o arquivo menor do "png" mais recente para o site
            scp_client.put(f'{dir_out}fdcf/fdcf.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/fdcf/fdcf.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}fdcf/fdcf_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/fdcf/fdcf_br.gif')

        # Se a variavel de controle de processamento do estado de sao paulo for True, realiza o processamento
        if s_sp:
            logging.info('')
            logging.info('Enviando produtos "SP"')
            # Contador para envio nas 16 bandas
            for x in range(1, 17):
                # Transforma o inteiro contador em string e com 2 digitos
                b = str(x).zfill(2)
                # Cria uma lista com os itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_sp.png$"
                ultima_sp = [name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_sp.png$', name)]
                # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
                ultima_sp.sort()
                # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
                ultima_sp.reverse()
                # Envia o arquivo "png" mais recente para o site, renomeando no destino
                scp_client.put(f'{dir_out}band{b}/{ultima_sp[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}_sp.png')
                # Envia o arquivo "gif" para o site
                scp_client.put(f'{dir_out}band{b}/band{b}_sp.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}_sp.gif')

            # Envio TrueColor
            # Cria uma lista com os itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_br.png$"
            ultima_sp = [name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_sp.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_sp.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_sp.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}truecolor/{ultima_sp[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/truecolor/truecolor_sp.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}truecolor/truecolor_sp.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/truecolor/truecolor_sp.gif')

            # Envio RRQPEF
            # Cria uma lista com os itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_sp.png$"
            ultima_sp = [name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_sp.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_sp.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_sp.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}rrqpef/{ultima_sp[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/rrqpef/rrqpef_sp.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}rrqpef/rrqpef_sp.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/rrqpef/rrqpef_sp.gif')

    except TimeoutError as e_timeout:
        logging.info('')
        logging.info(f'FALHA AO ENVIAR PRODUTOS PARA O SITE - {e_timeout}')
    except paramiko.SSHException as e_ssh:
        logging.info('')
        logging.info(f'FALHA AO ENVIAR PRODUTOS PARA O SITE - {e_ssh}')
    except scp.SCPException as e_scp:
        logging.info('')
        logging.info(f'FALHA AO ENVIAR PRODUTOS PARA O SITE - {e_scp}')



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

# configura o log
conf_log(arq_log)

# Pega a hora para fazer o calculo do tempo de processamento
start = time.time()

# Realizando a checagem de novas imagens para processamento
bands = check_images(bands, dir_in, dir_temp)

# Realiza etapas de processamento se houver alguma nova imagem
if bands["01"] or bands["02"] or bands["03"] or bands["04"] or bands["05"] or bands["06"] or bands["07"] or bands["08"] or \
        bands["09"] or bands["10"] or bands["11"] or bands["12"] or bands["13"] or bands["14"] or bands["15"] or bands["16"]:
   
    # Realiza processamento da imagens
    bands = processing(bands, br, sp)
    # Remove imagens processadas
    #remove_images(bands, br, sp)
    # Controle de quantidade de produtos devemos manter para produção do gif
    quantity_products()
    # Realiza processamento do gif
    process_gif(bands, br, sp)

else:
    logging.info("")
    logging.info("SEM NOVAS IMAGENS PARA PROCESSAMENTO")


# Finaliza o script
finalize_log_time(start)



# ======================================================================================================
