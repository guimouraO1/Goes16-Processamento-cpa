#!/opt/miniconda3/envs/goes/bin/python3
# -*- coding: utf-8 -*-

from osgeo import gdal  # Utilitario para a biblioteca GDAL
from modules.utilities import download_prod  # Funcao para download dos produtos do goes disponiveis na amazon
from modules.utilities import download_cmi
import datetime
from datetime import timedelta
import logging
import os
import time

gdal.PushErrorHandler('CPLQuietErrorHandler')   # Ignore GDAL warnings 

processing_start_time = time.time()

#Diretórios
dir_in = "/Scripts/goes16/processamento/goes/"
arq_log = "/Scripts/goes16/processamento/logs/" + str(datetime.date.today()) + ".log"
# Configurando log
logging.basicConfig(filename=arq_log, level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="ls%m/%Y %H:%M:%S")
# Capturando data/hora inicio
inicio = datetime.datetime.now().strftime("%d/%m/%Y-%H:%M:%S")
logging.info("")
logging.info("")
logging.info("=========================================================================================================")
logging.info("=                             DOWNLOAD IMAGENS GOES AS " + inicio + "                           =")
logging.info("=========================================================================================================")
logging.info("")

#Inicia contador para o tempo gasto de processamento
start = time.time()

#Captura hora atual em UTC para download no site da Amazon
data_hora_atual = datetime.datetime.utcnow()

#Atrasa 10 min para entrar em conformidade com Amazon
data_10_min = datetime.datetime.strftime(data_hora_atual-timedelta(minutes=10),'%Y%m%d%H%M')

#Correção para poder fazer download em qualquer horário
data_hora_download_file = data_10_min[0:11]+ '0'

# #Download band01
download_cmi(data_hora_download_file,1,dir_in +'band01')

# #Download band02
download_cmi(data_hora_download_file,2,dir_in +'band02')

# #Download band03
download_cmi(data_hora_download_file,3,dir_in +'band03')

# #Download band04
download_cmi(data_hora_download_file,4,dir_in +'band04')

# #Download band05
download_cmi(data_hora_download_file,5,dir_in +'band05')

# #Download band06
download_cmi(data_hora_download_file,6,dir_in +'band06')

# #Download band07
download_cmi(data_hora_download_file,7,dir_in +'band07')

# #Download band08
download_cmi(data_hora_download_file,8,dir_in +'band08')

# #Download band09
download_cmi(data_hora_download_file,9,dir_in +'band09')

# #Download band10
download_cmi(data_hora_download_file,10,dir_in +'band10')

# #Download band11
download_cmi(data_hora_download_file,11,dir_in +'band11')

# #Download band12
download_cmi(data_hora_download_file,12,dir_in +'band12')

# #Download band13
download_cmi(data_hora_download_file,13,dir_in +'band13')

# #Download band14
download_cmi(data_hora_download_file,14,dir_in +'band14')

# #Download band15
download_cmi(data_hora_download_file,15,dir_in +'band15')

# #Download band 16
download_cmi(data_hora_download_file,16,dir_in +'band16')

#Download band01
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band1 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,1,dir_in +'band01')
            break
        except:
            print("Sem imagens para Band1 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band1 no horário {startb1}")



#Download band02
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band2 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,2,dir_in +'band02')
            break
        except:
            print("Sem imagens para Band2 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band1 no horário {startb1}")



#Download band03
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band3 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,3,dir_in +'band03')
            break
        except:
            print("Sem imagens para Band3 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band3 no horário {startb1}")



#Download band04
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band4 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,4,dir_in +'band04')
            break
        except:
            print("Sem imagens para Band4 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band4 no horário {startb1}")



#Download band05
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band5 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,5,dir_in +'band05')
            break
        except:
            print("Sem imagens para Band5 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band5 no horário {startb1}")



#Download band06
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band6 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,6,dir_in +'band06')
            break
        except:
            print("Sem imagens para Band6 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band6 no horário {startb1}")



#Download band07
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band7 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,7,dir_in +'band07')
            break
        except:
            print("Sem imagens para Band7 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band7 no horário {startb1}")



#Download band08
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band8 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,8,dir_in +'band08')
            break
        except:
            print("Sem imagens para Band8 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band8 no horário {startb1}")



#Download band09
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band9 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,9,dir_in +'band09')
            break
        except:
            print("Sem imagens para Band9 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band9 no horário {startb1}")



#Download band10
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band10 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,10,dir_in +'band10')
            break
        except:
            print("Sem imagens para Band10 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band10 no horário {startb1}")



#Download band11
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band11 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,11,dir_in +'band11')
            break
        except:
            print("Sem imagens para Band11 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band11 no horário {startb1}")



#Download band12
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band12 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,12,dir_in +'band12')
            break
        except:
            print("Sem imagens para Band12 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band12 no horário {startb1}")



#Download band13
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band13 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,13,dir_in +'band13')
            break
        except:
            print("Sem imagens para Band13 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band13 no horário {startb1}")



#Download band14
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band14 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,14,dir_in +'band14')
            break
        except:
            print("Sem imagens para Band14 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band14 no horário {startb1}")



#Download band15
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band15 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,15,dir_in +'band15')
            break
        except:
            print("Sem imagens para Band15 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band15 no horário {startb1}")



#Download band16
logging.info("")
startb1 = datetime.datetime.now()
logging.info(f'Download Band16 as {startb1.strftime("%d/%m/%Y-%H:%M:%S")}')
for i in range(0, 3):
        try:
            download_cmi(data_hora_download_file,16,dir_in +'band16')
            break
        except:
            print("Sem imagens para Band16 no horário")
            print("Tentativa de Download {} de 3" .format(i + 1))
            logging.info(f"Sem imagens para Band16 no horário {startb1}")


def finalize(s):
    # Capturando data/hora final
    fim = datetime.datetime.now().strftime("%d/%m/%Y-%H:%M:%S")
    # Realiza o log do calculo do tempo de execucao
    logging.info("")
    logging.info("Tempo gasto " + str(round(time.time() - s, 4)) + ' segundos')
    logging.info("")
    logging.info("=========================================================================================================")
    logging.info("=                         PROCESSAMENTO ENCERRADO GOES AS " + fim + "                           =")
    logging.info("=========================================================================================================")
    logging.info("")
    logging.info("")

finalize(start)

print(round(time.time() - processing_start_time, 4) / 60)


