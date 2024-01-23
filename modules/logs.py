from datetime import datetime
import logging
import time

def conf_log(arq_log):
        # Configurando log
        logging.basicConfig(filename=arq_log, level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S")
        # Capturando data/hora inicio
        inicio = datetime.now().strftime("%d/%m/%Y-%H:%M:%S")
        logging.info("")
        logging.info("")
        logging.info("==============================================================================================")
        logging.info("=                         PROCESSANDO IMAGENS GOES AS " + inicio + "                         =")
        logging.info("==============================================================================================")
        logging.info("")

def finalize_log_time(s):
        # Capturando data/hora final
        fim = datetime.now().strftime("%d/%m/%Y-%H:%M:%S")
        # Realiza o log do calculo do tempo de execucao
        logging.info("")
        logging.info("Tempo gasto " + str(round(time.time() - s, 4)) + ' segundos')
        logging.info("")
        logging.info("===============================================================================================")
        logging.info("=                         PROCESSAMENTO ENCERRADO GOES AS " + fim + "                         =")
        logging.info("===============================================================================================")
        logging.info("")
        logging.info("")

def conf_log_D(arq_log):
        # Configurando log
        logging.basicConfig(filename=arq_log, level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%d/%m/%Y %H:%M:%S")
        # Capturando data/hora inicio
        inicio = datetime.now().strftime("%d/%m/%Y-%H:%M:%S")
        logging.info("")
        logging.info("")
        logging.info("===========================================================================================")
        logging.info("=                         DOWNLOAD IMAGENS GOES AS " + inicio + "                         =")
        logging.info("===========================================================================================")       
        logging.info("")


def finalize_log_time_D(s):
        # Capturando data/hora final
        fim = datetime.now().strftime("%d/%m/%Y-%H:%M:%S")
        # Realiza o log do calculo do tempo de execucao
        logging.info("")
        logging.info("Tempo gasto " + str(round(time.time() - s, 4)) + ' segundos')
        logging.info("")
        logging.info("==========================================================================================")
        logging.info("=                         DOWNLOAD ENCERRADO GOES AS " + fim + "                         =")
        logging.info("==========================================================================================")
        logging.info("")
        logging.info("")