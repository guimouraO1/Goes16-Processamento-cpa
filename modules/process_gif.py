import logging
import os
from multiprocessing import Process  # Utilitario para multiprocessamento


def process_gif(g_bands, g_br, g_sp, dir_out):

    # Cria lista vazia para controle do processamento paralelo
    gif_br = []
    # Cria lista vazia para controle do processamento paralelo
    gif_sp = []

    # Chamada de sistema para o software ffmpeg realizar a criacao do gif animado
    def create_gif_cmi(banda, roi):
        os.system(f'/usr/bin/ffmpeg -y -v warning -framerate 4 -pattern_type glob -i "{dir_out}band{banda}/band{banda}_*_*_{roi}.png" "{dir_out}band{banda}/band{banda}_{roi}.gif"')

    # Chamada de sistema para o software ffmpeg realizar a criacao do gif animado
    def create_gif(file_type, roi):
        os.system(f'/usr/bin/ffmpeg -y -v warning -framerate 4 -pattern_type glob -i "{dir_out}{file_type}/{file_type}_*_*_{roi}.png" "{dir_out}{file_type}/{file_type}_{roi}.gif"')

    # Se a variavel de controle de processamento do brasil for True, cria o gif
    if g_br:
        logging.info('')
        logging.info('CRIANDO GIF ANIMADO "BR"...')
        # Contador para gif nas 16 bandas
        for x in range(1, 17):
            # Transforma o inteiro contador em string e com 2 digitos
            b = str(x).zfill(2)
            # Se o dicionario de controle das bandas apontar True para essa banda, cria o gif
            if g_bands[b]:
                logging.info('Gif BR banda ' + b)
                # Cria o processo com a funcao gif
                process = Process(target=create_gif_cmi, args=(b, "br"))
                # Adiciona o processo na lista de controle do processamento paralelo
                gif_br.append(process)
                # Inicia o processo
                process.start()
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in gif_br:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()

    # Se a variavel de controle de processamento do estado de sao paulo for True, cria o gif
    if g_sp:
        logging.info('')
        logging.info('CRIANDO GIF ANIMADO "SP"...')
        # Contador para gif nas 16 bandas
        for x in range(1, 17):
            # Transforma o inteiro contador em string e com 2 digitos
            b = str(x).zfill(2)
            # Se o dicionario de controle das bandas apontar True para essa banda, cria o gif
            if g_bands[b]:
                logging.info('Gif SP banda ' + b)
                # Cria o processo com a funcao gif
                process = Process(target=create_gif_cmi, args=(b, "sp"))
                # Adiciona o processo na lista de controle do processamento paralelo
                gif_br.append(process)
                # Inicia o processo
                process.start()
        # Looping de controle que pausa o processamento principal ate que todos os processos da lista de controle do processamento paralelo sejam finalizados
        for process in gif_sp:
            # Bloqueia a execução do processo principal ate que o processo cujo metodo de join() é chamado termine
            process.join()

    # Se bands 17 for true processa True color gif br e sp
    if g_bands['17']:
        if g_br:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO TRUECOLOR "BR"...')
            # Cria o processo com a funcao gif
            create_gif("truecolor", "br")
        if g_sp:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO TRUECOLOR "SP"...')
            # Cria o processo com a funcao gif
            create_gif("truecolor", "sp")
    
    if g_bands['18']:
        if g_br:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO RRQPEF "BR"...')
            # Cria o processo com a funcao gif
            create_gif("rrqpef", "br")
        if g_br:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO RRQPEF "SP"...')
            # Cria o processo com a funcao gif
            create_gif("rrqpef", "sp")