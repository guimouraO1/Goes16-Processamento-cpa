import logging
import os
from multiprocessing import Process
from PIL import Image
import glob
from modules.dirs import get_dirs
import time

dirs = get_dirs()

# Acessando os diretórios usando as chaves do dicionário
dir_out = dirs['dir_out']

def create_gif(band, roi, dir_out):
    try:
        images = []
        for filename in sorted(glob.glob(f"{dir_out}{band}/{band}_*_*_{roi}.png")):
            img = Image.open(filename)
            images.append(img)
        # Salvar como um GIF
        images[0].save(f"{dir_out}{band}/{band}_{roi}.gif", save_all=True, append_images=images[1:], duration=400, loop=0)
    except Exception as e:
        logging.error(f'Erro ao criar GIF: {str(e)}')
        
def create_gif_resized(band, roi, dir_out, resize_factor=0.3):
    try:
        images = []
        for filename in sorted(glob.glob(f"{dir_out}{band}/{band}_*_*_{roi}.png")):
            img = Image.open(filename)
            # Resize the image to half of the original size
            width, height = img.size
            new_size = (int(width * resize_factor), int(height * resize_factor))
            img = img.resize(new_size, resample=Image.BICUBIC)  # Use BICUBIC as the default resampling filter
            images.append(img)
        # Save as a GIF
        images[0].save(f"{dir_out}index_gifs/{band}_{roi}.gif", save_all=True, append_images=images[1:], duration=400, loop=0, optimize=True)
    except Exception as e:
        logging.error(f'Error creating GIF: {str(e)}')
        
        
def process_gif(g_bands, g_br, g_sp, dir_out):
        
    start = time.time()  
    
    gif_br = []
    gif_sp = []

    if g_br:
        logging.info('')
        logging.info('CRIANDO GIF ANIMADO "BR"...')
        for x in range(1, 17):
            b = str(x).zfill(2)
            if g_bands[b]:
                logging.info('Gif BR banda ' + b)
                process = Process(target=create_gif, args=(f'band{b}', "br", dir_out))
                gif_br.append(process)
                process.start()
        for process in gif_br:
            process.join()

    if g_sp:
        logging.info('')
        logging.info('CRIANDO GIF ANIMADO "SP"...')
        for x in range(1, 17):
            b = str(x).zfill(2)
            if g_bands[b]:
                logging.info('Gif SP banda ' + b)
                process = Process(target=create_gif, args=(f'band{b}', "sp", dir_out))
                gif_sp.append(process)
                process.start()
        for process in gif_sp:
            process.join()

    if g_bands['17']:
        if g_br:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO TRUECOLOR "BR"...')
            create_gif("truecolor", "br", dir_out)
        if g_sp:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO TRUECOLOR "SP"...')
            create_gif("truecolor", "sp", dir_out)

    if g_bands['18']:
        if g_br:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO RRQPEF "BR"...')
            create_gif("rrqpef", "br", dir_out)
        if g_sp:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO RRQPEF "SP"...')
            create_gif("rrqpef", "sp", dir_out)

    if g_bands["19"]:
        if g_br:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO GLM "BR"...')
            create_gif("glm", "br", dir_out)

    if g_bands["20"]:
        try:
            if g_br:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO FDCF "BR"...')
                create_gif("ndvi", "br", dir_out)
        except:
            logging.info('Não existe imagens para processar GIF NDVI')

    if g_bands["21"]:
        try:
            if g_br:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO FDCF "BR"...')
                create_gif("fdcf", "br", dir_out)
        except:
            logging.info('Não existe imagens para processar GIF FDCF')

    if g_bands["22"]:
        try:
            if g_br:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO Airmass "BR"...')
                create_gif("airmass", "br", dir_out)
            if g_sp:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO Airmass "SP"...')
                create_gif("airmass", "sp", dir_out)
        except:
            logging.info('Não existe imagens para processar GIF Airmass')

    if g_bands["23"]:
        try:
            if g_br:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO lst "BR"...')
                create_gif("lst", "br", dir_out)
            if g_sp:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO lst "SP"...')
                create_gif("lst", "sp", dir_out)
        except:
            logging.info('Não existe imagens para processar GIF lst')
            
    if g_bands["24"]:
        try:
            if g_br:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO dmw "BR"...')
                create_gif("dmw", "br", dir_out)
            if g_sp:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO dmw "SP"...')
                create_gif("dmw", "sp", dir_out)
        except:
            logging.info('Não existe imagens para processar GIF lst')
            
    if g_bands['17']:
        if g_br:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO LOW QUALITY TRUECOLOR "BR"...')
            create_gif_resized("truecolor", "br", dir_out)

    if g_bands["22"]:
        try:
            if g_br:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO LOW QUALITY Airmass "BR"...')
                create_gif_resized("airmass", "br", dir_out)
        except:
            logging.info('Não existe imagens para processar GIF Airmass')

    if g_bands["23"]:
        try:
            if g_br:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO LOW QUALITY lst "BR"...')
                create_gif_resized("lst", "br", dir_out)
        except:
            logging.info('Não existe imagens para processar GIF lst')
    
    logging.info(f'Tempo de Processamento Gifs: {round((time.time() - start), 2)} segundos.')
