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
        images[0].save(f"{dir_out}index_gif/{band}_{roi}.gif", save_all=True, append_images=images[1:], duration=400, loop=0,optimize=True)
    except Exception as e:
        logging.error(f'Error creating GIF: {str(e)}')

def process_gif(g_bands, g_br, dir_out):
    
    start = time.time()  
    
    if g_bands['17']:
        if g_br:
            logging.info('')
            logging.info('CRIANDO GIF ANIMADO TRUECOLOR "BR"...')
            create_gif_resized("truecolor", "br", dir_out)

    if g_bands["22"]:
        try:
            if g_br:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO Airmass "BR"...')
                create_gif_resized("airmass", "br", dir_out)
        except:
            logging.info('Não existe imagens para processar GIF Airmass')

    if g_bands["23"]:
        try:
            if g_br:
                logging.info('')
                logging.info('CRIANDO GIF ANIMADO lst "BR"...')
                create_gif_resized("lst", "br", dir_out)
        except:
            logging.info('Não existe imagens para processar GIF lst')

    logging.info(f'Tempo de Processamento Gifs: {round((time.time() - start), 2)} segundos.')

bands = {}
# Todas as bandas da 01 a 21 recebem False      bands = {"01": False, "02": False......
for num in range(1, 24):
    b = str(num).zfill(2)
    bands[f'{b}'] = True
br = True
sp = True

process_gif(bands, br, dir_out)
