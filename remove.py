import logging
import os
import json

def openOld():
    with open('oldBands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
        return oldImages
    
def remove_images(dir_in):

    logging.info("")
    logging.info('REMOVENDO IMAGENS PROCESSADAS')

    # Contador para remover imagens nas 16 bandas
    for x in range(1, 17):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        oldBands = openOld()        
        try:
            # Remove a imagem
            os.remove(f'{dir_in}/band{b}/{oldBands[b]}')
            logging.info(f'imagem {oldBands[b]} removida com sucesso!')
        except:
            print(f'Sem imagens para remover {oldBands[b]}')

