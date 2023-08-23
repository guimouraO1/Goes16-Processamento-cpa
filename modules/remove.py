import logging
import os
import json

def openOld():
    with open('oldBands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
        return oldImages
    
def removeImagens(bands, dir_in):
    logging.info("")
    logging.info('REMOVENDO IMAGENS PROCESSADAS')
    # Contador para remover imagens nas 16 bandas
    for x in range(1, 17):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        oldBands = openOld()
        if bands[b] == True:       
            try:
                reprojbr = str(oldBands[b])
                reprojsp = str(oldBands[b])
                # Remove a imagem
                os.remove(f'{dir_in}/band{b}/{oldBands[b]}')
                os.remove(f'{dir_in}band{b}/{reprojbr.replace(".nc", "_reproj_br.nc")}')
                os.remove(f'{dir_in}band{b}/{reprojsp.replace(".nc", "_reproj_sp.nc")}')
                logging.info(f'Arquivo {oldBands[b]} removido com sucesso!')
            except:
                print(f'Sem imagens para remover {oldBands[b]}')