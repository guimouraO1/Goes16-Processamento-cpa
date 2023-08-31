import logging
import os
import json
import shutil
import datetime
import re

# Abre o json para pegar o nome dos arquivos
def abrir_old_json():
    with open('oldBands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
        return oldImages

# Remove os arquivos netCDF que foram processados
def remover_imagens(bands, dir_in):
    logging.info("")
    logging.info('REMOVENDO IMAGENS PROCESSADAS')
    # Contador para remover imagens nas 16 bandas
    for x in range(1, 17):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        oldBands = abrir_old_json()
        if bands[b] == True:       
            try:
                reprojbr = str(oldBands[b])
                reprojsp = str(oldBands[b])
                # Remove os arquivos netCDF que já foram processados
                os.remove(f'{dir_in}band{b}/{oldBands[b]}')
                os.remove(f'{dir_in}band{b}/{reprojbr.replace(".nc", "_reproj_br.nc")}')
                os.remove(f'{dir_in}band{b}/{reprojsp.replace(".nc", "_reproj_sp.nc")}')
                logging.info(f'Arquivo {oldBands[b]} removido com sucesso!')
            except:
                print(f'Sem imagens para remover {oldBands[b]}')


    if bands['18']:
        logging.info('Removendo netCDF-rrqpef processadas')
        shutil.rmtree(f'{dir_in}rrqpef/')
        os.mkdir(f'{dir_in}rrqpef/')

    # Controle de produtos de glm
    if bands['19'] or bands['13']:
        aux = False
        # Verificar produtos glm
        # Verifica a quantidade de itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_br.png$"
        result_br = len([name for name in os.listdir(f'{dir_in}glm') if os.path.isfile(os.path.join(f'{dir_in}glm', name)) and re.match('^OR_GLM-L2-LCFA_G16_.+_.+.nc$', name)])
        # Subtrai o limite de imagens BR
        result_br -= 48
        # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
        if result_br > 0:
            aux = True
            # Cria uma lista com os itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_br.png$"
            prod_br = [name for name in os.listdir(f'{dir_in}glm') if os.path.isfile(os.path.join(f'{dir_in}glm', name)) and re.match('^OR_GLM-L2-LCFA_G16_.+_.+.nc$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            prod_br.sort()
            # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
            for y in range(0, result_br):
                # Remove arquivo em excesso
                os.remove(f'{dir_in}glm/{prod_br[y]}')
        if aux:
            logging.info("Produtos em excesso removidos com sucesso")
        else:
            logging.info("Quantidade de produtos dentro do limite")