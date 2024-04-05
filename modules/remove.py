import logging
import os
import json
import re
from modules.dirs import get_dirs

# ============================================# Diretórios ========================================= #
dirs = get_dirs()

# Acessando os diretórios usando as chaves do dicionário
dir_main = dirs['dir_main']
dir_in = dirs['dir_in']

def abrir_old_json():
    global dir_main
    with open(f'{dir_main}old_bands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
        return oldImages

# Remove os arquivos netCDF que foram processados
def remover_imagens(bands, dir_in):
    
    logging.info("")
    logging.info('REMOVENDO IMAGENS PROCESSADAS')
    oldBands = abrir_old_json()
    
    # Contador para remover imagens nas 16 bandas
    for x in range(1, 17):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
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

    # Controle de produtos de rrqpef
    if bands['18']:
        ch18 = oldBands['18']
        logging.info(f'Removendo imagens rrqpef')
        try:
            os.remove(f'{dir_in}rrqpef/{ch18}')
            os.remove(f'{dir_in}rrqpef/{ch18}'.replace(".nc", "_reproj_br.nc"))
            os.remove(f'{dir_in}rrqpef/{ch18}'.replace(".nc", "_reproj_sp.nc"))
        except FileNotFoundError as e:
            pasta = os.listdir(f'{dir_in}rrqpef/')
            for item in pasta:
                os.remove(os.path.join(f'{dir_in}rrqpef/', item))
            logging.info(str(e))

    # Controle de produtos de glm
    if bands['19'] or bands['13']:
        logging.info(f'Removendo imagens GLM')
        try:
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
                logging.info("Produtos glm em excesso removidos com sucesso")
            else:
                logging.info("Quantidade de produtos glm dentro do limite")
        except:
            logging.info('Ocorreu um erro ao tentar apagar produtos glm')
            
   
    if bands['21']:
        ch21 = oldBands['21']
        logging.info(f'Removendo imagens FDCF')
        try:
            os.remove(f'{dir_in}fdcf/{ch21}')
        except FileNotFoundError as e:
            logging.info(str(e))
    
    
    if bands['23']:    
        ch23 = oldBands['23']
        logging.info(f'Removendo imagens LST')
        try:
            os.remove(f'{dir_in}lst/{ch23}')
        except FileNotFoundError as e:
            logging.info(str(e))
 
    #Remover DMW
    if bands['24']:    
        ch24 = oldBands['24']
        logging.info(f'Removendo imagens DMW')
        try:
            os.remove(f'{dir_in}dmw/{ch24}')
        except FileNotFoundError as e:
            logging.info(str(e))
            
    #Remover SST
    if bands['25']:    
        ch25 = oldBands['25']
        logging.info(f'Removendo imagens DMW')
        try:
            os.remove(f'{dir_in}sst/{ch25}')
        except FileNotFoundError as e:
            logging.info(str(e))