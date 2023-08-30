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
        
    if bands['19']:
        # pega lista de itens que sobraram em glm 
        images = [f for f in os.listdir(f'{dir_in}glm') if os.path.isfile(os.path.join(f'{dir_in}glm', f)) and re.match('^OR_GLM-L2-LCFA_G16_s.+_e.+_c.+.nc$', f)]
        images.sort()
        # band13 para comparação
        ch13 = oldBands['13']
        # Lista para excluir arquivos
        lista_exclusao = []
        # Pegando data inicio e fim para comparação
        ch13_data = (datetime.datetime.strptime(ch13[ch13.find("M6C13_G16_s") + 11:ch13.find("_e") - 1], '%Y%j%H%M%S'))
        date_ini = datetime.datetime(ch13_data.year, ch13_data.month, ch13_data.day, ch13_data.hour, ch13_data.minute)
        # Populando lista para exclusão de imagens fora do período
        for x in images:
            xtime = (datetime.datetime.strptime(x[x.find("GLM-L2-LCFA_G16_s") + 17:x.find("_e") - 1], '%Y%j%H%M%S'))
            if date_ini > xtime:
                lista_exclusao.append(x)
            else:
                continue
        # Exclui a lista de arquivos glm que já não vão ser mais processados pois estão antes da data
        pasta_glm = dir_in + 'glm/'
        try:
            [os.remove(os.path.join(pasta_glm, arquivo)) for arquivo in os.listdir(pasta_glm) if arquivo in lista_exclusao]
            logging.info('Arquivos da glm_lista antes da data {date_ini}foram excluídos com sucesso! ')
        except:
            logging.info('Não foi possível excluír arquivos glm fora da data limite ')