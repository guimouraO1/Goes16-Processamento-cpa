import os
import re
import json
import logging

# Cria o arquivo de processamento process.txt
def processTxt(newImagesName, dir_temp):
    logging.info(f'Criando arquivo process.txt para processamento') 
    # Cria o arquivo band??_process.txt com as imagens para processamento
    for banda, process_list in newImagesName.items():
        with open(f'{dir_temp}band{banda}_process.txt', 'w') as process:
            # Escreve as imagens da lista no arquivo, sendo cada uma em uma linha
            process.writelines(map(lambda f: f + '\n', process_list))

# Atualiza os dados em oldBands apenas as imagens que já foram processadas
def atualizarOldBands():
    logging.info(f'Atualizando oldBands.json') 
    
    # abre old bands e atribui a var oldImages
    with open('oldBands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
     # abre new bands e atribui a var newImages
    with open('newBands.json', 'r') as jsonNew:
        newImages = json.load(jsonNew)['newImagesName']
    # Atualiza em oldBands.json apenas as bands que têm mudanças
    for key, value in newImages.items():
        if key in oldImages and oldImages[key] != value:
            oldImages[key] = value
    # Salva as alterações de volta no arquivo oldBands.json
    with open('oldBands.json', 'w') as jsonOld:
        json.dump({'oldImagesName': oldImages}, jsonOld, indent=4)
        

def removerTodosExceto(nome_arquivo, pasta):
    [os.remove(os.path.join(pasta, arquivo)) for arquivo in os.listdir(pasta) if arquivo != nome_arquivo]


# checa se há imagens novas para processamento
def checarImagens(bands, dir_in, dir_temp):
    logging.info("VERIFICANDO NOVAS IMAGENS")
    # Cria o dicionario de dados
    newImagesName = {}
    # inicia o loop para comparar as imagens
    for x in range(1, 17):
        # Lista para verificação do arquivo mais recente
        band = []
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        # Pega todos as imagens que tiverem na pasta band{''}
        imagens = list(filter(lambda f: os.path.isfile(os.path.join(f'{dir_in}band{b}', f)) and re.match('^CG_ABI-L2-CMIPF-M[0-9]C[0-1][0-9]_G16_s.+_e.+_c.+.nc$', f), os.listdir(f'{dir_in}band{b}')))
        # Loop da band?? para popular lista e depois comparar
        if imagens:
            logging.info(f'Novas imagens para o dia band{b}')  
            for image in imagens:
                band.append(image)
            # Extraindo nome da band mais recente
            latestBand = max(band)
            # Removendo outras bands a mais
            removerTodosExceto(latestBand, f'{dir_in}band{b}/')
            bands[b] = True
            newImagesName[b] = latestBand
        else:
            logging.info(f'Sem imagens para o dia band{b}')
            bands[b] = False

    # Escreve o dicionario em um arquivo.json para salvar as novas imagens para processamento
    with open('newBands.json', 'w') as json_file:
        json.dump({'newImagesName': newImagesName}, json_file, indent=4)
    
    # Cria o arquivo process.txt para o processamento
    processTxt(newImagesName, dir_temp)
    
    # Atualiza os dados em oldBands apenas as imagens que já foram processadas
    atualizarOldBands()
    
    return bands



