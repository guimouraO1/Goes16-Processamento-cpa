import os
import re
import json

# Cria o arquivo de processamento process.txt
def processTxt(newImagesName, dir_temp):    
    # Cria o arquivo band??_process.txt com as imagens para processamento
    for banda, process_list in newImagesName.items():
        with open(f'{dir_temp}band{banda}_process.txt', 'w') as process:
            # Escreve as imagens da lista no arquivo, sendo cada uma em uma linha
            process.writelines(map(lambda f: f + '\n', process_list))

# Atualiza os dados em oldBands apenas as imagens que já foram processadas
def atualizarOldBands():
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
        
# checa se há imagens novas para processamento
def checarImagens(bands, dir_in, dir_temp):
    newImagesName = {}

    # Abre oldBands para checar as imagens antigas com as novas
    with open('oldBands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
    
    # inicia o loop para comparar as imagens
    for x in range(1, 17):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        # Pega todos as imagens que tiverem na pasta band{''}
        imagens = [f for f in os.listdir(f'{dir_in}band{b}') if os.path.isfile(os.path.join(f'{dir_in}band{b}', f)) and re.match('^CG_ABI-L2-CMIPF-M[0-9]C[0-1][0-9]_G16_s.+_e.+_c.+.nc$', f)]
        
       # ++++++++++++++++++++++++ arrumar por data e colocar condição de 1 arquivo por vez     
        
        # Se existir imagens e elas forem de uma data maior que a anterior
        if imagens and imagens > oldImages[b]:
            print(f'Novas imagens para a banda {b}')
            newImagesName[b] = imagens
            bands[b] = True
        else:
            print(f'Sem imagens para a banda {b}')
            bands[b] = False

    # Escreve em newBands.json as novas imagens
    with open('newBands.json', 'w') as json_file:
        json.dump({'newImagesName': newImagesName}, json_file, indent=4)
    
    # Cria o arquivo process.txt
    processTxt(newImagesName, dir_temp)
    
    # Atualiza os dados em oldBands apenas as imagens que já foram processadas
    atualizarOldBands()
    
    return bands

atualizarOldBands()