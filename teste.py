import json


# Função para abrir o arquivo "oldBands.json" e retornar a lista de "oldImagesName".
def openOld():
    with open('oldBands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
        return oldImages
    
    
    
old_bands = openOld()

print(type(old_bands['03']))