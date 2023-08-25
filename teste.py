import sys
import json
import zlib
# Função para abrir o arquivo.json
def openOld():
    with open('oldBands.json', 'r') as jsonOld:
        oldImages = json.load(jsonOld)['oldImagesName']
        return oldImages

old_img = openOld()

# Encontre o índice do elemento na lista
indice_elemento = old_img['01'].replace(".nc", "_reproj_br.nc")

# Calcule o tamanho em bits do elemento desejado
tamanho_elemento_bits = sys.getsizeof(indice_elemento) * 8

print(sys.getsizeof('CG_ABI-L2-CMIPF-M6C01_G16_s20232371200205_e20232371209513_c20232371211042_reproj_br.nc') * 8)

print(f'Tamanho do elemento {indice_elemento} em bits: {tamanho_elemento_bits} bits')

print(type(old_img['01']))


