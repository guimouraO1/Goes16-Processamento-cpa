import os
import datetime
import re
from modules.utilities import download_prod

ch13 = 'CG_ABI-L2-CMIPF-M6C13_G16_s20232421050206_e20232421059526_c20232421102264.nc'
dir_in = '/home/guimoura/download_amazon/goes/'


ftime = (datetime.datetime.strptime(ch13[ch13.find("M6C13_G16_s") + 11:ch13.find("_e") - 1], '%Y%j%H%M%S'))

download_prod(datetime.datetime.strftime(ftime, '%Y%m%d%H%M'), 'ABI-L2-RRQPEF', f'{dir_in}rrqpef/')

# Padrão regex para encontrar arquivos desejados
rrqpef = ch13.replace('CG_ABI-L2-CMIPF-M6C13_', 'OR_ABI-L2-RRQPEF-M6_')

print(rrqpef)

print(dir_in + f'/rrqpef/{rrqpef}')