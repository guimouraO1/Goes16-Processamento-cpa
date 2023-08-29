import os
import datetime
import re

ch13 = 'CG_ABI-L2-CMIPF-M6C13_G16_s20232411440206_e20232411449526_c20232411452272.nc'
dir_in = '/home/guimoura/download_amazon/goes/'



# Cria uma lista com os itens presentes no diretório da banda que são arquivos e terminam com ".nc"
images = [f for f in os.listdir(f'{dir_in}glm') if os.path.isfile(os.path.join(f'{dir_in}glm', f)) and re.match('^OR_GLM-L2-LCFA_G16_s.+_e.+_c.+.nc$', f)]
images.sort()
            
glm_list = []

ch13 = (datetime.datetime.strptime(ch13[ch13.find("M6C13_G16_s") + 11:ch13.find("_e") - 1], '%Y%j%H%M%S'))
date_ini = datetime.datetime(ch13.year, ch13.month, ch13.day, ch13.hour, ch13.minute)
date_end = datetime.datetime(ch13.year, ch13.month, ch13.day, ch13.hour, ch13.minute) + datetime.timedelta(minutes=9, seconds=59)


for x in images:
    xtime = (datetime.datetime.strptime(x[x.find("GLM-L2-LCFA_G16_s") + 17:x.find("_e") - 1], '%Y%j%H%M%S'))
    if date_ini <= xtime <= date_end:
        glm_list.append(x)
    else:
        continue

print(glm_list)