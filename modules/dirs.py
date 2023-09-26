import datetime

# Modifique aqui os diretórios para uso
def get_dirs():
    dir_main = '/home/guimoura/processamento/'
    dirs = {
        'dir_in': 'home/goes/',
        'dir_main': dir_main,
        'dir_out': dir_main + 'output/',
        'dir_shapefiles': dir_main + 'shapefiles/',
        'dir_colortables': dir_main + 'colortables/',
        'dir_logos': dir_main + 'logos/',
        'dir_temp': dir_main + 'temp/',
        'arq_log': dir_main + f'logs/Processamento-GOES_{str(datetime.date.today())}.log'
    }
    
    return dirs


