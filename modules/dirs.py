import datetime

# Modifique aqui os diretórios para uso
def get_dirs():
    dir_main = '/home/guimoura/download_amazon/'
    # dir_main = '/mnt/c/Users/Gui/OneDrive/Documentos/processamentopy/'
    dirs = {
        'dir_in': dir_main + 'goes/',
        'dir_main': dir_main,
        'dir_out': dir_main + 'output/',
        'dir_libs': dir_main + 'libs/',
        'dir_shapefiles': dir_main + 'shapefiles/',
        'dir_colortables': dir_main + 'colortables/',
        'dir_logos': dir_main + 'logos/',
        'dir_temp': dir_main + 'temp/',
        'arq_log': dir_main + 'logs/Processamento-GOES_' + str(datetime.date.today()) + '.log'
    }
    
    return dirs


