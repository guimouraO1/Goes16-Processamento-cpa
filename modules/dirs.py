import datetime

def get_dirs():
    dir_main = '/home/guimoura/download_amazon/'
    dirs = {
        'dir_in': '/home/guimoura/download_amazon/goes/',
        'dir_main': dir_main,
        'dir_out': dir_main + 'output/',
        'dir_libs': dir_main + 'libs/',
        'dir_shapefiles': dir_main + 'shapefiles/',
        'dir_colortables': dir_main + 'colortables/',
        'dir_logos': dir_main + 'logos/',
        'dir_temp': dir_main + 'temp/',
        'arq_log': '/home/guimoura/download_amazon/logs/Processamento-GOES_' + str(datetime.date.today()) + '.log'
    }
    
    return dirs