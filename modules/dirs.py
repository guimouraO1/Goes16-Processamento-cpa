from datetime import date

def get_dirs():
    dir_main = '/home/guimoura/processamento/'
    #dir_main = 'C:/Users/Gui/OneDrive/Documentos/processamentopy/'
    dirs = {
        'dir_main': dir_main,
        'dir_in': f'{dir_main}goes/',
        'dir_maps': f'{dir_main}maps/',
        'dir_out': f'{dir_main}output/',
        'dir_shapefiles': f'{dir_main}shapefiles/',
        'dir_colortables': f'{dir_main}colortables/',
        'dir_logos': f'{dir_main}logos/',
        'dir_temp': f'{dir_main}temp/',
        'arq_log': f'{dir_main}logs/Processamento-GOES_{str(date.today())}.log'
    }
    
    return dirs