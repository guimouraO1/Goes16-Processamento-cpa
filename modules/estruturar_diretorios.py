from dirs import get_dirs
import os

dirs = get_dirs()

# Acessando os diretórios usando as chaves do dicionário
dir_main = dirs['dir_main']

def criar_estrutura(dir_main):
    if not os.path.exists(f'{dir_main}logs'):
        # Cria o diretório base se não existir
        os.makedirs(f'{dir_main}logs')
        print(f'Diretório ''logs'' criado com sucesso.')
    else:
        print(f'Diretório logs já existe.')

    # Lista dos diretórios
    diretorios_base = ['goes', 'output']

    for diretorio_base in diretorios_base:
        diretorio_base_path = os.path.join(dir_main, diretorio_base)

        # Verifica se o diretório base já existe
        if not os.path.exists(diretorio_base_path):
            # Cria o diretório base se não existir
            os.makedirs(diretorio_base_path)
            print(f'Diretório "{diretorio_base_path}" criado com sucesso.')
        else:
            print(f'Diretório "{diretorio_base_path}" já existe.')

        # Subdiretórios band01 até band16 dentro do diretório base
        for i in range(1, 17):
            band_dir = os.path.join(diretorio_base_path, f"band{i:02d}")
            # Cria o subdiretório bandXX dentro do diretório base (se não existir)
            if not os.path.exists(band_dir):
                os.makedirs(band_dir)
                
        # Subdiretórios fdcf, rrqpef, glm, ndvi, truecolor dentro de cada bandXX
        subdirs = ["fdcf", "rrqpef", "glm", "ndvi", "truecolor", "clsm"]
        for subdir in subdirs:
            subdir_path = os.path.join(diretorio_base_path, subdir)
            # Cria o subdiretório (se não existir)
            if not os.path.exists(subdir_path):
                os.makedirs(subdir_path)
                print(f'Diretório "{subdir_path}" criado com sucesso.')


# Função para criar a estrutura
criar_estrutura(dir_main)
