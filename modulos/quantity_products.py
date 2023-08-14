import logging
import os
import re

def quantity_products(dir_out):

    logging.info("")
    logging.info("VERIFICANDO NUMERO DE PRODUTOS...")
    aux = False

    # Contador para verificar produtos nas 16 bandas
    for x in range(1, 17):
        # Transforma o inteiro contador em string e com 2 digitos
        b = str(x).zfill(2)
        # Verifica a quantidade de itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_br.png$"
        result_br = len([name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_br.png$', name)])
        # Verifica a quantidade de itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_sp.png$"
        result_sp = len([name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_sp.png$', name)])
        # Subtrai o limite de imagens BR
        result_br -= 48
        # Subtrai o limite de imagens SP
        result_sp -= 48
        # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
        if result_br > 0:
            aux = True
            # Cria uma lista com os itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_br.png$"
            prod_br = [name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            prod_br.sort()
            # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
            for y in range(0, result_br):
                # Remove arquivo em excesso
                os.remove(f'{dir_out}band{b}/{prod_br[y]}')
        # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
        if result_sp > 0:
            aux = True
            # Cria uma lista com os itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_sp.png$"
            prod_sp = [name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_sp.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            prod_sp.sort()
            # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
            for y in range(0, result_sp):
                # Remove arquivo em excesso
                os.remove(f'{dir_out}band{b}/{prod_sp[y]}')

    # Verificar produtos truecolor
    # Verifica a quantidade de itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_br.png$"
    result_br = len([name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_br.png$', name)])
    # Verifica a quantidade de itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_sp.png$"
    result_sp = len([name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_sp.png$', name)])
    # Subtrai o limite de imagens BR
    result_br -= 48
    # Subtrai o limite de imagens SP
    result_sp -= 48
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_br > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_br.png$"
        prod_br = [name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_br.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_br.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_br):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}truecolor/{prod_br[y]}')
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_sp > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_sp.png$"
        prod_sp = [name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_sp.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_sp.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_sp):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}truecolor/{prod_sp[y]}')

    # Verificar produtos rrqpef
    # Verifica a quantidade de itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_br.png$"
    result_br = len([name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_br.png$', name)])
    # Verifica a quantidade de itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_sp.png$"
    result_sp = len([name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_sp.png$', name)])
    # Subtrai o limite de imagens BR
    result_br -= 48
    # Subtrai o limite de imagens SP
    result_sp -= 48
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_br > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_br.png$"
        prod_br = [name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_br.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_br.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_br):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}rrqpef/{prod_br[y]}')
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_sp > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_sp.png$"
        prod_sp = [name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_sp.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_sp.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_sp):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}rrqpef/{prod_sp[y]}')

    # Verificar produtos glm
    # Verifica a quantidade de itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_br.png$"
    result_br = len([name for name in os.listdir(f'{dir_out}glm') if os.path.isfile(os.path.join(f'{dir_out}glm', name)) and re.match('^glm_.+_.+_br.png$', name)])
    # Verifica a quantidade de itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_sp.png$"
    result_sp = len([name for name in os.listdir(f'{dir_out}glm') if os.path.isfile(os.path.join(f'{dir_out}glm', name)) and re.match('^glm_.+_.+_sp.png$', name)])
    # Subtrai o limite de imagens BR
    result_br -= 48
    # Subtrai o limite de imagens SP
    result_sp -= 48
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_br > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_br.png$"
        prod_br = [name for name in os.listdir(f'{dir_out}glm') if os.path.isfile(os.path.join(f'{dir_out}glm', name)) and re.match('^glm_.+_.+_br.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_br.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_br):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}glm/{prod_br[y]}')
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_sp > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_sp.png$"
        prod_sp = [name for name in os.listdir(f'{dir_out}glm') if os.path.isfile(os.path.join(f'{dir_out}glm', name)) and re.match('^glm_.+_.+_sp.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_sp.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_sp):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}glm/{prod_sp[y]}')

    # Verificar produtos ndvi
    # Verifica a quantidade de itens no diretorio dos produtos ndvi que sao arquivos e se encaixa na expressao regular "^ndvi_.+_.+_br.png$"
    result_br = len([name for name in os.listdir(f'{dir_out}ndvi') if os.path.isfile(os.path.join(f'{dir_out}ndvi', name)) and re.match('^ndvi_.+_.+_br.png$', name)])
    # Subtrai o limite de imagens BR
    result_br -= 48
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_br > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos ndvi que sao arquivos e se encaixa na expressao regular "^ndvi_.+_.+_br.png$"
        prod_br = [name for name in os.listdir(f'{dir_out}ndvi') if os.path.isfile(os.path.join(f'{dir_out}ndvi', name)) and re.match('^ndvi_.+_.+_br.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_br.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_br):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}ndvi/{prod_br[y]}')

    # Verificar produtos fdcf
    # Verifica a quantidade de itens no diretorio dos produtos fdcf que sao arquivos e se encaixa na expressao regular "^fdcf_.+_.+_br.png$"
    result_br = len([name for name in os.listdir(f'{dir_out}fdcf') if os.path.isfile(os.path.join(f'{dir_out}fdcf', name)) and re.match('^fdcf_.+_.+_br.png$', name)])
    # Subtrai o limite de imagens BR
    result_br -= 48
    # Se o resultado for maior que zero, temos imagens em excesso, entao serao removidas
    if result_br > 0:
        aux = True
        # Cria uma lista com os itens no diretorio dos produtos fdcf que sao arquivos e se encaixa na expressao regular "^fdcf_.+_.+_br.png$"
        prod_br = [name for name in os.listdir(f'{dir_out}fdcf') if os.path.isfile(os.path.join(f'{dir_out}fdcf', name)) and re.match('^fdcf_.+_.+_br.png$', name)]
        # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
        prod_br.sort()
        # Looping do tamanho dos arquivos em excesso para remocao, removendo do inicio da lista os arquivos em excesso
        for y in range(0, result_br):
            # Remove arquivo em excesso
            os.remove(f'{dir_out}fdcf/{prod_br[y]}')

    if aux:
        logging.info("Produtos em excesso removidos com sucesso")
    else:
        logging.info("Quantidade de produtos dentro do limite")
