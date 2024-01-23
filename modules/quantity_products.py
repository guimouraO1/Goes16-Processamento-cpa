import os
import re
import logging

# Esta função remove produtos em excesso com base no tipo de produto e no diretório de saída.
def remove_excess_products(dir_out, product_type):
    try:
        # Conta quantos arquivos correspondem ao padrão especificado para o Brasil (br.png) e São Paulo (sp.png).
        result_br = len([name for name in os.listdir(f'{dir_out}{product_type}') if os.path.isfile(os.path.join(f'{dir_out}{product_type}', name)) and re.match(f'^{product_type}_.+_.+_br.png$', name)])
        result_sp = len([name for name in os.listdir(f'{dir_out}{product_type}') if os.path.isfile(os.path.join(f'{dir_out}{product_type}', name)) and re.match(f'^{product_type}_.+_.+_sp.png$', name)])
        
        # Subtrai 48 do resultado para determinar a quantidade máxima permitida.
        result_br -= 48
        result_sp -= 48

        # Se houver produtos em excesso, eles serão removidos.
        if result_br > 0 or result_sp > 0:
            prod_br = [name for name in os.listdir(f'{dir_out}{product_type}') if os.path.isfile(os.path.join(f'{dir_out}{product_type}', name)) and re.match(f'^{product_type}_.+_.+_br.png$', name)]
            prod_sp = [name for name in os.listdir(f'{dir_out}{product_type}') if os.path.isfile(os.path.join(f'{dir_out}{product_type}', name)) and re.match(f'^{product_type}_.+_.+_sp.png$', name)]
            prod_br.sort()
            prod_sp.sort()

            # Remove os produtos em excesso para o Brasil.
            for y in range(0, result_br):
                os.remove(f'{dir_out}{product_type}/{prod_br[y]}')

            # Remove os produtos em excesso para São Paulo.
            for y in range(0, result_sp):
                os.remove(f'{dir_out}{product_type}/{prod_sp[y]}')

            return True

        return False
    except Exception as e:
        logging.error(f"Erro ao processar {product_type}: {str(e)}")
        return False

# Esta função verifica a quantidade de produtos em excesso em várias categorias.
def quantity_products(dir_out):
    logging.info("")
    logging.info("VERIFICANDO NUMERO DE PRODUTOS...")
    aux = False
    try:
        # Loop para verificar várias categorias de produtos.
        for x in range(1, 17):
            b = str(x).zfill(2)
            if remove_excess_products(dir_out, f'band{b}'):
                aux = True

        # Verifica e remove produtos em excesso para outras categorias.
        if remove_excess_products(dir_out, 'truecolor'):
            aux = True

        if remove_excess_products(dir_out, 'rrqpef'):
            aux = True

        if remove_excess_products(dir_out, 'glm'):
            aux = True

        if remove_excess_products(dir_out, 'ndvi'):
            aux = True

        if remove_excess_products(dir_out, 'fdcf'):
            aux = True
        
        if remove_excess_products(dir_out, 'airmass'):
            aux = True
            
        if remove_excess_products(dir_out, 'lst'):
            aux = True

        if remove_excess_products(dir_out, 'dmw'):
            aux = True

        # Registra se produtos em excesso foram removidos com sucesso ou não.
        if aux:
            logging.info("Produtos em excesso removidos com sucesso")
        else:
            logging.info("Quantidade de produtos dentro do limite")
    except Exception as e:
        logging.info(f'Erro {e} ao remover imagens')