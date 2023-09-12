import logging
import os
import paramiko  # Utilitario para gerencia conexao SSH
import scp  # Utilitario para envio de arquivos com SCP
import re

def send_products(s_br, s_sp, dir_out):

    logging.info('')
    logging.info('ENVIANDO PRODUTOS PARA O SITE')
    try:
        # Criar objeto cliente SSH
        ssh = paramiko.SSHClient()
        # Carrega as chaves disponiveis no sistema
        ssh.load_system_host_keys() # find / -name known_hosts 2>/dev/null -----> acha aonde estão as chaves do sistema
        # Cria a conexao com o endereco informando, na porta informada e com o usuario informado
        ssh.connect('143.106.227.94', username="cpaunicamp", port=22000)
        # Cria o objeto cliente SCP com a conexao SSH
        scp_client = scp.SCPClient(ssh.get_transport())

        # Se a variavel de controle de processamento do brasil for True, realiza o processamento
        if s_br:
            logging.info('')
            logging.info('Enviando produtos "BR"')
            # Contador para envio nas 16 bandas
            for x in range(1, 17):
                # Transforma o inteiro contador em string e com 2 digitos
                b = str(x).zfill(2)
                # Cria uma lista com os itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_br.png$"
                ultima_br = [name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_br.png$', name)]
                # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
                ultima_br.sort()
                # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
                ultima_br.reverse()
                # Envia o arquivo "png" mais recente para o site, renomeando no destino
                scp_client.put(f'{dir_out}band{b}/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}_br.png')
                # Cria um arquivo menor do arquivo "png" mais recente
                os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}band{b}/{ultima_br[0]} -vf scale=448:321 {dir_out}band{b}/band{b}.png')
                # Envia o arquivo menor do "png" mais recente para o site
                scp_client.put(f'{dir_out}band{b}/band{b}.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}.png')
                # Envia o arquivo "gif" para o site
                scp_client.put(f'{dir_out}band{b}/band{b}_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}_br.gif')

        # Se a variavel de controle de processamento do estado de sao paulo for True, realiza o processamento
        if s_sp:
            logging.info('')
            logging.info('Enviando produtos "SP"')
            # Contador para envio nas 16 bandas
            for x in range(1, 17):
                # Transforma o inteiro contador em string e com 2 digitos
                b = str(x).zfill(2)
                # Cria uma lista com os itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_sp.png$"
                ultima_sp = [name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_sp.png$', name)]
                # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
                ultima_sp.sort()
                # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
                ultima_sp.reverse()
                # Envia o arquivo "png" mais recente para o site, renomeando no destino
                scp_client.put(f'{dir_out}band{b}/{ultima_sp[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}_sp.png')
                # Envia o arquivo "gif" para o site
                scp_client.put(f'{dir_out}band{b}/band{b}_sp.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}_sp.gif')


            # Envio TrueColor
            # Cria uma lista com os itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_br.png$"
            ultima_br = [name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_br.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_br.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}truecolor/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/truecolor/truecolor_br.png')
            # Cria um arquivo menor do arquivo "png" mais recente
            os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}truecolor/{ultima_br[0]} -vf scale=448:321 {dir_out}truecolor/truecolor.png')
            # Envia o arquivo menor do "png" mais recente para o site
            scp_client.put(f'{dir_out}truecolor/truecolor.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/truecolor/truecolor.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}truecolor/truecolor_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/truecolor/truecolor_br.gif')

            # Envio RRQPEF
            # Cria uma lista com os itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_br.png$"
            ultima_br = [name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_br.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_br.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}rrqpef/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/rrqpef/rrqpef_br.png')
            # Cria um arquivo menor do arquivo "png" mais recente
            os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}rrqpef/{ultima_br[0]} -vf scale=448:321 {dir_out}rrqpef/rrqpef.png')
            # Envia o arquivo menor do "png" mais recente para o site
            scp_client.put(f'{dir_out}rrqpef/rrqpef.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/rrqpef/rrqpef.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}rrqpef/rrqpef_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/rrqpef/rrqpef_br.gif')

            # Envio GLM
            # Cria uma lista com os itens no diretorio dos produtos glm que sao arquivos e se encaixa na expressao regular "^glm_.+_.+_br.png$"
            ultima_br = [name for name in os.listdir(f'{dir_out}glm') if os.path.isfile(os.path.join(f'{dir_out}glm', name)) and re.match('^glm_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_br.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_br.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}glm/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/glm/glm_br.png')
            # Cria um arquivo menor do arquivo "png" mais recente
            os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}glm/{ultima_br[0]} -vf scale=448:321 {dir_out}glm/glm.png')
            # Envia o arquivo menor do "png" mais recente para o site
            scp_client.put(f'{dir_out}glm/glm.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/glm/glm.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}glm/glm_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/glm/glm_br.gif')

            # Envio NDVI
            # Cria uma lista com os itens no diretorio dos produtos ndvi que sao arquivos e se encaixa na expressao regular "^ndvi_.+_.+_br.png$"
            ultima_br = [name for name in os.listdir(f'{dir_out}ndvi') if os.path.isfile(os.path.join(f'{dir_out}ndvi', name)) and re.match('^ndvi_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_br.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_br.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}ndvi/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/ndvi/ndvi_br.png')
            # Cria um arquivo menor do arquivo "png" mais recente
            os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}ndvi/{ultima_br[0]} -vf scale=448:321 {dir_out}ndvi/ndvi.png')
            # Envia o arquivo menor do "png" mais recente para o site
            scp_client.put(f'{dir_out}ndvi/ndvi.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/ndvi/ndvi.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}ndvi/ndvi_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/ndvi/ndvi_br.gif')

            # Envio FDCF
            # Cria uma lista com os itens no diretorio dos produtos fdcf que sao arquivos e se encaixa na expressao regular "^fdcf_.+_.+_br.png$"
            ultima_br = [name for name in os.listdir(f'{dir_out}fdcf') if os.path.isfile(os.path.join(f'{dir_out}fdcf', name)) and re.match('^fdcf_.+_.+_br.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_br.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_br.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}fdcf/{ultima_br[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/fdcf/fdcf_br.png')
            # Cria um arquivo menor do arquivo "png" mais recente
            os.system(f'/usr/bin/ffmpeg -y -v warning -i {dir_out}fdcf/{ultima_br[0]} -vf scale=448:321 {dir_out}fdcf/fdcf.png')
            # Envia o arquivo menor do "png" mais recente para o site
            scp_client.put(f'{dir_out}fdcf/fdcf.png', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/fdcf/fdcf.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}fdcf/fdcf_br.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/fdcf/fdcf_br.gif')

        # Se a variavel de controle de processamento do estado de sao paulo for True, realiza o processamento
        if s_sp:
            logging.info('')
            logging.info('Enviando produtos "SP"')
            # Contador para envio nas 16 bandas
            for x in range(1, 17):
                # Transforma o inteiro contador em string e com 2 digitos
                b = str(x).zfill(2)
                # Cria uma lista com os itens no diretorio dos produtos da banda que sao arquivos e se encaixa na expressao regular "^band[0-1][0-9]_.+_.+_sp.png$"
                ultima_sp = [name for name in os.listdir(f'{dir_out}band{b}') if os.path.isfile(os.path.join(f'{dir_out}band{b}', name)) and re.match('^band[0-1][0-9]_.+_.+_sp.png$', name)]
                # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
                ultima_sp.sort()
                # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
                ultima_sp.reverse()
                # Envia o arquivo "png" mais recente para o site, renomeando no destino
                scp_client.put(f'{dir_out}band{b}/{ultima_sp[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}_sp.png')
                # Envia o arquivo "gif" para o site
                scp_client.put(f'{dir_out}band{b}/band{b}_sp.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/band{b}/band{b}_sp.gif')

            # Envio TrueColor
            # Cria uma lista com os itens no diretorio dos produtos truecolor que sao arquivos e se encaixa na expressao regular "^truecolor_.+_.+_br.png$"
            ultima_sp = [name for name in os.listdir(f'{dir_out}truecolor') if os.path.isfile(os.path.join(f'{dir_out}truecolor', name)) and re.match('^truecolor_.+_.+_sp.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_sp.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_sp.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}truecolor/{ultima_sp[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/truecolor/truecolor_sp.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}truecolor/truecolor_sp.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/truecolor/truecolor_sp.gif')

            # Envio RRQPEF
            # Cria uma lista com os itens no diretorio dos produtos rrqpef que sao arquivos e se encaixa na expressao regular "^rrqpef_.+_.+_sp.png$"
            ultima_sp = [name for name in os.listdir(f'{dir_out}rrqpef') if os.path.isfile(os.path.join(f'{dir_out}rrqpef', name)) and re.match('^rrqpef_.+_.+_sp.png$', name)]
            # Ordena de forma alfabetica a lista, ficando assim os arquivos mais antigos no comeco
            ultima_sp.sort()
            # Realiza a inversao da lista, ficando assim os arquivos mais recentes no comeco
            ultima_sp.reverse()
            # Envia o arquivo "png" mais recente para o site, renomeando no destino
            scp_client.put(f'{dir_out}rrqpef/{ultima_sp[0]}', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/rrqpef/rrqpef_sp.png')
            # Envia o arquivo "gif" para o site
            scp_client.put(f'{dir_out}rrqpef/rrqpef_sp.gif', f'/var/www/html/cepagri/atualizacoes-regulares/goes16/rrqpef/rrqpef_sp.gif')

    except TimeoutError as e_timeout:
        logging.info('')
        logging.info(f'FALHA AO ENVIAR PRODUTOS PARA O SITE - {e_timeout}')
    except paramiko.SSHException as e_ssh:
        logging.info('')
        logging.info(f'FALHA AO ENVIAR PRODUTOS PARA O SITE - {e_ssh}')
    except scp.SCPException as e_scp:
        logging.info('')
        logging.info(f'FALHA AO ENVIAR PRODUTOS PARA O SITE - {e_scp}')
