bands = {}

# Todas as bandas da 01 a 21 recebem False
for num in range(1, 22):
    b = str(num).zfill(2)
    bands[f'{b}'] = True

br = True
sp = True

# Realiza etapas de processamento se houver alguma nova imagem
if bands["01"] or bands["02"] or bands["03"] or bands["04"] or bands["05"] or bands["06"] or bands["07"] or bands["08"] or \
        bands["09"] or bands["10"] or bands["11"] or bands["12"] or bands["13"] or bands["14"] or bands["15"] or bands["16"]:
    print(True)


if any(bands[key] for key in bands):
    print(False)