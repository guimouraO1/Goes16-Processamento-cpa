"""
    Retorna uma matriz RGBA normalizada correspondente a *x*.

    No caso normal, *x* é uma sequência 1D ou 2D de escalares, e
    a correspondente matriz `~numpy.ndarray` de valores RGBA será retornada,
    com base na normalização e mapa de cores definido para este objeto ScalarMappable.

    Existe um caso especial, para lidar com imagens que já são
    RGB ou RGBA, como poderia ter sido lido de um arquivo de imagem.
    Se *x* for uma `~numpy.ndarray` com 3 dimensões,
    e a última dimensão for 3 ou 4, então ela será
    tratada como uma matriz RGB ou RGBA, e nenhuma conversão será feita.
    A matriz pode ser do tipo `~numpy.uint8`, ou pode ser números de ponto flutuante com
    valores na faixa de 0-1; caso contrário, uma ValueError será levantada.
    Se for uma matriz mascarada, quaisquer elementos mascarados serão definidos como alfa 0.
    Se a última dimensão for 3, o argumento *alpha* (com valor padrão de 1)
    será usado para preencher a transparência. Se a última dimensão
    for 4, o argumento *alpha* é ignorado; ele não substitui o alfa preexistente.
    Uma ValueError será levantada
    se a terceira dimensão for diferente de 3 ou 4.

    Em ambos os casos, se *bytes* for *False* (padrão), a matriz RGBA
    será de números de ponto flutuante na faixa de 0-1; se for *True*,
    a matriz RGBA retornada será do tipo `~numpy.uint8` na faixa de 0 a 255.

    Se o argumento *norm* for False, nenhuma normalização dos dados de entrada será
    realizada e será assumido que eles estão na faixa (0-1).

"""
