import numpy as np

def calcular_penalidade_gap(k):
    return 2 + (k - 1)

def alinhar(sequencias):
    if len(sequencias) == 1:
        return sequencias[0]
    
    mid = len(sequencias) // 2
    left_align = alinhar(sequencias[:mid])
    right_align = alinhar(sequencias[mid:])
    
    return alinhar_duas(left_align, right_align)

def alinhar_duas(seq1, seq2):
    len1, len2 = len(seq1), len(seq2)
    matriz = np.zeros((len1 + 1, len2 + 1))
    
    for i in range(1, len1 + 1):
        matriz[i][0] = calcular_penalidade_gap(i)
    for j in range(1, len2 + 1):
        matriz[0][j] = calcular_penalidade_gap(j)
    
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            match = matriz[i-1][j-1] + (0 if seq1[i-1] == seq2[j-1] else 1)
            delete = matriz[i-1][j] + calcular_penalidade_gap(1)
            insert = matriz[i][j-1] + calcular_penalidade_gap(1)
            matriz[i][j] = min(match, delete, insert)
    
    # Imprimir a matriz de pontuação
    print("Matriz de Pontuação:")
    print(matriz)
    
    alinhamento1 = []
    alinhamento2 = []
    i, j = len1, len2
    while i > 0 and j > 0:
        score = matriz[i][j]
        score_diag = matriz[i-1][j-1]
        score_up = matriz[i][j-1]
        score_left = matriz[i-1][j]
        
        if score == score_diag + (0 if seq1[i-1] == seq2[j-1] else 1):
            alinhamento1.append(seq1[i-1])
            alinhamento2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif score == score_left + calcular_penalidade_gap(1):
            alinhamento1.append(seq1[i-1])
            alinhamento2.append('-')
            i -= 1
        elif score == score_up + calcular_penalidade_gap(1):
            alinhamento1.append('-')
            alinhamento2.append(seq2[j-1])
            j -= 1
    
    while i > 0:
        alinhamento1.append(seq1[i-1])
        alinhamento2.append('-')
        i -= 1
    while j > 0:
        alinhamento1.append('-')
        alinhamento2.append(seq2[j-1])
        j -= 1
    
    return ''.join(reversed(alinhamento1)), ''.join(reversed(alinhamento2))

# Exemplo de uso
sequencias = ["GATTACA", "GCATGCU"]
alinhamento = alinhar(sequencias)
print(f"Alinhamento resultante: {alinhamento}")