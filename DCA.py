import numpy as np
def score(a, b):
    return 0 if a == b else 1  

def affine_gap_penalty(length, open_cost=2, extend_cost=1):
    return open_cost + (length - 1) * extend_cost

def needleman_wunsch(s1, s2, gap_open=2, gap_extend=1):
    len_s1, len_s2 = len(s1), len(s2)
    score_matrix = np.zeros((len_s1 + 1, len_s2 + 1))
    
    for i in range(1, len_s1 + 1):
        score_matrix[i][0] = affine_gap_penalty(i, gap_open, gap_extend)
    for j in range(1, len_s2 + 1):
        score_matrix[0][j] = affine_gap_penalty(j, gap_open, gap_extend)
    
    for i in range(1, len_s1 + 1):
        for j in range(1, len_s2 + 1):
            match = score_matrix[i-1][j-1] + score(s1[i-1], s2[j-1])
            delete = score_matrix[i-1][j] + affine_gap_penalty(1, gap_open, gap_extend)
            insert = score_matrix[i][j-1] + affine_gap_penalty(1, gap_open, gap_extend)
            score_matrix[i][j] = min(match, delete, insert)
    
    return score_matrix

def pad_matrix(matrix, target_size, pad_value=10):
    current_size = matrix.shape[0]
    if current_size < target_size:
        padding = np.full((target_size - current_size, matrix.shape[1]), pad_value)
        matrix = np.vstack((matrix, padding))
    return matrix


def dca(s1, s2, gap_open=2, gap_extend=1, recursion_stop=5):

    if len(s1) <= recursion_stop or len(s2) <= recursion_stop:
        return needleman_wunsch(s1, s2, gap_open, gap_extend)

    mid_s1 = len(s1) // 2
    left_half = dca(s1[:mid_s1], s2, gap_open, gap_extend, recursion_stop)
    right_half = dca(s1[mid_s1:], s2, gap_open, gap_extend, recursion_stop)
    max_rows = max(left_half.shape[0], right_half.shape[0])
    left_half = pad_matrix(left_half, max_rows)
    right_half = pad_matrix(right_half, max_rows)
    
    return np.concatenate((left_half, right_half), axis=1)

s1 = "GATTACA"
s2 = "GCATGCU"

alignment_matrix = dca(s1, s2)
print(alignment_matrix)
