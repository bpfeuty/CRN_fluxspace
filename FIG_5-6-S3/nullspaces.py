import numpy as np
import sympy as sp

def find_integer_right_null_vectors(matrix):
    sym_matrix = sp.Matrix(matrix)
    null_space = sym_matrix.nullspace() 
    #SymPy often uses the RREF to determine the null space of a matrix
    if not null_space:
        raise ValueError("The matrix has a trivial null space (only the zero vector).")
 
    # Ensure eigenvectors have integer values by finding a common multiple
    integer_right_null_vectors = []
    for vec in null_space:
        # Find the least common multiple (LCM) of the denominators to convert to integers
        lcm_denominator = sp.lcm([term.q for term in vec])
        integer_vec = [int(term * lcm_denominator) for term in vec]
        integer_right_null_vectors.append(integer_vec)

    return integer_right_null_vectors

def right_null_space_dimension(matrix):
    # Convert the input matrix to a SymPy Matrix
    sym_matrix = sp.Matrix(matrix)
    
    # Find the null space of the matrix (right null space)
    null_space = sym_matrix.nullspace()

    # The dimension of the null space is the number of vectors in it
    return len(null_space)

def find_integer_left_null_vectors(matrix):
    # Convert the input matrix to a SymPy Matrix
    sym_matrix = sp.Matrix(matrix)
    
    # Compute the transpose of the matrix
    sym_matrix_T = sym_matrix.transpose()
    
    # Find the null space of the transpose matrix (left null space of the original matrix)
    null_space_T = sym_matrix_T.nullspace()

    if not null_space_T:
        raise ValueError("The matrix has a trivial left null space (only the zero vector).")
    
    # Ensure null vectors have integer values by finding a common multiple
    integer_left_null_vectors = []
    for vec in null_space_T:
        # Find the least common multiple (LCM) of the denominators to convert to integers
        lcm_denominator = sp.lcm([term.q for term in vec])
        integer_vec = [int(term * lcm_denominator) for term in vec]
        integer_left_null_vectors.append(integer_vec)

    return integer_left_null_vectors

def left_null_space_dimension(matrix):
    # Convert the input matrix to a SymPy Matrix
    sym_matrix = sp.Matrix(matrix)
    
    # Compute the transpose of the matrix
    sym_matrix_T = sym_matrix.transpose()
    
    # Find the null space of the transpose matrix (left null space of the original matrix)
    null_space_T = sym_matrix_T.nullspace()

    # The dimension of the null space is the number of vectors in it
    return len(null_space_T)

# Example usage
matrix = [
    [1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 2, 0],
    [0, 0, 0, 0, 1, 0, 0,10, 0,-10,13,-1],
    [0, 0, 0, 0, 0, 1, 0, 0, 0,  0,-3, 0],
    [0, 0, 0, 0, 0, 0, 1, 8, 1, 0, -4, 0],
    [0, 0, 0, 0, 0, 0, 0, 2,-1, -1,-1, 0],
    [0, 0, 0, 0, 0, 0, 0,-10,-1,10, -10, 1],
    [0, 0, 0, 0, 0, 0, 0,10, 1,-10, 10, -1],
    [-1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0],
    [0, -1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0],
    [0, 0, -1, 0, 0, 0,  0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1, 0, 0,  0, 0, 0, 0, 0, 0],
    [0, 0, 0,  0, -1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0,  0, 0, -1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0,  0, 0,  0, -1, 0, 0, 0, 0, 0],
]

matrix2 = [
    [1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 2, 0],
    [0, 0, 0, 0, 1, 0, 0,10, 0,-10,13,-1],
    [0, 0, 0, 0, 0, 1, 0, 0, 0,  0,-3, 0],
    [0, 0, 0, 0, 0, 0, 1, 8, 1, 0, -4, 0],
    [0, 0, 0, 0, 0, 0, 0, 2,-1, -1,-1, 0],
    [0, 0, 0, 0, 0, 0, 0,-10,-1,10, -10, 1],
    [0, 0, 0, 0, 0, 0, 0,10,1,-10, 10, -1],
]


dimcoker = left_null_space_dimension(matrix)
dimcokerx = left_null_space_dimension(matrix2)
print("dimensionality of the left null space of S:", dimcoker)
print("dimensionality of the left null space of S^x (n_uc):", dimcokerx)
print("number of broken conservation laws:", dimcoker-dimcokerx)
dimkerx = right_null_space_dimension(matrix2)
dimker = right_null_space_dimension(matrix)
print("dimensionality of the right null space of S (n_cc):", dimker)
print("dimensionality of the right null space of Sx (n_ec+n_cc):", dimkerx)
print("number of emergent cycles:", dimkerx-dimker)
left_null_vectors = find_integer_left_null_vectors(matrix)
right_null_vectors = find_integer_right_null_vectors(matrix2)
print("conservation vectors:")
for vec in left_null_vectors:
    print(vec,np.dot(vec, matrix))

print("reaction cycle vectors:")
for vec in right_null_vectors:
    print(vec,np.dot(matrix2, vec))
