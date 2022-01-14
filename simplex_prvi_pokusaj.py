# from scipy.optimize import linprog 

# c = [-1, -1, -1, -1]
# A_ub = [[-1, -1, -1, -1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]
# b_ub = [-1, 0, 0, 0]
# A_eq = [[1, 1, 1, 1]]
# b_eq = [1]
# res = linprog(c, A_ub, b_ub, A_eq, b_eq, method='simplex')
# print(res)



# linprog(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None, method='simplex', callback=None, options={'maxiter': 5000, 'disp': False, 'presolve': True, 'tol': 1e-12, 'autoscale': False, 'rr': True, 'bland': False}, x0=None)

# minimize: min(-2x_1 - 3x_2 - 4x_3)
# input c = [-2, -3, -4]^T
# input c = [-2, -3, -4, 0, 0]
# with: 
# 3x_1 + 2x_2 + 1x_3 <= 10
# 2x_1 + 5x_2 + 3x_3 <= 15
# x_1, x_2, x_3 >= 0
# which then becomes: 
# 3x_1 + 2x_2 + 1x_3 + x_4 = 10
# 2x_1 + 5x_2 + 3x_3 + x_5 = 15
# meaning A = [[3, 2, 1, 1, 0], [2, 5, 3, 0, 1]]
# and b = [10, 15]



from cmath import pi
from re import T
import traceback
import numpy as np

def simplex(c, A, b):



    # dodati kolone, ukoliko je to potrebno
    # broj ogranicenja
    n = A.shape[0]
    # tu gde je 1 u base-u ce znaciti da je to bazna promenljiva
    base = np.zeros(c.shape[0], dtype=int)
    A_1 = np.identity(n)
    
    for i in range(n):
        c = np.append(c, 0)    
        base = np.append(base, 1)
    A = np.append(A, A_1, axis=1)

    # broj promenljivih
    m = A.shape[1]

    B = A_1
    B_inv = np.linalg.inv(B)

    table = np.zeros((n + 1, n + 3))

    free_column = table[:, n + 2]
    bool_free_column = False
    for i in free_column:
        if i < 0:
            bool_free_column = True
            break

    bool_free_column = True
    while (bool_free_column): 
        bool_free_column = False
        j = 0
        k = 0

        # u prvu kolonu postavlja koja je zapravo po redu promenljiva
        # TODO: stavlja poziciju + 1
        # ako se j += 1 stavi ispod if-a onda ce staviti na poziciju, ali videcu jos kako cu uraditi
        for i in base:
            j += 1
            if i == 1:
                k += 1
                table[k: k + 1, 0] = j

    table[1:, 1:-2] = B_inv

    # omega
    c_B = np.array([])
    for i in range(m):
        if base[i] == 1:
            c_B = np.append(c_B, c[i])
    w = c_B.dot(B_inv)
    table[0, 1:-2] = w

    b_B = B_inv @ b

    table[1:, -2] = b_B.transpose()

    c_Bb_B = c_B.dot(b_B)
    table[0, -2] = c_Bb_B

    # msm da bi ovde trebao for

    adaaada = 0

    while adaaada < 10:
        print(adaaada)
        adaaada += 1

        finding_min = np.array([])
        for j in range(1, m + 1):
            if j not in table[1:, 0:1]:
                temp = w.dot(A[:, j - 1]) - c[j - 1]
                finding_min = np.append(finding_min, temp)


        min_index = 0
        min_elem = finding_min[0]
        for (index, i) in enumerate(finding_min): 
            if i < min_elem:
                min_elem = i
                min_index = index

        if finding_min[min_index] > 0:
            # TODO: ovde je return 
            # prva vrednost je vrednost max z
            return (table[0, -2])
        



        # TODO: ovde treba + 1 eventualno

        table[0, -1] = min_elem

        
        

        pivot_row_index = 0
        pivot_column_index = 0
        pivot_value = 0

        pivot_column_index = n + 2
        pivot_row_index = np.argmin(table[1:, pivot_column_index]) + 1
        print("+" * 20)
        print(table)
        print("+" * 20)
        print(table[pivot_row_index, pivot_column_index])
        if (table[pivot_row_index, pivot_column_index] < 0):
            print("error")
            return 
        pivot_value = table[pivot_row_index, pivot_column_index]
        

        for i in range(table.shape[0]):
            if i != pivot_row_index:
                for k in range (table.shape[1] - 2):
                    table[i, k + 1] = table[i, k + 1] - table[i, pivot_column_index] * table[pivot_row_index, k + 1] / pivot_value

        table[pivot_row_index, 1:-1] = table[pivot_row_index, 1:-1] / pivot_value
        table[pivot_row_index, 0] = min_index + 1
        print(f"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa {table[pivot_row_index, 0]}")
        
        # a_j = A[:, min_index]
        to_replace = table[pivot_row_index, 0]

        # TODO: ovde stao 

        print(int(to_replace - 1))
        B_inv = table[1:, 1:-2]
        print(A)
        print(A[:, int(to_replace - 1)])
        a_j = A[:, int(to_replace - 1)]
        print(a_j)
        table[1:, -1] = a_j @ B_inv
        print("6" * 100)
        print(a_j)
        print(B_inv)
        print(a_j @ B_inv)
        print("6" * 100)
        print(table)
        print("+" * 20)
        
        

        base[min_index] = 1
        base[int(to_replace - 1)] = 0
        print("-" * 50)
        print(min_index)
        print(pivot_row_index)
        print(base)
        print("-" * 50)

       
        

        print(table)

        c_B = np.array([])
        for i in range(m):
            if base[i] == 1:
                c_B = np.append(c_B, c[i])
        print(c_B)
        print("aaaaaaaaaaaaaa")
        w = c_B.dot(B_inv)
        print(w)
        table[0, 1:-2] = w
        print(table)

        b_B = B_inv @ b
        print(b_B)

        table[1:, -2] = b_B.transpose()
        print(table)

        c_Bb_B = c_B.dot(b_B)
        print(c_Bb_B)
        table[0, -2] = c_Bb_B

        print(table)






    
c = np.array([6, 14, 13])
A = np.array([[0.5, 2, 1], [1, 2, 4]])
b = np.array([[24], [60]])
print(simplex(c = c, A = A, b = b))

# c = np.array([6, 14, 13, 14])
# A = np.array([[0.5, 2, 1, 4], [1, 2, 4, 3], [3, 3, 4, 6]])
# b = np.array([[24], [60], [13]])
# simplex(c = c, A = A, b = b)

# c = np.array([6, 14, 13, 14, 17])
# A = np.array([[0.5, 2, 1, 4, 6], [1, 2, 4, 3, 6], [3, 3, 4, 6, 1], [1, 1, 1, 6, 6]])
# b = np.array([[24], [60], [13], [5]])
# simplex(c = c, A = A, b = b)