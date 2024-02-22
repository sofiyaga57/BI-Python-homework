import numpy as np


class User:
    pass


class Booking:
    pass


class LabEquipment:
    pass


class GenCodeInterpreter:
    pass


def meet_the_dunders():
    res = 0

    matrix = []
    for idx in range(0, 100, 10):
        matrix += [list(range(idx, idx + 10))]

    def func_1(x):
        return x in range(1, 5, 2)
        
    def func_2(x):
        return [x[col] for col in selected_columns_indices]
        
    selected_columns_indices = list(filter(func_1, range(len(matrix))))
    selected_columns = map(func_2, matrix)

    arr = np.array(list(selected_columns))

    mask = arr[:, 1] % 3 == 0
    new_arr = arr[mask]

    product = new_arr @ new_arr.T

    if (product[0] < 1000).all() and (product[2] > 1000).any():
        res = int(product.mean() // 10 % 100)
    return res