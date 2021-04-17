import numpy as np


def print_hi(name, age):
    text = 'Hi, %s!' % name
    Age = '%f' % age

    return text


def avg_calc(*args):
    stat = (sum(list(args)))/len(list(args))

    return stat


if __name__ == '__main__':
    print(print_hi('Alice', 22))
    print(avg_calc(7, 8, 6))


