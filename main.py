import numpy as np


def print_hi(name):
    text = 'Hi, %s!' % name

    return text


def avg_calc(*args):
    stat = (sum(list(args)))/len(list(args))

    return stat


if __name__ == '__main__':
    print(print_hi('Alice'))
    print(avg_calc(7, 8, 6))

