__author__ = 'david'

from scipy.special import hyp2f1, binom
import numpy as np
from scipy.optimize import brentq


def a_2m(m, z):
    fac1 = (1 + z) ** (2 * m)
    return fac1 - fac1 * (z ** m) * (m + 1) / m * binom(2 * m, m - 1) * hyp2f1(2 * m + 1, m, m + 1, -z)


def a_2m_der(m, z):
    return 2 * m * a_2m(m, z) / (1 + z) - (m + 1) * binom(2 * m, m - 1) * (z ** (m - 1)) / (1 + z)


def a_2m_1(m, z):
    fac1 = (1 + z) ** (2 * m + 1)
    return fac1 - binom(2 * m + 1, m) * fac1 * (z ** (m + 1)) * hyp2f1(2 * m + 2, m + 1, m + 2, -z)


def a_2m_1_der(m, z):
    return (2 * m + 1) * a_2m_1(m, z) / (1 + z) - (m + 1) * binom(2 * m + 1, m) * (z ** m) / (1 + z)


def init_all_exp_2m(temp, alpha, m):
    exp_1 = np.exp(-(alpha - 1) / temp)
    exp_1_2m = np.exp((alpha - 1) / temp / m)
    exp_2 = np.exp(alpha / temp)
    exp_2_2m = np.exp(-alpha / temp / m)
    return exp_1, exp_1_2m, exp_2, exp_2_2m


def init_all_exp_2m_1(temp, alpha, m):
    exp_1 = np.exp(-(alpha - 1) / temp)
    exp_1_2m_1 = np.exp(2 * (alpha - 1) / temp / (2 * m + 1))
    exp_2 = np.exp(alpha / temp)
    exp_2_2m_1 = np.exp(-2 * alpha / temp / (2 * m + 1))
    return exp_1, exp_1_2m_1, exp_2, exp_2_2m_1


def z_2m(m, exp_1, exp_1_2m, exp_2, exp_2_2m):
    return exp_1 * a_2m(m, exp_1_2m) + binom(2 * m, m) + exp_2 * a_2m(m, exp_2_2m)


def z_2m_1(m, exp_1, exp_1_2m_1, exp_2, exp_2_2m_1):
    return exp_1 * a_2m_1(m, exp_1_2m_1) + exp_2 * a_2m_1(m, exp_2_2m_1)


def q_2m(m, exp_1, exp_1_2m, exp_2, exp_2_2m):
    return exp_1 * exp_1_2m * a_2m_der(m, exp_1_2m) + m * binom(2 * m, m) + \
           exp_2 * 2 * m * a_2m(m, exp_2_2m) - exp_2 * exp_2_2m * a_2m_der(m, exp_2_2m)


def q_2m_1(m, exp_1, exp_1_2m_1, exp_2, exp_2_2m_1):
    return exp_1 * exp_1_2m_1 * a_2m_1_der(m, exp_1_2m_1) + \
           exp_2 * (2 * m + 1) * a_2m_1(m, exp_2_2m_1) - exp_2 * exp_2_2m_1 * a_2m_1_der(m, exp_2_2m_1)


def func_2m(alpha, temp, m, mean_gamma):
    exp_1, exp_1_2m, exp_2, exp_2_2m = init_all_exp_2m(temp, alpha, m)
    return q_2m(m, exp_1, exp_1_2m, exp_2, exp_2_2m) / z_2m(m, exp_1, exp_1_2m, exp_2, exp_2_2m) - \
           m * (1 + mean_gamma) / mean_gamma


def func_2m_1(alpha, temp, m, mean_gamma):
    exp_1, exp_1_2m_1, exp_2, exp_2_2m_1 = init_all_exp_2m_1(temp, alpha, m)
    return q_2m_1(m, exp_1, exp_1_2m_1, exp_2, exp_2_2m_1) / z_2m_1(m, exp_1, exp_1_2m_1, exp_2, exp_2_2m_1) - \
           (2 * m + 1) * (1 + mean_gamma) / 2 / mean_gamma


def find_root(temp, alpha0=0, alpha1=1, g=2, mean_gamma=2):
    if g % 2 == 0:
        return brentq(lambda x: func_2m(x, temp, int(g / 2), mean_gamma), alpha0, alpha1)
    else:
        return brentq(lambda x: func_2m_1(x, temp, int((g - 1) / 2), mean_gamma), alpha0, alpha1)


def construct_phase_diagram(temp_list, g, mean_gamma, fileout):
    w = open(fileout, "w")
    for temp in temp_list:
        try:
            alpha = find_root(temp, g=g, mean_gamma=mean_gamma)
            w.write(str(temp) + "\t" + str(alpha) + "\n")
            print(temp, alpha)
        except ValueError:
            print("There is no transition for T > " + str(temp))
            break
    w.close()


def get_temp_list(temp_f, dtemp):
    temp_list = np.arange(0, temp_f, dtemp)
    temp_list[0] = 0.005
    return temp_list


def main():
    g = 5
    c = 5
    temp_f = 2
    dtemp = 0.01

    mean_gamma = c - 1
    graph_type = 'RGfc'
    fileout = "Phase_Diagram_Homophily_" + graph_type + "_c_" + str(c) + "_G_" + str(g) + ".txt"
    temp_list = get_temp_list(temp_f, dtemp)

    construct_phase_diagram(temp_list, g, mean_gamma, fileout)

    return 0


if __name__ == '__main__':
    main()
