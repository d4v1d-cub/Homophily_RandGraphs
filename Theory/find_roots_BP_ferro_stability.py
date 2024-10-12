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


def init_all_exp_2m(temp, alpha, m, gamma):
    exp_1 = np.exp(-(alpha - 1) * gamma / temp)
    exp_1_2m = np.exp((alpha - 1) / temp / m)
    exp_1_2m_neg = 1.0 / exp_1_2m
    exp_1_2m_gamma = np.exp((alpha - 1) * gamma / temp / m)
    sinh_1_2m = (exp_1_2m - exp_1_2m_neg) / 2
    exp_2 = np.exp(alpha * gamma / temp)
    exp_2_2m = np.exp(-alpha / temp / m)
    exp_2_2m_pos = 1.0 / exp_2_2m
    exp_2_2m_gamma = np.exp(-alpha * gamma / temp / m)
    sinh_2_2m = (exp_2_2m_pos - exp_2_2m) / 2
    return exp_1, exp_1_2m, exp_1_2m_neg, exp_1_2m_gamma, sinh_1_2m, \
           exp_2, exp_2_2m, exp_2_2m_pos, exp_2_2m_gamma, sinh_2_2m


def r_2m(m, exp_1, exp_1_2m, exp_1_2m_gamma, exp_2, exp_2_2m_pos, exp_2_2m_gamma):
    return 2 * m * (1 - exp_1_2m) * exp_1 * a_2m(m, exp_1_2m_gamma) + \
           2 * m * binom(2 * m, m) * (1 - exp_2_2m_pos) + \
           2 * m * (1 - exp_2_2m_pos) * exp_2 * a_2m(m, exp_2_2m_gamma)


def q_2m(m, exp_1, exp_1_2m_neg, exp_1_2m_gamma, sinh_1_2m, 
            exp_2, exp_2_2m_pos, exp_2_2m_gamma, sinh_2_2m):
    return 2 * sinh_1_2m * exp_1 * exp_1_2m_gamma * a_2m_der(m, exp_1_2m_gamma) + \
           m * binom(2 * m, m) * (exp_2_2m_pos - exp_1_2m_neg) + \
           4 * m * sinh_2_2m * exp_2 * a_2m(m, exp_2_2m_gamma) - \
           2 * sinh_2_2m * exp_2 * exp_2_2m_gamma * a_2m_der(m, exp_2_2m_gamma)


def func_2m(alpha, temp, m, gamma):
    exp_1, exp_1_2m, exp_1_2m_neg, exp_1_2m_gamma, sinh_1_2m, \
    exp_2, exp_2_2m, exp_2_2m_pos, exp_2_2m_gamma, sinh_2_2m = init_all_exp_2m(temp, alpha, m, gamma)
    return 2 * m - 2 * gamma * m * (1 - exp_2_2m) + gamma / exp_2 * \
           (r_2m(m, exp_1, exp_1_2m, exp_1_2m_gamma, exp_2, exp_2_2m_pos, exp_2_2m_gamma) +
            q_2m(m, exp_1, exp_1_2m_neg, exp_1_2m_gamma, sinh_1_2m, exp_2, exp_2_2m_pos, exp_2_2m_gamma, sinh_2_2m))



def init_all_exp_2m_1(temp, alpha, m, gamma):
    exp_1 = np.exp(-(alpha - 1) * gamma / temp)
    exp_1_2m_1 = np.exp(2 * (alpha - 1) / temp / (2 * m + 1))
    exp_1_2m_1_neg = 1.0 / exp_1_2m_1
    exp_1_2m_1_gamma = np.exp(2 * (alpha - 1) * gamma / temp / (2 * m + 1))
    exp_1_2m_1_gamma_1 = exp_1_2m_1 / exp_1_2m_1_gamma
    sinh_1_2m_1 = (exp_1_2m_1 - exp_1_2m_1_neg) / 2
    sqrt_exp_1_2m_1_pos = np.exp((alpha - 1) / temp / (2 * m + 1))
    sqrt_exp_1_2m_1_neg = 1.0 / sqrt_exp_1_2m_1_pos
    exp_2 = np.exp(alpha * gamma / temp)
    exp_2_2m_1 = np.exp(-2 * alpha / temp / (2 * m + 1))
    exp_2_2m_1_pos = 1.0 / exp_2_2m_1
    exp_2_2m_1_gamma = np.exp(-2 * alpha * gamma / temp / (2 * m + 1))
    exp_2_2m_1_gamma_1 = exp_2_2m_1 / exp_2_2m_1_gamma
    sinh_2_2m_1 = (exp_2_2m_1_pos - exp_2_2m_1) / 2
    sqrt_exp_2_2m_1_pos = np.exp(alpha / temp / (2 * m + 1))
    sqrt_exp_2_2m_1_neg = 1.0 / sqrt_exp_2_2m_1_pos
    return exp_1, exp_1_2m_1, exp_1_2m_1_neg, exp_1_2m_1_gamma, exp_1_2m_1_gamma_1, sinh_1_2m_1, sqrt_exp_1_2m_1_pos, sqrt_exp_1_2m_1_neg, \
           exp_2, exp_2_2m_1, exp_2_2m_1_pos, exp_2_2m_1_gamma, exp_2_2m_1_gamma_1, sinh_2_2m_1, sqrt_exp_2_2m_1_pos, sqrt_exp_2_2m_1_neg



def r_2m_1(m, exp_1, exp_1_2m_1, exp_1_2m_1_gamma, exp_1_2m_1_gamma_1, sqrt_exp_1_2m_1_pos, 
              exp_2, exp_2_2m_1_pos, exp_2_2m_1_gamma, sqrt_exp_2_2m_1_pos):
    return (2 * m + 1) * (1 - exp_1_2m_1) * exp_1 * a_2m_1(m, exp_1_2m_1_gamma) + \
           (2 * m + 1) * binom(2 * m + 1, m) * exp_1_2m_1_gamma_1 * (sqrt_exp_1_2m_1_pos - sqrt_exp_2_2m_1_pos) + \
           (2 * m + 1) * (1 - exp_2_2m_1_pos) * exp_2 * a_2m_1(m, exp_2_2m_1_gamma)


def q_2m_1(m, exp_1, exp_1_2m_1_gamma, exp_1_2m_1_gamma_1, sinh_1_2m_1, sqrt_exp_1_2m_1_pos, sqrt_exp_1_2m_1_neg,
              exp_2, exp_2_2m_1_gamma, exp_2_2m_1_gamma_1, sinh_2_2m_1, sqrt_exp_2_2m_1_pos, sqrt_exp_2_2m_1_neg):
    return 2 * sinh_1_2m_1 * exp_1 * exp_1_2m_1_gamma * a_2m_1_der(m, exp_1_2m_1_gamma) + \
           m * binom(2 * m + 1, m) * exp_1_2m_1_gamma_1 * (sqrt_exp_2_2m_1_pos - sqrt_exp_1_2m_1_pos) + \
           (m + 1) * binom(2 * m + 1, m + 1) * exp_2_2m_1_gamma_1 * (sqrt_exp_2_2m_1_neg - sqrt_exp_1_2m_1_neg) + \
           2 * (2 * m + 1) * sinh_2_2m_1 * exp_2 * a_2m_1(m, exp_2_2m_1_gamma) - \
           2 * sinh_2_2m_1 * exp_2 * exp_2_2m_1_gamma * a_2m_1_der(m, exp_2_2m_1_gamma)


def func_2m_1(alpha, temp, m, gamma):
    exp_1, exp_1_2m_1, exp_1_2m_1_neg, exp_1_2m_1_gamma, exp_1_2m_1_gamma_1, sinh_1_2m_1, sqrt_exp_1_2m_1_pos, sqrt_exp_1_2m_1_neg, \
    exp_2, exp_2_2m_1, exp_2_2m_1_pos, exp_2_2m_1_gamma, exp_2_2m_1_gamma_1, sinh_2_2m_1, sqrt_exp_2_2m_1_pos, sqrt_exp_2_2m_1_neg = init_all_exp_2m_1(temp, alpha, m, gamma)
    return 2 * m + 1 - gamma * (2 * m + 1) * (1 - exp_2_2m_1) + gamma / exp_2 * \
           (r_2m_1(m, exp_1, exp_1_2m_1, exp_1_2m_1_gamma, exp_1_2m_1_gamma_1, sqrt_exp_1_2m_1_pos, 
                      exp_2, exp_2_2m_1_pos, exp_2_2m_1_gamma, sqrt_exp_2_2m_1_pos) +
            q_2m_1(m, exp_1, exp_1_2m_1_gamma, exp_1_2m_1_gamma_1, sinh_1_2m_1, sqrt_exp_1_2m_1_pos, sqrt_exp_1_2m_1_neg, 
                      exp_2, exp_2_2m_1_gamma, exp_2_2m_1_gamma_1, sinh_2_2m_1, sqrt_exp_2_2m_1_pos, sqrt_exp_2_2m_1_neg))


def find_root(temp, alpha0=0, alpha1=1, g=2, gamma=2):
    if g % 2 == 0:
        return brentq(lambda x: func_2m(x, temp, int(g / 2), gamma), alpha0, alpha1)
    else:
        return brentq(lambda x: func_2m_1(x, temp, int((g - 1) / 2), gamma), alpha0, alpha1)


def construct_phase_diagram(temp_list, g, gamma, fileout):
    w = open(fileout, "w")
    for temp in temp_list:
        try:
            alpha = find_root(temp, g=g, gamma=gamma)
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
    g = 4
    c = 4
    temp_f = 2
    dtemp = 0.01

    mean_gamma = c - 1
    graph_type = 'RGfc'
    fileout = "Stability_Ferro_Homophily_" + graph_type + "_c_" + str(c) + "_G_" + str(g) + ".txt"
    temp_list = get_temp_list(temp_f, dtemp)

    construct_phase_diagram(temp_list, g, mean_gamma, fileout)

    return 0


if __name__ == '__main__':
    main()
