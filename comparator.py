#! /usr/bin/python3

import getopt
import sys

import matplotlib.pyplot as plt
import numpy as np


def read_file(file: str) -> dict:
    """ Reading a file with n-meres and its z-score/counts [and optionaly with hydrophobicity]"""

    S = {}
    f = open(file).readlines()
    for i in f:
        i = i.strip().split("\t")
        S[i[0]] = [float(x) for x in i[1:]]
    return S


def compareOUT(s1: dict, s2: dict, hydroph: bool, out: str) -> None:
    """ Saving to a file compared data from two files"""

    with open(out, "w") as f:
        cnt = 0
        for perm, z in s1.items():
            try:
                if hydroph:
                    s = perm + "\t" + str(z) + "\t"
                else:
                    s = perm + "\t" + str(z) + "\t" + str(s2[perm]) + "\n"
            except KeyError:
                s = perm + "\t" + str(z) + "\t" + str(None) + "\n"
            f.write(s)
            cnt += 1
            if cnt % 50000 == 0:
                print(round(cnt / len(s1.items()) * 100))


def compare(s1: dict, s2: dict, hydroph: bool) -> dict:
    S = {}
    for perm, data in s1.items():
        try:
            if hydroph:
                # sequence of n-mer; z-score for seq in i file; z-score in j file; hydroph. of seq
                S[perm] = [data[0], s2[perm][0], data[1]]
            else:
                # sequence of n-mer : z-score for seq in i file; z-score in j file
                S[perm] = [data[0], s2[perm][0]]
        except KeyError:
            pass
    return S


def draw(comp: dict, hydroph: bool, name: str) -> None:
    fig = plt.figure(figsize=(15, 15))
    X = np.empty(shape=len(comp.items()))
    Y = np.empty(shape=len(comp.items()))
    if hydroph:
        v = np.empty(shape=len(comp.items()))
    iter = 0
    for perm, z in comp.items():
        X[iter] = z[0]
        Y[iter] = z[1]
        if hydroph:
            v[iter] = z[2]
        iter += 1
    print("Plotting.")
    plt.figure(figsize=(25, 20))
    if hydroph:
        plt.scatter(X, Y, alpha=0.75, s=12, c=v, cmap="jet")
    else:
        plt.scatter(X, Y, alpha=0.75, s=12)
    plt.title("Comparision of z-scores", fontsize=23)
    plt.xlabel("IDP ZScores", fontsize=20)
    plt.ylabel("nonIDP ZScores", fontsize=20)
    plt.ylim(-20, 100)
    plt.xlim(-20, 100)
    if hydroph:
        plt.colorbar()
    print("Saving.")
    plt.savefig(name)


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hi:j:o:p:b")
    except getopt.GetoptError:
        print(
            "ERROR! Pass correct arguemnts\n",
            argv[0],
            "-i <inputfile1> -j <inputfile2> [-o <outplutfile> -p <plot filename> -b]",
        )
        sys.exit()
    if "-h" in argv:
        print(
            "This program comapres two files with n-meres, their z-scores and hydrophobicity\n"
            "The result is a plot comparing two sets of proteins\n",
            argv[0],
            "-i <inputfile1> -j <inputfile2> [-o <outplutfile> -p <plot filename> -b]\n"
            "-h for help\n"
            "-b if files contain hydrophobicity data\n",
        )
        sys.exit()
    if "-i" not in argv or "-j" not in argv or "-h" in argv:
        print(
            "ERROR! Pass correct arguemnts\n",
            argv[0],
            "-i <inputfile1> -j <inputfile2> [-o <outplutfile>]",
        )
        sys.exit()

    s1, s2 = None, None
    hydroph = False
    plot_it = False  # flag indicating to make a plot

    for opt, arg in opts:
        if opt == "-i":
            s1 = read_file(arg)
        if opt == "-j":
            s2 = read_file(arg)
        if opt == "-o":
            outputfile = arg
        if opt == "-b":
            hydroph = True
        if opt == "-p":
            plot_it = arg

    if "-o" not in argv:
        C = compare(s1, s2, hydroph)
    else:
        compareOUT(s1, s2, hydroph, outputfile)
    if plot_it:
        draw(C, hydroph, plot_it)


if __name__ == "__main__":
    main(sys.argv)
