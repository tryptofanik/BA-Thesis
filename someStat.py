#! /usr/bin/python3

import collections
import getopt
import itertools
import math
import operator as op
import sys
from functools import reduce

import numpy as np


def ncr(n: int, r: int) -> float:
    """Counts Newton."""

    r = min(r, n - r)
    numer = reduce(op.mul, range(n, n - r, -1), 1)
    denom = reduce(op.mul, range(1, r + 1), 1)
    return numer / denom


def permutation_groups(n: int) -> dict:
    """ It returns all possible permutation groups for given length."""

    x = itertools.combinations_with_replacement("ACDEFGHIJKLMNPQRSTWY", r=n)
    S = {}
    for i in x:
        S["".join(i)] = 0
    return S


def read_data(file: str) -> dict:
    """Reading file with n-meres and their counts."""

    f = open(file).readlines()
    S = {}
    for i in f:
        i = i.split("\t")
        S[i[0]] = int(i[1])
    return S


def count_group_len(data: dict, n: int) -> collections.Counter:
    """It counts total occurances of all seqs in permutation group."""

    G = permutation_groups(n)
    C = collections.Counter()
    for i in G:
        C[i] = 0
    for k, v in data.items():
        k = "".join(sorted(k))
        C[k] += v
    return C


def calculateM(s: str) -> int:
    """ For given string calculates all possible permutations (m)."""

    a = list(s)
    n = len(a)
    cnts = [1 for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if a[i] == a[j]:
                cnts[i] += 1
    m = 1
    j = 0
    while n > 1:
        m *= ncr(n, cnts[j])
        n -= cnts[j]
        j += 1
    return m


def calculate_z_score(n_obs: int, n: int, p: int) -> float:
    """Calculate one z-score."""

    return round((n_obs - n * p) * 1.0 / math.sqrt(n * p * (1 - p)), 2)


def z_scores(data: dict, n: int) -> dict:
    """Calculate z-scores for all given sequences."""

    C = count_group_len(data, n)
    Z = {}
    cnt = 0
    for perm, Nobs in data.items():
        cnt += 1
        if cnt % 100000 == 0:
            print(round(cnt / len(data.items()) * 100, 1), "%")
        m = calculateM(perm)
        p = 1 / m
        group_len = C["".join(sorted(perm))]
        if n == 0 or p == 0 or p == 1:
            continue
        z = calculate_z_score(Nobs, group_len, p)
        Z[perm] = z
    return Z


def main(argv):
    try:
        opts, _ = getopt.getopt(argv[1:], "hi:n:o:")
    except getopt.GetoptError:
        print(
            "ERROR! Pass correct arguments\n",
            argv[0],
            "-i <inputfile> -n <n_meres> [-o <outplutfile>]",
        )
        return
    if "-h" in argv:
        print(
            """This program is for anylysis of short n-meres-peptides in proteins.\n
                    It returns z-scores of n-meres-peptides based on the size of its permutation group\n
                    Input file should contain nmeres with its counts in dataset.\n
                    Output file is n-mer with its z-score"""
        )
        return

    if "-n" not in argv or "-i" not in argv:
        print(
            "ERROR! Pass correct arguments\n",
            argv[0],
            "-i <inputfile> -n <n_meres> [-o <outplutfile>]",
        )
        return

    outputfile = argv[0] + "_ZScores"

    for opt, arg in opts:
        if opt == "-i":
            inputfile = arg
        if opt == "-n":
            n = int(arg)
        if opt == "-o":
            outputfile = arg

    data = read_data(inputfile)
    Z = z_scores(data, n)
    o = open(outputfile, "w")
    for i in sorted(Z.items(), key=lambda kv: kv[1], reverse=True):
        o.write(i[0])
        o.write("\t")
        o.write(str(i[1]))
        o.write("\n")
    # G = countGroupLen(data, 5)
    # for k,v in G.items():
    #     print(k, '\t', v)


if __name__ == "__main__":
    main(sys.argv)
