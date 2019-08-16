#! /usr/bin/python3
import collections
import getopt
import operator
import sys

import numpy as np
from Bio import SeqIO


def readData(file: str) -> dict:
    f = SeqIO.parse(file, "fasta")
    S = {}
    for rec in f:
        S[rec.id] = str(rec.seq)
    return S


def how_many_aa(data: dict) -> np.array:
    aa = np.zeros(data.size)
    for i in range(len(data)):
        aa[i] = len(data[i])
    return aa


def divide_sequence(seq: str, n: int) -> np.chararray:
    if len(seq) - n + 1 <= 0:
        raise Exception("Too short sequence")
    mers = np.chararray(shape=(len(seq) - n + 1), itemsize=n)
    if len(seq) < n:
        return

    for i in range(len(seq) - n + 1):
        mers[i] = seq[i : i + n]  # to imporve from O(n^2) -> O(n)
    return mers


def slices(data: dict, n: int) -> collections.Counter:
    X = collections.Counter()
    for id, seq in data.items():
        try:
            mers = divide_sequence(seq, n)
            for mer in mers:
                X[mer] += 1
        except Exception:
            continue
    return X


def main(argv):
    try:
        opts, _ = getopt.getopt(argv[1:], "hi:n:o:")
    except getopt.GetoptError:
        print(
            "ERROR! Pass correct arguemnts\n",
            argv[0],
            "-i <inputfile> -n <n_meres> [-o <outplutfile>]",
        )
        return

    outplutfile = ""

    for opt, arg in opts:
        if opt == "-h":
            print(
                """This program divides protein sequences into short n-meres-peptides\nIt accepts one .fasta file which is obligatory.\nYou also have to specify n value"""
            )
            print(argv[0], "-i <inputfile> -n <n_meres> [-o <outplutfile>]")
            return
        if opt == "-i":
            inputfile = arg
        if opt == "-n":
            n = int(arg)
        if opt == "-o":
            outplutfile = arg

    if "-n" not in argv or "-i" not in argv:
        print(
            "ERROR! Pass correct arguemnts\n",
            argv[0],
            "-i <inputfile> -n <n_meres> [-o <outplutfile>]",
        )
        return

    data = readData(inputfile)
    Slices = slices(data, n)

    sortedSlices = sorted(Slices.items(), key=operator.itemgetter(1), reverse=True)

    if outplutfile != "":
        f = open(outplutfile, "w")
    else:
        outplutfile = inputfile.split(".")[0] + "_" + str(n) + "-meres"
        f = open(outplutfile, "w")

    for k, v in sortedSlices:
        f.write(k.decode("utf-8"))
        f.write("\t")
        f.write(str(v))
        f.write("\n")


if __name__ == "__main__":
    main(sys.argv)
