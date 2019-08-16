#! /usr/bin/python3

import getopt
import os
import re
import sys
import time

import numpy as np
import pandas as pd


def runningmean(elements: list, n: int) -> np.array:
    """Function implementing running mean with increasing and decresing slicing 
    window of length n at the beginning and the end of a list."""

    if n % 2 == 0:
        raise AttributeError("n must be odd")
    R = np.empty((len(elements)))
    sum = 0
    for i in range(int(n / 2) + 1):
        sum += elements[i]
    for i in range(len(elements) - int(n / 2) - 1):
        R[i] = sum / min(n, i + int(n / 2) + 1)
        sum += elements[i + int(n / 2) + 1]
        if i - int(n / 2) >= 0:
            sum -= elements[i - int(n / 2)]
    for i in range(len(elements) - int(n / 2) - 1, len(elements)):
        R[i] = sum / min(n, len(elements) - i + int(n / 2))
        sum -= elements[i - int(n / 2)]

    return R


def slicer(M: np.array, anchordivide: bool, cutoff: int = 15) -> list:
    """Slice the sequence of proteins into idp-like and nonidp fragments.
    Parameters:
        M - a matrix from iupred; M=[aa,iu,anch] where aa is aminoacid, 
            iu - idp score; anch - achore score for predicting binding site of idp
        anchordivide - bool value which specify if sequeces should be divided using anchor score.
        cutoff - minimal length of the sequence
    """

    idp = ""
    non = ""
    currentidp = ""
    currentnon = ""

    toidp = False
    tonon = False

    if not anchordivide:
        runmean = runningmean(M[1], 15)
    else:
        runmean = runningmean(M[2], 15)

    for i in range(len(runmean)):
        if runmean[i] > 0.5:
            currentidp += M[0][i]
            if tonon and len(currentnon) > cutoff:
                non += currentnon + "\n"
                currentnon = ""
            toidp = True
            tonon = False
        elif runmean[i] < 0.5:
            currentnon += M[0][i]
            if toidp and len(currentidp) > cutoff:
                idp += currentidp + "\n"
                currentidp = ""
            tonon = True
            toidp = False
    if len(currentnon) > 10:
        non += currentnon + "\n"
    if len(currentidp) > 10:
        idp += currentidp + "\n"
    return [idp, non]


def divider(data: dict, out: str, anchordivide: bool, cutoff: int) -> None:
    # this divide sequences to idp and nonidp given data from iupred
    like = open("IDP_" + out + ".fasta", "a+")
    nonlike = open("nonIDP_" + out + ".fasta", "a+")

    for protId, M in data.items():
        if protId == "#":
            continue
        IDP, nonIDP = slicer(M, anchordivide, cutoff)
        like.write(protId + "\n")
        nonlike.write(protId + "\n")
        like.write(IDP + "\n")
        nonlike.write(nonIDP + "\n")

    like.close()
    nonlike.close()


def readData(file: str) -> dict:
    """Returns the dictionary where key is the name of a protein and value a list of data from iupred"""

    with open(file) as f:
        f = open(file)
        text = f.read()
    results = re.split("###\n", text)
    R = {}
    cnt = 0
    del results[-1]
    for rec in results:
        cnt += 1
        # teraz zajmuje sie jednym białkiem
        rec = rec.split("\n")
        id = rec[0]
        del rec[0]
        del rec[-1]
        aa = np.empty([len(rec)], dtype=str)
        iu = np.empty([len(rec)], dtype=float)
        anch = np.empty([len(rec)], dtype=float)
        for i in range(len(rec)):
            aaline = rec[i].split("\t")
            aa[i] = aaline[1]
            iu[i] = aaline[2]
            anch[i] = aaline[3]
        M = [aa, iu, anch]
        R[id] = M
    return R  #


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hi:o:ac:")
    except getopt.GetoptError:
        print(
            "ERROR! Pass correct arguemnts\n",
            argv[0],
            "-i <inputfile> -o <outplutfile> -c <cutoff> -a <divide by anchor score>",
        )
        sys.exit()

    outputfile = "seqs"
    data = None
    anchorDivide = False
    cutoff = 10

    for opt, arg in opts:
        if opt == "-h":
            print(
                "Script takes the data from iupred2a and divide it to two files. In one there are fragments of proteins which were classified as IDP-like,"
                "and in the other there are fragments classified as nonIDP. One protein sequence can be splited many time, in particular into very small fragments."
                "So the script takes advantage of running mean to make the predictive function smooth."
            )
            print(
                argv[0],
                "-i <inputfile> -o <outplutfile> [-c <cutoff> -a <divide by anchor score>]",
            )
            sys.exit()
        if opt == "-i":
            inputfile = os.path.join(os.getcwd(), arg)
            data = readData(inputfile)
        if opt == "-o":
            outputfile = arg
        if opt == "-a":
            anchorDivide = True
        if opt == "c":
            cutoff = arg

    divider(data, outputfile, anchorDivide, cutoff)


if __name__ == "__main__":
    start = time.time()
    main(sys.argv)
    end = time.time()
    print(end - start)
