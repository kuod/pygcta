import numpy as np
import pdb
import sys

from libs.kernel import *
from libs.phenotypes import *
from libs.genotypes import *
from libs.linearmodel import *

fn_bed = "data/test.bed"
fn_phen = "data/test.phen"
if __name__ == "__main__":

    ### get genotypes, phenotypes and one kernel
    GT = genotypes()
    GT.read_plink(fn_bed)
    K = GT.getK()

    Y = phenotypes(fn_phen)
    Y.read_phen(fn_phen = fn_phen)

    ### init linear model
    LMM = pygcta(Y = Y, K = [GT.getK()])
    LMM.optimize()

