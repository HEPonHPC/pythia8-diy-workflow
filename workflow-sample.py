#!/usr/bin/env python3

"""
    %prog steering.dat parameter.ranges  [options]
"""
import apprentice as app
import numpy as np

class Sobol(object):

    def __init__(self, dim, initialSeed=0):
        self._seed=initialSeed
        self._dim=dim

    def shoot(self):
        try:
            import sobol
        except ImportError:
            raise Exception("sobol not available, try pip install sobol")
        p, newseed = sobol.i4_sobol(self._dim, self._seed)
        self._seed = newseed
        return p

    def __call__(self):
        return self.shoot()

def writeFiles(proccard, pnames, points, outdir, fnamep="params.dat", fnameg="generator.cmd"):
    from os.path import join, exists
    for num, p in enumerate(points):
        npad = "{}".format(num).zfill(1+int(np.ceil(np.log10(len(points)))))
        outd = join(outdir, npad)
        if not exists(outd):
            import os
            os.makedirs(outd)

        outfparams = join(outd, fnamep)
        with open(outfparams, "w") as pf:
            for k, v in zip(pnames, p):
                pf.write("{name} {val:e}\n".format(name=k, val=v))

        outfgenerator = join(outd, fnameg)
        with open(outfgenerator, "w") as pg:
            for l in proccard:
                pg.write(l+"\n")
            for k, v in zip(pnames, p):
                pg.write("{name} = {val:e}\n".format(name=k, val=v))


def sample(ranges, npoints, algorithm="uni", seed=1234):

    names = [r[0] for r in ranges]
    DIM = len(names)
    pmin  = np.array([r[1] for r in ranges])
    pmax  = np.array([r[2] for r in ranges])
    LEN = pmax - pmin

    if algorithm=="uni":
        RSP = np.random.uniform(low=pmin, high=pmax,size=(npoints, DIM))

    elif algorithm=="lhs":
        try:
            import pyDOE2
        except ImportError:
            raise Exception("pyDOE2 not available, try pip install pyDOE2")
        RSP = pmin + LEN*pyDOE2.lhs(DIM, samples=npoints)

    elif algorithm=="sobol":
        s = Sobol(DIM, seed)
        S = np.zeros((npoints, DIM))
        for i in range(npoints): S[i] = s()
        RSP = pmin + LEN*S
    else:
        raise Exception("Sampling algorithm {} is not a valid choice.".format(algorithm))

    return names, RSP



def readProcessCard(fname):
    with open(fname) as f:
        L = [l.strip() for l in f]
    return L


def readRanges(fname):
    with open(fname) as f:
        L = [l.strip() for l in f]
    L = [l for l in L if not l.startswith("#") and not len(l)==0]
    if len(L) == 0:
        raise Exception("There are no free parameters")

    RR = []

    for l in L:
        name, pmin, pmax = l.split()
        pmin = float(pmin)
        pmax = float(pmax)
        RR.append((name, pmin, pmax))
    return RR

if __name__=="__main__":

    import optparse, os, sys
    op = optparse.OptionParser(usage=__doc__)
    op.add_option("-v", "--debug", dest="DEBUG"  , action="store_true", default=False, help="Turn on some debug messages")
    op.add_option("--force"      , dest="FORCE"  , action="store_true", default=False, help="Allow overwriting of output directory.")
    op.add_option("-o"           , dest="OUTDIR" , default="scan", help="Output directory (default: %default)")
    op.add_option("--seed"       , dest="SEED"   , default=1234, type=int, help="The base random seed (default: %default)")
    op.add_option("--order"      , dest="ORDER"  , default="3,0", type=str, help="The polynomial orders to use --- this, together with the dimension of the parameter space determines the number of points sampled. (default: %default)")
    op.add_option("--sampler"    , dest="SAMPLER", default="uni", type=str, help="The sampling algorithm uni|sobol|lhs. (default: %default)")
    op.add_option("--fidelity"   , dest="FIDEL"  , default=2, type=float, help="The degree of oversampling. The minimal number of points required is determined by the dimension of the parameter space and the polynomial orders. The degree of oversampling acts as a multiplier on that number (default: %default)")
    op.add_option("--pfname"     , dest="PFNAME" , default="params.dat", type=str, help="The name for parameter files. (default: %default)")
    op.add_option("--gfname"     , dest="GFNAME" , default="generator.cmd", type=str, help="The name for generator input cards. (default: %default)")
    opts, args = op.parse_args()


    if (opts.FIDEL) <1:
        print("Error, oversampling must be greater or equals to 1, you supplied {}, exiting...".format(opts.FIDEL))
        sys.exit(1)

    if len(args)!=2:
        print("Error, invalid number of arguments. Must be 2, you povided {}, exiting...".format(len(args)))
        sys.exit(1)

    if opts.SAMPLER not in ["uni", "lhs", "sobol"]:
        print("Invalid sampling algorithm specified: {}, exiting...".format(opts.SAMPLER))
        exit(1)

    try:
        pc = readProcessCard(args[0])
        rr = readRanges(args[1])
    except Exception as e:
        print("An error occured when reading the input files: {}".format(e))
        exit(1)

    if os.path.exists(opts.OUTDIR) and not opts.FORCE:
        print("Output directory {} exists. Use --force to allow overwriting or chose a different location with -o, exiting...".format(opts.OUTDIR))
        sys.exit(1)


    DIM = len(rr)
    M,N = opts.ORDER.split(",")
    M = int(M)
    N = int(N)
    ncoeff  = app.tools.numCoeffsRapp(DIM, (M,N))
    npoints = ncoeff*opts.FIDEL
    print("Sampling {} points with {} algorithm.".format(npoints, opts.SAMPLER))

    pnames, sampledpoints = sample(rr, npoints, opts.SAMPLER, opts.SEED)


    writeFiles(pc, pnames, sampledpoints, opts.OUTDIR, opts.PFNAME, opts.GFNAME)
