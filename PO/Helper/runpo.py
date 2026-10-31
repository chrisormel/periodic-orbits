import rebound
import pandas as pd
import pickle
import pocore
import numpy as np

class POsol:
    def __init__ (self, chain, delta, *qList, model='vector'):
        #essential initialization...
        pocore.posol_init(self, chain, delta, *qList)

#attach function from pocore to POsol class
POsol.x_to_orbital = pocore.x_to_orbital
POsol.setup_sim = pocore.setup_sim


def load_solution (ftable, cid, key):
    df = pd.read_pickle(ftable)
    chain = df.attrs['longchain']

    if cid not in df.index:
        raise SystemExit(f'Grid point cid={cid} does not exist in {ftable}.')
    solutions = df.at[cid,'xsols']
    if solutions is None or not solutions:
        raise SystemExit(f'No PO solutions exist at cid={cid} in {ftable}.')
    if key not in solutions:
        raise SystemExit(
            f'Solution {key!r} does not exist at cid={cid}. '
            f'Available keys: {", ".join(solutions)}'
        )

    #the grid and solution vectors
    xv = df.at[cid,'xvec']      #PO grid point  [delta, q1, q2...]
    xs = solutions[key]       #PO solution    [see notes]
    tpo = df.at[cid,'tpos'][key]#PO period

    print('The selected cube grid point:')
    print(f'  (cid,key): ({cid},{key}), resonance chain: {chain}')
    print(f'  delta: {xv[0]:.10g}')
    for i, mass in enumerate(xv[1:], start=1):
        print(f'  mass body {i} (relative to central mass): {mass:.10g}')

    print(f'\nsolution vector, xs:\n {xs} \n')
    return xv, xs, tpo, chain


def run_solution (chain, xv, xs, tarr, tpo):
    """integrates a sample solution"""

    posim = POsol (chain, *xv)      #PO simulation object
    dpartL = posim.x_to_orbital(xs) #convert PO-solution to orbital elements list
    sim = posim.setup_sim (xs)      #rebound simulation object

    ## integrate the orbit


    ps = sim.particles

    #print some output (arbitrary)
    print('{:>8} {:>8} {:>8} {:>10} {:>10} {:>8} {:>8}'.format('time', 'time/tpo', 'M1', 'e1', 'e2', 'w2-w1', 'θ2-θ1' ))
    print(' '.join('-' * width for width in (8, 8, 8, 10, 10, 8, 8)))
    for k,time in enumerate(tarr):
        sim.integrate(time)

        print('{:8.3f} {:8.3f} {:8.2f} {:10.2e} {:10.2e} {:8.2f} {:8.2f}'.format(
            time, time/tpo, ps[1].M, ps[1].e, ps[2].e, 
            ps[2].pomega -ps[1].pomega, ps[2].theta -ps[1].theta))


if __name__=="__main__":

    ftable = '../Cubes/321low.pkl'

    #solutions are identified by their cube-grid identifier "cid"
    #and multiplicity number "key"
    cid = 25000
    key = 'A0'

    xv, xs, tpo, chain = load_solution (ftable, cid, key)

    npo = 1    #run for 1 PO-cycle
    npt = 21
    tend = npo*tpo
    tarr = np.linspace(0, tend, npt)

    run_solution (chain, xv, xs, tarr, tpo)
