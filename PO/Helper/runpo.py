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


cubedir = '../Cubes/'
ftable = '21low.pkl'

df = pd.read_pickle(cubedir+ftable)
chain = df.attrs['longchain']

#solutions are identified by their cube-grid identifier "cid"
#and multiplicity number "key"
cid = 1000
key = 'A0'

#the grid and solution vectors
xv = df.at[cid,'xvec']      #PO grid point  [delta, q1, q2...]
xs = df.at[cid,'xsols'][key]#PO solution    [see notes]
tpo = df.at[cid,'tpos'][key]#PO period


posim = POsol (chain, *xv)   #PO simulation object
dpartL = posim.x_to_orbital(xs)
sim = posim.setup_sim (xs) #rebound simulation object

## integrate the orbit

npo = 1         #number of PO-cycles to run
npt = npo*100
tend = npo*tpo
tarr = np.linspace(sim.t, tend, npt)

ps = sim.particles
for k,time in enumerate(tarr):
    sim.integrate(time)

    print('{:10.3f} {:10.2f} {:10.2f}'.format(time/tpo, ps[1].M, ps[1].theta))


