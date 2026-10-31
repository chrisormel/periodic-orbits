

import rebound
import numpy as np
import pofun


def posol_init (self, chain, delta, *qList):
    """
    self = posol object

    This initialize a posol object

    - chain     :resonance chain (e.g., 3:2 or 4:2:1)
    - delta     :period ratio parameter
    - qList     :planet masses
    """

    self.np = len(qList)        #number of planets
    self.qarr = np.array((1.0,)+qList)
    self.cumq = np.cumsum(self.qarr)
    self.xvec = [delta]+list(qList)     #input parameters
    self.chain = chain  #e.g., 8:4:2:1
    self.nums = np.float64(chain.split(':'))
    self.prat21 = (delta+1) *self.nums[0]/self.nums[1] #period ratio 2nd-to-1st planet

    #solution types // e=ecc, ev=ecc.vector, theta=true longitude
    self.type_sol = ['e', 'ev', 'ev', 'theta']
    for i in range(self.np-2):
        self.type_sol += ['d','ev','ev','theta']#d=Δ (our def)


def x_to_orbital (self, xs):
    """
    returns a particle property list //P,e,m,omega//
    which REBOUND can use to initialize particles

    - npl       :number of particles
    - xs        :solution variables
    - type_xs   :solution type
    - xvec      :grid parameters (Δ, {masses})

    We're building a list of rebound input dicts:
        partL = [{planet 1}, {planet 2}, ...]


    Note: there may be some overlap of convert_to_rebound in genfun
    """

    #the first (inner-most) planet is always:
    partL = []
    dpart1 = {'P':2*np.pi, 'omega':0.0, 'theta':0.0, 'm':self.xvec[1]}
    dpart1['e'] = xs[0] #eccetricity
    partL.append(dpart1)

    #other planets
    for i in range(1,self.np):
        dpart = {'m':self.xvec[i+1]}
        if i==1:
            prat = self.prat21 #period of 2nd planet
            dpart['P'] = dpart1['P'] *prat

        if i==1:
            k0 = 1; ke=4
        else:
            k0 = 4*(i-1)
            ke = 4*i

        for k in range(k0,ke):
            if self.type_sol[k] in ['e', 'theta', 'l', 'omega']:
                dpart[self.type_sol[k]] = xs[k]
            elif self.type_sol[k]=='ev' and self.type_sol[k+1]=='ev': #eccentr. vector
                e,w = pofun.evec_to_ew(*xs[k:k+2])
                dpart['e'] = e
                dpart['omega'] = w
            elif self.type_sol[k]=='ev':
                pass

            elif self.type_sol[k]=='d':#delta
                pq = self.nums[i-1]/self.nums[i]
                prat *= (1+xs[k]) *pq #period ith planet
                dpart['P'] = dpart1['P'] *prat
            else:
                print('[pocore]:This should not occurr')
                sys.exit()

        partL.append(dpart)

    return partL


def setup_sim (self, x, integrator='ias15', itgopts={}):
    """
    Setup simulation. Jacobean orbital elements are used...

    x -- solution variables
    """

    if not np.isfinite(x).all():
        return None
    partL = self.x_to_orbital (x)

    if any(not 0<=part['e']<1 or not np.isfinite(part['P']) or part['P']<=0
           for part in partL):
        return None

    if integrator in ['whfast', 'whfast512']:
        ## scaling with (1-e)**1.5 to resolve pericenter passages.
        pscale = min([part['P']*(1-part['e'])**1.5/(2*np.pi) for part in partL])

    #add particles inside-out
    sim = rebound.Simulation()
    sim.add(m=1) #central mass

    for part in partL:
        sim.add(**part)

    #important!
    sim.move_to_com()

    if integrator=='ias15':
        #the absolute lower limit...
        #[25.11.01]: increased to speed up some...
        sim.integrator.min_dt = 1e-2 #rebound v 5.2.0 // should be enough for PO!
    elif integrator in ['whfast', 'whfast512']:
        ps = sim.particles
        amax = max([p.a*(1+p.e) for p in ps[1:]])
        sim.integrator = integrator
        sim.integrator.corrector = 17
        if integrator=='whfast':
            sim.integrator.safe_mode = 0
#       sim.exact_finish_time = 0       #now as argument apparently
        sim.exit_min_distance = 0.001    #failsafe criterion...
        sim.exit_max_distance = 1.5*amax #escape criterion...
        sim.dt = 0.02 *2*np.pi *pscale #pericenter scaling
    else:
        print('[pocore.setup_sim]:no valid integrator')
        sys.exit()

    #[26.09.27] allow non-default options
    for opt,val in itgopts.items():
        setattr(sim, opt, val)

    return sim
