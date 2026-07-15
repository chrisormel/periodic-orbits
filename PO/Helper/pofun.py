

import rebound
import numpy as np


def evec_to_ew (he, ke):

    e2 = np.sqrt(he**2 +ke**2)
    w2 = np.arctan2(ke,he)

    return e2, w2
