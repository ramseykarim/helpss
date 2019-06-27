import numpy as np
import mpy_utils as mpu
from Greybody import Greybody
from scipy.optimize import minimize


ITER = {'a': 0}

# Standard guesses for Tc, Nh, Nc
# Nh, Nc are log10
standard_x0 = [10, 22, 23]

def goodness_of_fit_f_3p(x, dusts, obs, err, instr, Th, dof):
    ITER['a'] += 1
    # x is [Tc, Nh, Nc] (Ns in log10)
    src = Greybody([Th, x[0]], [10**N for N in x[1:]], dusts)
    # Computes reduced chi^2 given a source model and observations/errors
    return sum((d.detect(src) - o)**2 / (e*e) for d, o, e in zip(instr, obs, err))/dof


def fit_source_3p(observations, errors, detectors, dusts, Th=15., dof=1.):
    result = minimize(goodness_of_fit_f_3p,
        x0=standard_x0,
        args=(dusts, observations, errors, detectors, Th, dof),
        bounds=((0, None), (18, 25), (18, 25)),
        options={'maxiter': 50}
    )
    # print("RESULT:", result.message)
    return result.x


def goodness_of_fit_f_2p(x, dust, obs, err, instr, dof):
    ITER['a'] += 1
    # x is [Tc, Nc] (N in log10)
    src = Greybody(x[0], 10**x[1], dust)
    # Computes reduced chi^2 given a source model and observations/errors
    return sum((d.detect(src) - o)**2 / (e*e) for d, o, e in zip(instr, obs, err))/dof

def fit_source_2p(observations, errors, detectors, dust, dof=2.):
    result = minimize(goodness_of_fit_f_2p,
        x0=[standard_x0[0], standard_x0[2]],
        args=(dust, observations, errors, detectors, dof),
        bounds=((0, None), (18, 25)),
        options={'maxiter': 50}
    )
    # print("RESULT:", result.message)
    return result.x


def goodness_of_fit_f_1p(x, dust, obs, err, instr, Th, dof):
    ITER['a'] += 1
    # x is [Tc, Nc] (N in log10)
    src = Greybody(Th, 10**x[0], dust)
    # Computes reduced chi^2 given a source model and observations/errors
    return sum((d.detect(src) - o)**2 / (e*e) for d, o, e in zip(instr, obs, err))/dof

def fit_source_1p(observations, errors, detectors, dust, Th=15., dof=3.):
    result = minimize(goodness_of_fit_f_1p,
        x0=[20,],
        args=(dust, observations, errors, detectors, Th, dof),
        bounds=((18, 25),),
        options={'maxiter': 50}
    )
    return result.x


def bootstrap_errors(observations, errors, detectors, dusts,
    niter=30, fit_f=fit_source_2p, **kwargs):
    observations = np.array(observations)
    errors = np.array(errors)
    results = mpu.deque()
    for i in range(niter):
        obs_perturbed = np.random.normal(loc=observations,
            scale=errors)
        current_result = fit_f(obs_perturbed, errors, detectors, dusts, **kwargs)
        current_result[1:] = [10**x for x in current_result[1:]]
        results.append(current_result)
    fitted_param_sets = list(results)
    param_uncertainties = [np.std(x) for x in zip(*results)]
    return fitted_param_sets, param_uncertainties
