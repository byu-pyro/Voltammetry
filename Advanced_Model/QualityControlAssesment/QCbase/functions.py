import pickle
from matplotlib import pyplot as plt
import numpy as np
import scipy.optimize as opt
import time

# constants
R = 8.314  # J/mol-K
F = 96485  # C/mol e-


# a function to compute the density of the salts at a given temperature
def get_density_salts(T, density_parameters):
    """
    :param T: Temperature in Kelvin
    :param density_parameters: 2D array of the density parameters (salts, parameters)
    :return: returns the density of each salt from a linear equation
    """
    return density_parameters[:, 0] + density_parameters[:, 1] * T


# a function to compute the diffusivity of ions at a given temperature
def get_D_ions(T, diffusivity_parameters):
    """
    :param T: temperature in Kelvin
    :param diffusivity_parameters: 2d array of the density parameters (ions, parameters)
    :return:  returns the diffusivity of each ion in mol/cm^2
    """
    return diffusivity_parameters[:, 0] * np.exp(diffusivity_parameters[:, 1] / R / T)


# a function to get the desired grid space
def get_gridspace(L, n_x, alpha):
    """
    :param L: Distance in cm to the point that concentration is considered bulk
    :param n_x: the number of grid points wanted
    :param alpha: Increase each grid space by a factor of alpha
    :return: 1D array of the gridspace
    """
    dx = L / np.sum(alpha ** np.array(range(n_x)))
    x = np.array([0.0])
    for j in range(n_x):
        x = np.append(x, alpha ** j * dx + x[-1])
    return x


# a function to get the average molecular weight
def get_MW_ave(y_salts, MW_salts):
    """
    :param y_salts: mass fraction of the salts
    :param MW_salts: molecular weights of the salts (g/mol)
    :return: the average molecular weight (g/mol)
    """
    return 1 / np.sum(y_salts / MW_salts)


def get_x_salts(y_salts, MW_salts):
    """
    :param y_salts: mass fraction of the salts
    :param MW_salts: molecular weights of the salts (g/mol)
    :return: the mole fractions of the salts
    """
    return y_salts * get_MW_ave(y_salts, MW_salts) / MW_salts


def get_C_ions(T, y_salts, MW_salts, density_parameters, stoic_ions):
    """
    :param T: temperature in Kelvin
    :param y_salts: mass fraction of the salts
    :param MW_salts: molecular weight of the salts (g/mol)
    :param density_parameters: density parameters to get density as a function of temperature
    :param stoic_ions: stoichiometry of ions in the salts.
    :return: concentration of the ions in mol/cm^3
    """

    # get density of salts
    density_salts = get_density_salts(T, density_parameters)

    # get salt mole fraction
    x_salts = get_x_salts(y_salts, MW_salts)

    # molar average of each salt density
    density_ave = np.sum(density_salts * x_salts)  # g/cm^3

    # average molecular weight of the salts
    MW_ave = get_MW_ave(y_salts, MW_salts)  # g/mol^3

    # concentration of the ions
    C_ions = np.dot(x_salts, stoic_ions) / MW_ave * density_ave  # mol/cm^3

    return C_ions


# a function to get the mobility of ions
def get_u_ions(T, z_ions, diffusivity_parameters):
    """
    :param T: temperature in Kelvin
    :param z_ions: the charge of the ions
    :param diffusivity_parameters: the diffusivity parameters of the ions
    :return: the mobilities of the ions
    """
    D_ions = get_D_ions(T, diffusivity_parameters)
    # mobilities
    return np.abs(z_ions) * F * D_ions / R / T


# a function to calculate ohmic resistance based on bulk conditions.
def get_R_ohmic(L_RE, C_ions, T, z_ions, diffusivity_parameters):
    """
    :param L_RE: The distance to the reference electrode
    :param C_ions: the concentration of ions
    :param T: the temperature in Kelvin
    :param z_ions: the charge of the ions
    :param diffusivity_parameters: the diffusivity parameters of the ions
    :return: The ohmic resistance of the system
    """
    # get the mobility of ions
    u_ions = get_u_ions(T, z_ions, diffusivity_parameters)

    # bulk conductivity of the salt
    conductivity = F * np.sum(np.abs(z_ions) * u_ions * C_ions)

    return L_RE / conductivity


# a function to calculate the solid activity based on the amount of deposit
def solid_activity(N, ideal_solid_activity, No):
    """
    :param N: the amount of deposit (mol/cm^2)
    :param ideal_solid_activity: boolean to include solid activity or not
    :param No: amount of deposit to have full surface coverage (mol/cm^2)
    :return: the solid activity (unitless)
    """

    if ideal_solid_activity:
        return 1
    else:
        coverage = N / No
        return 1 - np.exp(-3.91e-3 * coverage - 6.9 * coverage ** 10.5)


# a function to determine the ion activity from concentration
def ion_activity(C, n, ideal_ion_activity, dEo, T):
    """
    :param C: The concentration of the analyte in mol/cm^2
    :param n: the number of electrons transferred
    :param ideal_ion_activity: boolean of the ion activity inclusion
    :param dEo: the difference in pure and dilute standard potentials
    :param T: the temperature in K
    :return: the ion activity/unitless
    """
    if ideal_ion_activity:
        return C / 1
    else:
        f = 1 - np.exp(-3.26295483e+01 * C - 2.29456177e+06 * C ** 2.94736811e+00)
        gamma = np.exp(dEo * n * F / R / T * f)
        return gamma * C / 1


# a function to determine the equilibrium potential from the Nernst equation
def Nernst(N, C, Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No):
    """
    :param N: the amount of deposit (mol/cm^2)
    :param C: the concentration of the analyte at the surface (M)
    :param Eo: the standard potential (V)
    :param n: the number of electrons transferred
    :param T: the temperature (K)
    :param ideal_solid_activity: boolean of solid activity inclusion
    :param ideal_ion_activity: boolean of ion activity inclusion
    :param dEo: the difference in pure and dilute standard potentials (V)
    :param No: the amount of deposit to have full surface coverage (mol/cm^2)
    :return: the equilibrium potential (V)
    """
    a_solid = solid_activity(N, ideal_solid_activity, No)
    a_ion = ion_activity(C, n, ideal_ion_activity, dEo, T)
    return Eo - R * T / n / F * np.log(a_solid / a_ion)


# a function that gets the potential based on ideal solid activity
def get_E_star(C, Eo, n, T, ideal_ion_activity, dEo):
    """
    :param C: the concentration of the analyte at the surface (M)
    :param Eo: the standard potential (V)
    :param n: the number of electrons transferred
    :param T: the temperature (K)
    :param ideal_ion_activity: boolean of ion activity inclusion
    :param dEo: the difference in pure and dilute standard potentials (V)
    :return: the equilibrium potential with ideal solid activity (V)
    """
    return Nernst(1, C, Eo, n, T, True, ideal_ion_activity, dEo, 1)


# a function to get the amount of deposit based on an equilibrium potential and concentration
def get_N(N_guess, E, C, Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No):
    """
    :param N_guess: A guess at the amount of deposit (mol/cm^2)
    :param E: the equilibrium potential (V)
    :param C: the concentration of the analyte at the surface (M)
    :param Eo: the standard potential (V)
    :param n: the number of electrons transferred
    :param T: the temperature (K)
    :param ideal_solid_activity: boolean of solid activity inclusion
    :param ideal_ion_activity: boolean of ion activity inclusion
    :param dEo: the difference in pure and dilute standard potentials (V)
    :param No: the amount of deposit to have full surface coverage (mol/cm^2)
    :return: the amount of deposit at the surface (mol/cm^2)
    """

    if ideal_solid_activity:
        return 0.0

    def solver(N):
        return E - Nernst(N, C, Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No)

    solution = opt.fsolve(solver, np.array([N_guess]), full_output=True)
    return solution[0][0]


def get_C(C_guess, E, N, Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No):
    """
    :param C_guess: A guess of the surface concentration of the analyte (M)
    :param E: the equilibrium potential (V)
    :param N: the amount of surface deposit (mol/cm^2)
    :param Eo: the standard potential (V)
    :param n: the number of electrons transferred
    :param T: the temperature (K)
    :param ideal_solid_activity: boolean of solid activity inclusion
    :param ideal_ion_activity: boolean of ion activity inclusion
    :param dEo: the difference in pure and dilute standard potentials (V)
    :param No: the amount of deposit to have full surface coverage (mol/cm^2)
    :return: The concentration of the analyte at the surface (M)
    """

    def solver(C):
        return E - Nernst(N, C, Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No)

    solution = opt.fsolve(solver, np.array([C_guess]), full_output=True)
    return solution[0]


def get_diffusion(C, D, x):
    return -D[:, None] * (C[:, 1:] - C[:, :-1]) / (x[1:] - x[:-1])


def get_potential_gradient(C, i, z, D, T, x):
    term1 = F * np.sum(z[:, None] * D[:, None] * (C[:, 1:] - C[:, :-1]) / (x[1:] - x[:-1]), axis=0)
    term2 = F * np.sum(z[:, None] ** 2 * F * D[:, None] / R / T * (C[:, 1:] + C[:, :-1]) / 2, axis=0)
    return -(i + term1) / term2


def get_migration(C, i, z, D, T, x):
    return -z[:, None] * F * D[:, None] / R / T * (C[:, 1:] + C[:, :-1]) / 2 * get_potential_gradient(C, i, z, D, T, x)


def get_flux(C, i, no_migration, z, D, T, x):
    if no_migration:
        return get_diffusion(C, D, x)
    else:
        return get_diffusion(C, D, x) + get_migration(C, i, z, D, T, x)


def get_ohmic_losses(C, i, no_ohmic_losses, z, D, T, x):
    if no_ohmic_losses:
        E_ohmic = 0
    else:
        E_ohmic = -np.sum(get_potential_gradient(C, i, z, D, T, x) * (x[1:] - x[:-1]))
    return E_ohmic


def LSV(t, v, E_start):
    E = E_start - v * t
    return E


def open_sim(sim):
    t = sim['data']['t'][-1]
    E = sim['data']['E'][-1]
    i = sim['data']['i'][-1]
    C = sim['data']['C'][-1, :, :]
    N = sim['data']['N'][-1]
    capture_rate = sim['data']['capture_rate']
    compute_time = sim['data']['compute_time'][-1]

    x = sim['parameters']['x']
    T = sim['parameters']['T']
    v = sim['parameters']['v']
    n = sim['parameters']['n']
    Eo = sim['parameters']['Eo']
    D = get_D_ions(T, sim['parameters']['diffusivity_parameters_ions'])
    z = sim['parameters']['z_ions']
    stoic_rxn = sim['parameters']['stoic_rxn']
    j_analyte = sim['parameters']['j_analyte']
    E_start = sim['parameters']['E_start']
    E_star = sim['parameters']['E_star']

    ideal_solid_activity = sim['assumptions']['ideal_solid_activity']
    ideal_ion_activity = sim['assumptions']['ideal_ion_activity']
    no_migration = sim['assumptions']['no_migration']
    no_ohmic_losses = sim['assumptions']['no_ohmic_losses']
    iR_compensation = sim['assumptions']['iR_compensation']
    No = sim['assumptions']['No']
    dEo = sim['assumptions']['dEo']
    L_RE = sim['assumptions']['L_RE']
    R_ohmic = sim['assumptions']['R_ohmic']

    return t, E, i, C, N, compute_time, x, T, v, n, Eo, D, z, stoic_rxn, capture_rate, j_analyte, E_start, E_star, \
        ideal_solid_activity, ideal_ion_activity, no_migration, no_ohmic_losses, iR_compensation, No, dEo, L_RE, R_ohmic


def export(sim, debug):
    exporting = True
    while exporting:
        try:
            export_time = time.time()
            sim['export_time'] = export_time
            file = open(f'{sim["path"]}{sim["file_name"]}', 'wb')
            pickle.dump(sim, file)
            file.close()
            exporting = False
        except:
            if debug: print("error while exporting, trying again")
    return export_time


def get_current(N_i, C, dt, i_i, E_app, J, stoic_rxn, n, T, D, dEo, No, z, x, ideal_solid_activity,
                ideal_ion_activity, no_ohmic_losses, Eo, debug):
    # solver for N, i, and surface concentrations
    def solver(param):
        N, i, *C_surf = param

        # material balance for the deposit
        eq = np.array([N_i - i / n / F * dt - N])

        # equilibrium concentration at surface
        E_ohmic = get_ohmic_losses(C, i, no_ohmic_losses, z, D, T, x)
        E_WE = E_app - E_ohmic
        eq = np.append(eq, Nernst(N, C_surf[2], Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No) - E_WE)

        # material balance on ions
        eq = np.append(eq, C[:, 0] + (-stoic_rxn * i / n / F - J[:, 0]) / (x[1] / 2) * dt - C_surf)
        return eq

    # get guess values
    # ------------------------------------------------------------------------------------------------------------------
    E_ohmic = get_ohmic_losses(C, i_i, no_ohmic_losses, z, D, T, x)
    E_WE = E_app - E_ohmic
    E_eq = Nernst(N_i, C[2, 0], Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No)
    dx = x[1] / 2

    if N_i < No and not ideal_solid_activity:  # reaction is maily dependent on solid activity
        if debug: print("Guess based on ion")
        N = get_N(N_i, E_WE, C[2, 0], Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No)
        i = -(N - N_i) * n * F / dt
        C_surf = C[:, 0] + (-stoic_rxn * i / n / F - J[:, 0]) / (x[1] / 2) * dt
    else:  # reaction is mainly dependent on ion activity
        if debug: print("Guess based on solid")
        C_U = get_C(C[2, 0], E_WE, N_i, Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No)
        i = ((C_U - C[2, 0]) / dt * dx + J[2, 0]) * F * n / -stoic_rxn[2]
        C_surf = C[:, 0] + (-stoic_rxn * i / n / F - J[:, 0]) / dx * dt
        C_surf[2] = C_U
        N = N_i - i / n / F * dt
        if N < 0.0:  # if the guess value is oxidizing more than the amount of deposit
            if debug: print("not enough deposit")
            N = N_i
            i = 0.0
            C_surf = C[:, 0]

    # if i != 0.0: debug = True

    # start printing guess values and solutions if guess values don't make sense.
    # if N < 0.0 or np.any(np.array(C_surf) <= 0.0):
    #     debug = True

    if debug: print(f'Guess Values - \n\tN:{N}, \n\ti:{i}, \n\tC_surf:{C_surf}, \n\tdt: {dt}, \n\tE_app: {E_app}'
                    f'\n\tE_WE: {E_WE}, \n\tN_i: {N_i}, \n\tC_i: {C[:, 0]}, \n\ti_i: {i_i}, \n\tJ: {J[:, 0]}')
    # if debug: input("wait")
    # ------------------------------------------------------------------------------------------------------------------

    # No need to run solver if oxidizing and no deposit
    if N == 0.0 and E_WE > E_eq:
        return N, i, *C_surf

    # run the solver
    solution = opt.fsolve(solver, np.append(np.array([N, i]), C_surf), full_output=True)
    N, i, *C_surf = solution[0]

    if debug: print(f'Solution - N:{N}, i:{i}, C_surf:{C_surf}')

    return N, i, *C_surf


def run(sim, debug):

    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel(r'Potential/V vs. $E^*$')
    ax.set_ylabel(r'Current Density/mA $cm^2$')
    ax.set_title(sim['file_name'])
    line, = ax.plot([], [])
    fig.show()
    plt.ion()

    # start a timer
    start_time = time.time()

    # mark the sim as running
    sim['status'] = 'running'

    # export the sim with updated status
    last_export_time = export(sim, debug)

    # save all the parameters to variables
    (t, E, i, C, N, compute_time, x, T, v, n, Eo, D, z, stoic_rxn, capture_rate, j_analyte, E_start, E_star,
     ideal_solid_activity, ideal_ion_activity, no_migration, no_ohmic_losses, iR_compensation, No, dEo, L_RE, R_ohmic) \
        = open_sim(sim)

    def attempt_timestep(dt, t, N, C, i, n_timesteps):

        # This code ensure the time steps are 10 times smaller than the max timestep for stability
        # -------------------------------------------------------------------------------------------------------------

        # check the max timestep for stability
        max_dt = x[1] ** 2 / D[j_analyte] / 2
        max_dt = max_dt / 10

        # if the time step is too big
        if dt > max_dt:

            # calculate how many times the time step should be divided
            n_div = int(dt / max_dt + 1)

            # loop through that many times to do smaller timesteps
            for foo in range(n_div):
                t, N, C, i, n_timesteps = attempt_timestep(dt / n_div, t, N, C, i, n_timesteps)
            return t, N, C, i, n_timesteps
        # -------------------------------------------------------------------------------------------------------------

        # solve the surface conditions
        # -------------------------------------------------------------------------------------------------------------

        # get the flux
        J = get_flux(C, i, no_migration, z, D, T, x)

        # get the potential from the waveform
        E_waveform = LSV(t + dt, v, E_start)

        # correct the potential for ohmic losses
        E_app = E_waveform + iR_compensation * i * R_ohmic

        N_temp, i_temp, *C_temp = get_current(N, C, dt, i, E_app, J, stoic_rxn, n, T, D, dEo, No,
                                              z, x, ideal_solid_activity, ideal_ion_activity, no_ohmic_losses, Eo,
                                              debug)
        # -------------------------------------------------------------------------------------------------------------

        # check if a realistic answer
        # -------------------------------------------------------------------------------------------------------------
        success = True
        if N_temp < 0.0 or np.any(np.array(C_temp) <= 0.0):
            if debug: print("Negative Result")
            success = False

        # Add any other success criteria here

        if not success:
            n_div = 2
            for foo in range(n_div):
                t, N, C, i, n_timesteps = attempt_timestep(dt / n_div, t, N, C, i, n_timesteps)
            if debug: print(f'after cutting - N:{N}, i:{i}, C_surf:{C[:, 0]}')
        else:
            N, i, C[:, 0] = N_temp, i_temp, C_temp
            C[:, 1:-1] += (J[:, :-1] - J[:, 1:]) / ((x[2:] + x[1:-1]) / 2 - (x[1:-1] + x[:-2]) / 2) * dt
            t += dt
            n_timesteps += 1
        # -------------------------------------------------------------------------------------------------------------
        return t, N, C, i, n_timesteps

    # get the maximum dt based on the capture rate
    dt = capture_rate / np.abs(v)

    # iteration counter
    iterations = 0

    # a stop flag
    stop = False
    while not stop:
        iterations += 1  # increment the counter

        n_timesteps = 0  # a counter for the number of recursive steps taken
        t, N, C, i, n_timesteps = attempt_timestep(dt, t, N, C, i, n_timesteps)
        E = LSV(t, v, E_start)


        # update sim
        sim['data']['t'] = np.append(sim['data']['t'], t)
        sim['data']['N'] = np.append(sim['data']['N'], N)
        sim['data']['C'] = np.append(sim['data']['C'], C[None, :,:], axis=0)
        sim['data']['i'] = np.append(sim['data']['i'], i)
        sim['data']['E'] = np.append(sim['data']['E'], E)

        if debug: print(f'E: {E}, E_star: {E_star}')
        # exit condition(s)
        if E < E_star:  # only start looking if we passed E_star
            i_peak = min(sim['data']['i'][np.where(sim['data']['E'] < E_star)])  # Find the peak
            if i > i_peak / 2:  # stop when current crops to half the peak
                stop = True

        if time.time() - last_export_time > 10:  # if it has been more than 10 seconds since last export
            last_export_time = export(sim, debug)
            print(f'Data Points: {iterations}, Computational Rate: '
                  f'{round(t*v*1000*60/(time.time()-start_time),2)} mV/min')
            iterations = 0
            #line.set_data(sim['data']['E'] - E_star, sim['data']['i'])
            #ax.relim()
            #ax.autoscale_view(True, True, True)
            #fig.canvas.draw()
            #fig.canvas.flush_events()
            #fig.tight_layout()

    last_export_time = export(sim, debug)
    
    # Plotting at the end of the full simulation
    #fig, ax = plt.subplots(1, 1)
    #ax.set_xlabel(r'Potential/V vs. $E^*$')
    #ax.set_ylabel(r'Current Density/mA $cm^2$')
    #ax.set_title(sim['file_name'])
    #line, = ax.plot([], [])
    #line.set_data(sim['data']['E'] - E_star, sim['data']['i'])
    #ax.relim()
    #ax.autoscale_view(True, True, True)
    #fig.canvas.draw()
    #fig.canvas.flush_events()
    #fig.tight_layout()
    #fig.show()

    # Saving the plot to a png with the same name as this file
    #fig.savefig(sim['file_name'] + '.png')



