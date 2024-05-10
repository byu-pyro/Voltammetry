from functions import *
import subprocess

# Parameters for the physical system
# ----------------------------------------------------------------------------------------------------------------
for j in range(10):
    # salts used
    salts = np.array(["LiCl", "KCl", "UCl3"])

    # salt mass fractions
    y_salts = np.array([0.287, 0.703, 0.01])

    # molecular weight of each salt
    MW_salts = np.array([42.394, 74.551, 344.388])  # g/mol

    # density parameters of the salt based on a linear function (density at 0K, density change per degree K)
    density_parameters_salts = np.array([
        [1.88, -4.32E-04],  # LiCl
        [2.14, -5.83E-04],  # KCl
        [13.66, -79.4E-04, ],  # UCl3
    ])

    # ions in the salts
    ions = np.array(["Li+", "K-", "[UCl6]3-", "Cl-"])

    # charge of the ions
    z_ions = np.array([1, 1, -3, -1])

    # stoichiometry of the ions in the salt
    stoic_ions = np.array([
        # Li+,      K+,     [UCl6]3-,   Cl+
        [1,         0,      0,          1],  # LiCl
        [0,         1,      0,          1],  # KCl
        [0,         0,      1,          -3],  # UCl3, this has -3 for chlorine because it forms the complex
    ])

    # diffusivity parameters of the salt based on the arrheneus equation(pre-exponential factor, activation energy)
    diffusivity_parameters_ions = np.array([
        [1.65E-3, -2.516E4],  # Li+
        [1.93E-3, -2.812E4],  # K+
        [1.37E-3, -2.420E4],  # Cl-
        [5.92E-4, -2.044E4],  # [UCl6]3-
    ])

    # Temperature
    T_C = 500  # C
    T = T_C + 273.15  # K

    # The scan rate
    v = 0.02 #V/s

    # number of electrons transferred
    n = 3

    # standard potential
    Eo = 0.0

    # the index of the analyte
    j_analyte = 2

    # stoichiometry of the reduction  reaction
    stoic_rxn = np.array([0, 0, -1, 6])

    # Grid Space
    x = get_gridspace(L=0.5, n_x=70, alpha=1.1)

    # how often do you want a data point (in V)
    capture_rate = 0.0001

    # Simulation Assumptions
    # ---------------------------------------------------------------------------------------------------------------------
    ideal_solid_activity = False #assume ideal solid activity?
    ideal_ion_activity = False
    no_migration = False #assume no migration?
    no_ohmic_losses = False

    # percent of estimated ohmic resistance to account for
    iR_compensation = 0.85 #make sure no_ohmic losses is False if not 0

    # amount of deposit to have full coverage of the electrode
    No = 4e-9  # mol/cm^2

    # difference between the pure and dilute reference state
    dEo = -0.113  # V

    # distance between the working electrode and reference electrode
    L_RE = 0.25 # cm



    # Calculatin Input
    # --------------------------------------------------------------------------------------------------------


    # initial concentration of all the ions
    C_ions = get_C_ions(T, y_salts, MW_salts, density_parameters_salts, stoic_ions)

    # resistance between the working electrode and the reference electrode
    R_ohmic = get_R_ohmic(L_RE, C_ions, T, z_ions, diffusivity_parameters_ions)

    # the reference potential for the plots
    E_star = get_E_star(C_ions[j_analyte], Eo, n, T, ideal_ion_activity, dEo)

    # The starting potential for the scans
    E_start = E_star + 0.3

    # get the amount of deposit at the surface
    N = np.linspace(No / 1e8, No, 10000)  # an array of guess values
    # E at all of these N
    N_guess = N[np.argmin(np.abs(Nernst(N, C_ions[j_analyte], Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No)
                                 - E_start))]  # closest guess
    N_start = get_N(N_guess, E_start, C_ions[j_analyte], Eo, n, T, ideal_solid_activity, ideal_ion_activity, dEo, No)

    #Initial concentraiton profile
    C_start = np.ones((len(ions), len(x))) * C_ions[:, None]

    # Initialize Simulation File
    # --------------------------------------------------------------------------------------------------------------------
    sim = {
        'status': 'not started',
        'file_name': f'test{j}',
        'path': 'data/',
        'parameters': {
            'T': T,
            'salts': salts,
            'y_salts': y_salts,
            'MW_salts': MW_salts,
            'density_parameters_salts': density_parameters_salts,
            'ions': ions,
            'z_ions': z_ions,
            'stoic_ions': stoic_ions,
            'diffusivity_parameters_ions': diffusivity_parameters_ions,
            'j_analyte': j_analyte,
            'E_start': E_start,
            'v': v,
            'n': n,
            'Eo': Eo,
            'stoic_rxn': stoic_rxn,
            'x': x,
            'E_star': E_star},
        'assumptions': {
            'ideal_solid_activity': ideal_solid_activity,
            'ideal_ion_activity': ideal_ion_activity,
            'no_migration': no_migration,
            'no_ohmic_losses': no_ohmic_losses,
            'iR_compensation': iR_compensation,
            'No': No,
            'dEo': dEo,
            'L_RE': L_RE,
            'R_ohmic': R_ohmic},
        'data': {
            'capture_rate': capture_rate,
            't': np.array([0.0]),  # 1D array for time
            'E': np.array([E_start]),  # 1D array from potential
            'i': np.array([0.0]),  # 1D array for current
            'C': np.array([C_start]),  # 3D array of concentration (time, ions, distance)
            'N': np.array([N_start]),
            'compute_time': np.array([0.0])}}  # 1D array from amount of deposit

    file = open(f'{sim["path"]}{sim["file_name"]}', 'wb')
    pickle.dump(sim, file)



def run_python_script():
    try:
        # The command you want to execute
        command = "bash submision_script.sh"
        
        # Use subprocess to execute the command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        
        # Get the output
        output, error = process.communicate()
        
        # Print the output
        print('Output: ', output.decode())
        
        # Print the error if there is any
        if error:
            print('Error: ', error.decode())
            
    except Exception as e:
        print("An error occurred: ", str(e))

# Call the function
run_python_script()
