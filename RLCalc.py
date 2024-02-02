import math
import cmath
import sys

def circuit_properties(circuit_type, Xc, Xl, Z, Y, phi, angular_f, resonant_f, quality_factor, resonant_ang_f):
    formatted_Xc = format_values_r(abs(Xc))
    formatted_xl = format_values_r(abs(Xl))
    formatted_Z = format_values_r(abs(Z))
    formatted_resonant_f = format_values_f(abs(resonant_f))
    print(f"{circuit_type} Circuit:")
    if circuit_type == "Series RC" or circuit_type == "Parallel RC":
        print("Angular frequency(\u03C9) = {} rad/s ".format(angular_f))
        print("Capacitive reactance(Xc) = {} ".format(formatted_Xc))
        print("Total Impedance = {}".format(formatted_Z))
        print("Admittance = {:.6f} siemens".format(abs(Y)))
        print("Phase difference = {:.6f}\u00B0".format(phi * 180 / cmath.pi))
        print("Phase difference = {:.6f} rad".format(phi))
    elif circuit_type == "Series RL" or circuit_type == "Parallel RL":
        print("Angular frequency(\u03C9) = {} rad/s ".format(angular_f))
        print("Inductive reactance(Xl) = {}".format(formatted_xl))
        print("Total Impedance = {}".format(formatted_Z))
        print("Admittance = {:.6f} siemens".format(abs(Y)))
        print("Phase difference = {:.6f}\u00B0".format(phi * 180 / cmath.pi))
        print("Phase difference = {:.6f} rad".format(phi))
    elif circuit_type == "Series RLC" or circuit_type == "Parallel RLC":
        print("Angular frequency(\u03C9) = {} rad/s ".format(angular_f))
        print("Capacitive reactance(Xc) = {} ".format(formatted_Xc))
        print("Inductive reactance(Xl) = {}".format(formatted_xl))
        print("Total Impedance = {}".format(formatted_Z))
        print("Admittance = {:.6f} siemens".format(abs(Y)))
        print("Phase difference = {:.6f}\u00B0".format(phi * 180 / cmath.pi))
        print("Phase difference = {:.6f} rad".format(phi))
        print("Quality factor: {}".format(quality_factor))
        print("Resonant frequency = {} ".format(formatted_resonant_f))
        print("Resonant angular frequency = {:.6f} rad/s".format(resonant_ang_f))
    else:
        print("If you see this David made an oopsie")
        sys.exit()

calculate_omega = lambda f:2*cmath.pi*f
calculate_Xc = lambda C, f:1/((calculate_omega(f))*C)
calculate_Xl = lambda L, f:(calculate_omega(f))*L
calculate_impedance_S_RC = lambda R, Xc: cmath.sqrt((R ** 2) + (1 / (2 * cmath.pi * f * C)**2))
calculate_impedance_P_RC = lambda R, Xc: R*Xc/cmath.sqrt(R**2 + Xc**2)
calculate_impedance_S_RL = lambda R, Xl: cmath.sqrt(R ** 2 + Xl ** 2)
calculate_impedance_P_RL = lambda R, Xl: cmath.sqrt(1 / ((1 / R ** 2) + (1 / Xl ** 2)))
calculate_impedance_P_RLC = lambda R, Xc, Xl: 1 / cmath.sqrt((1 / R ** 2) + ((1 / Xl) - (1 / Xc))**2)
calculate_impedance_S_RLC = lambda R, Xc, Xl: cmath.sqrt(R **2 + (Xl - Xc)**2)
calculate_admittance = lambda Z: 1/Z
calculate_phase_angle_S_RC = lambda R, C,: math.atan(-(1 / (R * (calculate_omega(f)) * C)))
calculate_phase_angle_P_RC = lambda R, C,: math.atan(-(calculate_omega(f) * C * R))
calculate_phase_angle_S_RL = lambda R, L : math.atan(calculate_omega(f) * L / R)
calculate_phase_angle_P_RL = lambda R, L : math.atan(R / (calculate_omega(f) * L))
calculate_phase_angle_P_RLC = lambda R, Xc, Xl : math.atan((R / Xl) - (R / Xc))
calculate_phase_angle_S_RLC = lambda R, Xc, Xl: math.atan((Xl - Xc) / R)
calculate_AngularFrequency = lambda f : 2*cmath.pi*f
calculate_resonant_frequency = lambda L, C : 1 / (2*cmath.pi*math.sqrt(L*C))
calculate_quality_factor_S_RLC = lambda R, L, C : (1/R)*math.sqrt(L/C)
calculate_quality_factor_P_RLC = lambda R, L, C : R*math.sqrt(C/L)
calculate_resonant_angular_frequency = lambda L, C : 1 /math.sqrt(L*C)

def parallel_rc_calc(R, C, f,):
    angular_f = calculate_AngularFrequency(f)
    Xc = calculate_Xc(C, f)
    Z = calculate_impedance_P_RC(R, Xc)
    Y = calculate_admittance(Z)
    phi = calculate_phase_angle_P_RC(R, C)
    circuit_properties("Parallel RC", Xc, 0, Z, Y, phi, angular_f, 0, 0, 0)
def series_rc_calc(R, C, f):
    angular_f = calculate_AngularFrequency(f)
    Xc = calculate_Xc(C, f)
    Z = calculate_impedance_S_RC(R, Xc)
    Y = calculate_admittance(Z)
    phi = calculate_phase_angle_S_RC(R, C)
    circuit_properties("Series RC", Xc, 0, Z, Y, phi, angular_f, 0, 0, 0)
def series_rl_calc(R, L, f):
    angular_f = calculate_AngularFrequency(f)
    Xl = calculate_Xl(L, f)
    Z = calculate_impedance_S_RL(R, Xl)
    Y = calculate_admittance(Z)
    phi = calculate_phase_angle_S_RL(R, L)
    circuit_properties("Series RL", 0, Xl, Z, Y, phi, angular_f, 0, 0, 0)
def parallel_rl_calc(R, L, f):
    angular_f = calculate_AngularFrequency(f)
    Xl = calculate_Xl(L, f)
    Z = calculate_impedance_P_RL(R, Xl)
    Y = calculate_admittance(Z)
    phi = calculate_phase_angle_P_RL(R, L)
    circuit_properties("Parallel RL", 0, Xl, Z, Y, phi, angular_f, 0, 0, 0)
def parallel_rlc_calc(R, C, f, L):
    angular_f = calculate_AngularFrequency(f)
    Xc = calculate_Xc(C, f)
    Xl = calculate_Xl(L, f)
    Z = calculate_impedance_P_RLC(R, Xc, Xl)
    Y = calculate_admittance(Z)
    phi = calculate_phase_angle_P_RLC(R, Xc, Xl)
    resonant_f = calculate_resonant_frequency(L, C)
    resonant_angular_f = calculate_resonant_angular_frequency(L, C)
    q_factor = calculate_quality_factor_P_RLC(R, L, C)
    circuit_properties("Parallel RLC", Xc, Xl, Z, Y, phi, angular_f, resonant_f, q_factor, resonant_angular_f)

def series_rlc_calc(R, C, f, L):
    angular_f = calculate_AngularFrequency(f)
    Xc = calculate_Xc(C, f)
    Xl = calculate_Xl(L, f)
    Z = calculate_impedance_S_RLC(R, Xc, Xl)
    Y = calculate_admittance(Z)
    phi = calculate_phase_angle_S_RLC(R, Xc, Xl)
    resonant_f = calculate_resonant_frequency(L, C)
    resonant_angular_f = calculate_resonant_angular_frequency(L, C)
    q_factor = calculate_quality_factor_S_RLC(R, L, C)
    circuit_properties("Series RLC", Xc, Xl, Z, Y, phi, angular_f, resonant_f, q_factor, resonant_angular_f)

# Gets a valid unit
def get_valid_unit(prompt, valid_units):
    while True:
        user_input = input(prompt)
        check_quit_input(user_input)
        if user_input in valid_units:
            return user_input
        else:
            print("Invalid unit!")
# Gets valid value
def get_valid_value(prompt):
    while True:
        user_input = input(prompt)
        check_quit_input(user_input)
        try:
            value = float(user_input)
            if value > 0:
                return value
            else:
                print("Please enter a value greater than 0.")
        except ValueError:
            print("Invalid Value!")

# Converts_units
def convert_units(x, unit):
    if unit in ["Ohm", "Hz", "F", "H"]:
        return x
    elif unit in ["kOhm", "kHz"]:
        return x * 1000.0
    elif unit in ["MOhm", "MHz"]:
        return x * 1000000.0
    elif unit == "GHz":
        return x * 1000000000.0
    elif unit in ["mOhm", "mHz", "mF", "mH"]:
        return x * 0.001
    elif unit in ["uF", "uH"]:
        return x * 0.000001
    elif unit == "nF":
        return x * 0.000000001
    elif unit == "pF":
        return x * 0.000000000001

def format_values_r(x):
    if x >= 1e6:
        return "{:.6f} M\u03A9".format(x / 1e6)
    elif x >= 1e3:
        return "{:.6f} k\u03A9".format(x / 1e3)
    elif x < 0.01:
        return "{:.6f} m\u03A9".format(x / 0.001)
    else:
        return "{:.6f} \u03A9".format(x)

def format_values_f(x):
    if x >= 1e6:
        return "{:.6f} MHz".format(x / 1e6)
    elif x >= 1e3:
        return "{:.6f} kHz".format(x / 1e3)
    elif x < 0.01:
        return "{:.6f} mHz".format(x / 0.001)
    else:
        return "{:.6f} Hz".format(x)
# quit program function
def check_quit_input(x):
    while True:
        if x in quit_commands:
            sys.exit()
        else:
            break
# Valid units
valid_unit_R = ["mOhm", "Ohm", "kOhm", "MOhm"]
valid_unit_f = ["mHz", "Hz", "kHz", "MHz", "GHz"]
valid_unit_L = ["H", "mH", "uH"]
valid_unit_C = ["F", "mF", "uF", "nF", "pF"]
valid_circuit_type = ["RL", "RC", "RLC"]
valid_circuit_kind = ["S", "P"]
quit_commands = ["quit", "q", "exit", "QUIT", "Q", "EXIT"]

# Quiz
while True:
    circuit_type = input("Select circuit type (RL, RC, or RLC): ")
    if check_quit_input(circuit_type):
        break
    elif circuit_type.upper() not in valid_circuit_type:
        print("Invalid circuit type. Please select RL, RC, or RLC")
    else:
        break

while True:
    circuit_kind = input("Select circuit type (P for parallel and S for serial): ")
    if check_quit_input(circuit_kind):
        break
    elif circuit_kind.upper() not in valid_circuit_kind:
        print("Invalid circuit type. ")
    else:
        break

    # Series_RC_circuit_values
if circuit_type.upper() == "RC" and circuit_kind.upper() == "S":
    # Resistence input
    R_unit = get_valid_unit("Enter the resistence unit (mOhm, Ohm, kOhm, MOhm): ", valid_unit_R)
    R = get_valid_value("Enter a resistence value: ")
    R = convert_units(R, R_unit)
    # Frequency input
    f_unit = get_valid_unit("Enter frequency unit (mHz, Hz, kHz, MHz, GHz): ", valid_unit_f)
    f = get_valid_value("Enter frequency value: ")
    f = convert_units(f, f_unit)
    # Capacitance input
    C_unit = get_valid_unit("Enter a capacitance unit (F, mF, uF, nF, pF): ", valid_unit_C)
    C = get_valid_value("Enter a capacitance value: ")
    C = convert_units(C, C_unit)
    series_rc_calc(R, C, f)

# Parallel_RC_circuit_values
elif circuit_type.upper() == "RC" and circuit_kind.upper() == "P":
    # Resistence input
    R_unit = get_valid_unit("Enter the resistence unit (mOhm, Ohm, kOhm, MOhm): ", valid_unit_R)
    R = get_valid_value("Enter a resistence value: ")
    R = convert_units(R, R_unit)
    # Frequency input
    f_unit = get_valid_unit("Enter frequency unit (mHz, Hz, kHz, MHz, GHz): ", valid_unit_f)
    f = get_valid_value("Enter frequency value: ")
    f = convert_units(f, f_unit)
    # Capacitance input
    C_unit = get_valid_unit("Enter a capacitance unit (F, mF, uF, nF, pF): ", valid_unit_C)
    C = get_valid_value("Enter a capacitance value: ")
    C = convert_units(C, C_unit)
    parallel_rc_calc(R, C, f)

    # Series_RL_circuit_values
elif circuit_type.upper() == "RL" and circuit_kind.upper() == "S":
    # Resistence input
    R_unit = get_valid_unit("Enter the resistence unit (mOhm, Ohm, kOhm, MOhm): ", valid_unit_R)
    R = get_valid_value("Enter a resistence value: ")
    R = convert_units(R, R_unit)
    # Frequency input
    f_unit = get_valid_unit("Enter frequency unit (mHz, Hz, kHz, MHz, GHz): ", valid_unit_f)
    f = get_valid_value("Enter frequency value: ")
    f = convert_units(f, f_unit)
    # Inductance input
    L_unit = get_valid_unit("Enter Inductance unit (H, mH, uH): ", valid_unit_L)
    L = get_valid_value("Enter inductance value: ")
    L = convert_units(L, L_unit)
    series_rl_calc(R, L, f)

    # Parallel_RL_circuit_values
elif circuit_type.upper() == "RL" and circuit_kind.upper() == "P":
    # Resistence input
    R_unit = get_valid_unit("Enter the resistence unit (mOhm, Ohm, kOhm, MOhm): ", valid_unit_R)
    R = get_valid_value("Enter a resistence value: ")
    R = convert_units(R, R_unit)
    # Frequency input
    f_unit = get_valid_unit("Enter frequency unit (mHz, Hz, kHz, MHz, GHz): ", valid_unit_f)
    f = get_valid_value("Enter frequency value: ")
    f = convert_units(f, f_unit)
    # Inductance input
    L_unit = get_valid_unit("Enter Inductance unit (H, mH, uH): ", valid_unit_L)
    L = get_valid_value("Enter inductance value: ")
    L = convert_units(L, L_unit)
    parallel_rl_calc(R, L, f)

    # Series_RLC_circuit_values
elif circuit_type.upper() == "RLC" and circuit_kind.upper() == "S":
    # Resistence input
    R_unit = get_valid_unit("Enter the resistence unit (mOhm, Ohm, kOhm, MOhm): ", valid_unit_R)
    R = get_valid_value("Enter a resistence value: ")
    R = convert_units(R, R_unit)
    # Frequency input
    f_unit = get_valid_unit("Enter frequency unit (mHz, Hz, kHz, MHz, GHz): ", valid_unit_f)
    f = get_valid_value("Enter frequency value: ")
    f = convert_units(f, f_unit)
    # Inductance input
    L_unit = get_valid_unit("Enter Inductance unit (H, mH, uH): ", valid_unit_L)
    L = get_valid_value("Enter inductance value: ")
    L = convert_units(L, L_unit)
    # Capacitance input
    C_unit = get_valid_unit("Enter a capacitance unit (F, mF, uF, nF, pF): ", valid_unit_C)
    C = get_valid_value("Enter a capacitance value: ")
    C = convert_units(C, C_unit)
    series_rlc_calc(R, C, f, L)

    # Parallel_RLC_Circuit_values
else:
    # Resistence input
    R_unit = get_valid_unit("Enter the resistence unit (mOhm, Ohm, kOhm, MOhm): ", valid_unit_R)
    R = get_valid_value("Enter a resistence value: ")
    R = convert_units(R, R_unit)
    # Frequency input
    f_unit = get_valid_unit("Enter frequency unit (mHz, Hz, kHz, MHz, GHz): ", valid_unit_f)
    f = get_valid_value("Enter frequency value: ")
    f = convert_units(f, f_unit)
    # Inductance input
    L_unit = get_valid_unit("Enter Inductance unit (H, mH, uH): ", valid_unit_L)
    L = get_valid_value("Enter inductance value: ")
    L = convert_units(L, L_unit)
    # Capacitance input
    C_unit = get_valid_unit("Enter a capacitance unit (F, mF, uF, nF, pF): ", valid_unit_C)
    C = get_valid_value("Enter a capacitance value: ")
    C = convert_units(C, C_unit)
    parallel_rlc_calc(R, C, f, L)