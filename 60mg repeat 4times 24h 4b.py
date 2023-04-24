import numpy as np
import matplotlib.pyplot as plt

# Parameters
DOSE = 60  # mg
BW = 60  # kg
t_max = 2  # h
duration = 16  # h
aftereffect = 26  # h
interval = 4  # h
onset = 0.05
off = 4  # h
pH = 7.4
F = 0.95
Vd = 4

# Functions
def calculate_concentration_curves(DOSE, BW, t_max, duration, aftereffect, onset, off, pH, F, Vd):
    ka = (F * DOSE / Vd) / (1 + 10 ** ((6.56 - pH) * np.log10(np.exp(1)) / np.log10(10)) * np.exp(-1 * (F * DOSE / Vd) * onset / BW))
    ke = np.log(2) / (duration + aftereffect)
    t = np.arange(0, 72, 0.1)
    Cp = (ka / (ka - ke)) * (np.exp(-1 * ke * t) - np.exp(-1 * ka * t)) * DOSE / Vd
    Cbrain = Cp / 3.35
    occupancy = Cbrain ** 0.6 / (Cbrain ** 0.6 + 10 ** 0.6) / 0.5 * 100
    return t, Cp, occupancy

def plot_concentration_and_occupancy(DOSE, BW, t_max, duration, aftereffect, onset, off, pH, F, Vd):
    t, Cp, occupancy = calculate_concentration_curves(DOSE, BW, t_max, duration, aftereffect, onset, off, pH, F, Vd)
    
    fig, ax1 = plt.subplots()
    ax1.plot(t, Cp, label='Blood Concentration', color='blue')
    ax1.set_ylabel('Blood Concentration (mg/L)')
    ax1.set_ylim(0, 120)

    ax2 = ax1.twinx()
    ax2.plot(t, occupancy, label='DAT Occupancy', color='red')
    ax2.axhline(y=50, color='black', linestyle='--', label='50% Occupancy')
    ax2.set_ylabel('DAT Occupancy (%)')
    ax2.set_ylim(0, 100)

    ax1.set_xlabel('Time (h)')
    ax1.set_xlim(0, 72)
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')

    plt.title(f"Blood Concentration and DAT Occupancy with Dose: {DOSE} mg")
    plt.tight_layout()
    plt.show()

# Main function
def main():
    print("4'-Br-4-methylaminorex Concentration and DAT Occupancy Calculator")
    print(f"Parameters are preset with the following values:\nDose: {DOSE} mg\nDosing interval: {interval} h")

    plot_concentration_and_occupancy(DOSE, BW, t_max, duration, aftereffect, onset, off, pH, F, Vd)

if __name__ == "__main__":
    main()
