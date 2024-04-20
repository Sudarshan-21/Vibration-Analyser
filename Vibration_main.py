import numpy as np
import matplotlib.pyplot as plt
import argparse

def free_undamped_waveform(time_array, x0, v0, natural_frequency):
    """Generate a natural undamped vibration waveform."""
    # Calculate R and phi
    R = np.sqrt(x0**2 + (v0 / natural_frequency)**2)
    phi = np.arctan((x0 * natural_frequency) / v0)
    
    # Calculate waveform
    waveform = R * np.sin(natural_frequency * time_array + phi)
    return waveform

def free_damped_waveform(time_array, x0, v0, natural_frequency, damping_ratio):
    """Generate a damped vibration waveform."""
    # Calculate damping factor
    
    if damping_ratio < 1:
        # Underdamped vibration
        omega_d = natural_frequency * np.sqrt(1 - damping_ratio**2)
        A = np.sqrt(x0**2 * natural_frequency**2 + v0**2 + 2 * damping_ratio * natural_frequency * v0 * x0)
        phi = np.arctan((v0 + damping_ratio * natural_frequency * x0) / (omega_d * x0))
        
        waveform = A * np.cos(omega_d * time_array - phi) * np.exp(-damping_ratio * natural_frequency * time_array)
    
    elif damping_ratio == 1:
        # Critically damped vibration
        A1 = x0
        A2 = v0 + x0 * natural_frequency
        waveform = (A1 + A2 * time_array) * np.exp(-natural_frequency * time_array)
    
    else:
        # Overdamped vibration
        omega_d = natural_frequency * np.sqrt(damping_ratio**2 - 1)
        alpha1 = -damping_ratio * natural_frequency + omega_d
        alpha2 = -damping_ratio * natural_frequency - omega_d
        A1 = (v0 - alpha2 * x0) / (alpha1 - alpha2)
        A2 = (alpha1 * x0 - v0) / (alpha1 - alpha2)
        waveform = A1 * np.exp(alpha1 * time_array) + A2 * np.exp(alpha2 * time_array)
    
    # Calculate waveform
    return waveform

def forced_undamped_waveform(time_array, natural_frequency, F0, stiffness_coefficient, omega):
    """Generate a forced undamped vibration waveform."""
    print(f"F0: {F0}")
    print(f"Stiff_coeff: {stiffness_coefficient}")
    print(f"nat freq: {natural_frequency}")
    print(f"omega: {omega}")
    P = (F0 / stiffness_coefficient) / (1 - (omega / natural_frequency)**2)
    
    waveform = P * np.cos(omega * time_array )
    return waveform

def forced_damped_waveform(time_array, natural_frequency, damping_ratio, F0, stiffness_constant, damping_coefficient, omega, mass):
    """Generate a forced damped vibration waveform."""
   # Calculate phase angle (phi)
    phi = np.arctan((damping_coefficient * omega) / (stiffness_constant - mass * omega**2))

    # Calculate amplitude (R)
    R = (F0 / stiffness_constant) / np.sqrt((1 - (omega / natural_frequency)**2)**2 + (2 * damping_ratio * (omega / natural_frequency))**2)

    # Calculate waveform
    waveform = R * np.cos(omega * time_array )
    return waveform

def plot_waveform(time_array, waveform, vibration_type):
    """
    Plot the waveform.

    Parameters:
    - time_array: Array of time values.
    - waveform: Array of waveform values.
    - vibration_type: Type of vibration ('natural_undamped', 'natural_damped', 'forced_undamped', 'forced_damped').
    """
    plt.figure(figsize=(10, 4))
    plt.plot(time_array, waveform)
    plt.title(f'{vibration_type.capitalize()} Vibration Waveform')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Amplitude')
    plt.grid(True)
    plt.show()

def main(args):
    # Set parameters
    duration = args.duration
    sampling_rate = args.sampling_rate
    num_samples = int(duration * sampling_rate)
    time_array = np.linspace(0, duration, num_samples)

    # Calculate natural frequency from provided stiffness (k) and mass (m)
    natural_frequency = np.sqrt(args.stiffness_constant / args.mass)

    # Determine damping ratio (zeta) if damping coefficient (C) > 0
    if args.damping_coeff > 0:
        critical_damping_coefficient = 2 * np.sqrt(args.stiffness_constant * args.mass)
        damping_ratio = args.damping_coeff / critical_damping_coefficient
    else:
        damping_ratio = 0.0  # Default damping ratio for undamped cases

    # Generate waveform based on the selected vibration type
    if args.vibration_type == 'free_undamped':
        waveform = free_undamped_waveform(time_array, args.x0, args.v0, natural_frequency)
    elif args.vibration_type == 'free_damped':
        waveform = free_damped_waveform(time_array, args.x0, args.v0, natural_frequency, damping_ratio, )
    elif args.vibration_type == 'forced_undamped':
        waveform = forced_undamped_waveform(time_array, natural_frequency, args.F0, args.stiffness_constant, args.omega)
    elif args.vibration_type == 'forced_damped':
        waveform = forced_damped_waveform(time_array, natural_frequency, damping_ratio, args.F0, args.stiffness_constant, args.damping_coeff, args.omega, args.mass)
    else:
        print("Invalid vibration type. This script supports 'natural_undamped', 'natural_damped', 'forced_undamped', or 'forced_damped' vibrations.")
        return

    # Plot waveform
    plot_waveform(time_array, waveform, args.vibration_type)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Simulate vibration waveform.')
    parser.add_argument('--vibration_type', type=str, default='natural_undamped',
                        help='Type of vibration (default: natural_undamped)')
    parser.add_argument('--stiffness_constant', type=float, required=True,
                        help='Stiffness constant (k) in N/m')
    parser.add_argument('--mass', type=float, required=True,
                        help='Mass (m) in kg')
    parser.add_argument('--x0', type=float, required=True,
                        help='Initial displacement (x0) in meters')
    parser.add_argument('--v0', type=float, required=True,
                        help='Initial velocity (v0) in m/s')
    parser.add_argument('--damping_coeff', type=float, default=0.0,
                        help='Damping coefficient (C) in Ns/m (default: 0.0)')
    parser.add_argument('--F0', type=float, default=0.0,
                        help='External force (F0) in N (default: 0.0)')
    parser.add_argument('--omega', type=float, default=0.0,
                        help='Frequency of External force (F0) in rad/s (default: 0.0)')
    parser.add_argument('--duration', type=float, default=10.0,
                        help='Duration of the waveform in seconds (default: 5.0)')
    parser.add_argument('--sampling_rate', type=int, default=1000,
                        help='Sampling rate of the waveform in Hz (default: 1000)')
    args = parser.parse_args()

    main(args)
