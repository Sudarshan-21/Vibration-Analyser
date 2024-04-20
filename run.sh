#!/bin/bash
echo "Enter stiffness constant (k) in N/m:"
read stiffness_constant

# Prompt user to enter mass (m) in kg
echo "Enter mass (m) in kg:"
read mass

# Calculate natural frequency using the formula sqrt(k/m)
natural_frequency=$(echo "scale=6; sqrt($stiffness_constant / $mass)" | bc)

echo "Natural Frequency (ω₀) = $natural_frequency rad/s"

echo "Enter initial displacement (x0) in meters:"
read x0

# Prompt user to enter initial velocity (v0) in m/s
echo "Enter initial velocity (v0) in m/s:"
read v0

# Prompt user to select a vibration type
echo "Enter your choice (1.free_undamped, 2.free_damped, 3.forced_undamped, or 4.forced_damped):"
read choice

# Determine the vibration type based on user choice
case $choice in
    1)
        vibration_type="free_undamped"
        ;;
    2)
        vibration_type="free_damped"
        ;;
    3)
        vibration_type="forced_undamped"
        ;;
    4)
        vibration_type="forced_damped"
        ;;
    *)
        echo "Invalid choice. Please enter a valid input"
        exit 1
        ;;
esac

# Set default values for duration and sampling rate
duration=2.0
sampling_rate=1000

# Prompt for damping coefficient and validate if necessary
while true; do
    if [ "$vibration_type" == "free_damped" ] || [ "$vibration_type" == "forced_damped" ]; then
        echo "Enter damping coefficient (must be positive):"
        read damping_coeff
        if (( $(echo "$damping_coeff >= 0" | bc -l) )); then
            break  # Valid positive value, exit loop
        else
            echo "Invalid damping coefficient. Please enter a positive value."
        fi
    else
        # Default damping coefficient for non-damped vibration types
        damping_coeff=0.0
        break  # No damping coefficient needed, exit loop
    fi
done

if [ "$vibration_type" == "forced_undamped" ] || [ "$vibration_type" == "forced_damped" ]; then
    echo "Enter external force (F0) in N:"
    read F0
    echo "Enter frequency of external force (omega) in rad/s:"
    read omega
else
    # Default external force value for non-forced vibration types
    F0=0.0
    omega=0.0
fi

echo "External force (F0): $F0"
echo "Frequency of external force (omega): $omega"


# Run the Python script with selected parameters
python3 Vibration_main.py --vibration_type "$vibration_type" --stiffness_constant $stiffness_constant --mass $mass --x0 $x0 --v0 $v0 --damping_coeff $damping_coeff --F0 $F0 --omega $omega