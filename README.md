# README.md

## Mechanical Fatigue Simulation of Molecular Bonds

This code simulates the mechanical fatigue of a system of $N_B$ bonds arranged in a column, subjected to a force perpendicular to the column. The simulation tracks the breaking and potential rebinding of these bonds over cycles. The goal is to understand how different models of force distribution upon bond breakage affect the overall lifetime of the system (time until all bonds rupture).

### File Descriptions

* **`main.py`**: This is the main script that runs the simulations. It defines the simulation parameters, iterates through different stress levels (by changing the initial breaking probability), and for each stress level, runs multiple repetitions. It utilizes the `bonds` class from `Bonds.py` to model the bond system and records the number of cycles until complete failure. It also includes options to save the detailed breaking/rebinding events, the cycles to failure for each repetition, and the summary statistics (mean, median, standard deviation) at each stress level. Additionally, it has options to generate histograms of failure cycles and plots of probability vs. median failure cycles.

* **`Bonds.py`**: This file contains the `bonds` class, which encapsulates the state and behavior of the system of bonds. The class initializes the bonds, their breaking probabilities, and handles the logic for bond breaking according to the chosen model (`p_const`, `p_increase`, or `zipper`). It also implements the optional bond rebinding mechanism. The file also includes helper functions:
    * `add_bonds_to_dic`: Records bond breaking/rebinding events.
    * `plot_hist`: Generates a histogram of the number of cycles to failure.
    * `plot_p_n`: Generates a log-log plot of breaking probability vs. median cycles to failure.

### Simulation Models

The `model` parameter in `main.py` determines how the probability of bond breaking evolves during the simulation:

1.  **`p_const` (Constant Probability)**: The probability of breaking ($p_1$) is the same for all bonds and remains constant throughout the simulation, even after bonds break.
2.  **`p_increase` (Increasing Probability)**: Initially, all bonds have the same breaking probability ($p_1$). However, after each bond breaks, the total applied force is assumed to be redistributed equally among the remaining intact bonds, increasing the force on each and thus their breaking probability.
3.  **`zipper`**: Bonds break sequentially, one at a time. The entire applied force is concentrated on the first unbroken bond. Once it breaks, the force shifts to the next unbroken bond, and so on.

### Parameters in `main.py`

* `model`: The bond breaking model to use ('p\_const', 'p\_increase', 'zipper').
* `rebind_on`: A boolean indicating whether broken bonds can rebind.
* `Rebind_P1`: The probability of a broken bond rebinding in a single cycle.
* `number_of_bonds`: The total number of bonds in the system.
* `number_of_s_changes`: The number of different stress levels (or initial breaking probabilities) to simulate.
* `number_of_repetitions`: The number of simulation runs to perform at each stress level to gather statistics.
* `starting_p0`: The base probability of a bond breaking with zero external force.
* `f_change`: The step size for changing the force (used to vary the initial breaking probability).
* `f_const0`: The initial total force-related constant.
* `directory`: The directory where the results will be saved.
* `date`: A date string to include in the filenames.
* `write_bonds`: If `True`, detailed information about each bond's breaking/rebinding events is saved to a file. 
* `file_name_b`: A unique string to add to the bond event filename.
* `write_ns`: If `True`, the number of cycles until all bonds break for each repetition at each stress level is saved.
* `file_name_alln`: A unique string to add to the all-bonds-broken filename.
* `file_name_s`: A unique string to add to the summary statistics filename (mean, median, std).
* `to_plot_hist`: If `True`, a histogram of the number of cycles to failure is plotted for each stress level.
* `to_plot_p_n`: If `True`, a plot of the initial breaking probability vs. the median number of cycles to failure is generated.

### How to Run

1.  Ensure you have Python 3 installed.
2.  Save the code for the `bonds` class in a file named `Bonds.py` in the same directory as `main.py`.
3.  Run the simulation by executing `python main.py` from your terminal.

The results will be saved in the `results/<model_name>/` directory, including CSV files containing the simulation data and, optionally, plots if the `to_plot_hist` or `to_plot_p_n` flags are set to `True`.
