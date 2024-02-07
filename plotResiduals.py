import os
import matplotlib.pyplot as plt
import numpy as np
import re

# Function to extract residuals from the log file
def extract_residuals(log_file, variable):
    """
    Function to extract residuals from the log file.

    Args:
        log_file (str): The path to the log file.
        variable (str): The variable to extract residuals for.

    Returns:
        numpy.ndarray: An array of residuals.
    """
    if not isinstance(variable, str) or variable == '':
        raise ValueError("Variable must be a non-empty string.")

    if not os.path.isfile(log_file):
        raise FileNotFoundError(f"Log file '{log_file}' does not exist.")

    residuals = []
    pattern = re.compile(r"[-+]?\d*\.\d+|\d+")
    with open(log_file, 'r') as file:
        for line in file:
            if f'Solving for {variable}' in line:
                residual = pattern.search(line.split(',')[1])
                if residual:
                    residuals.append(float(residual.group()))
    return np.array(residuals)


def save_residuals(residuals: dict, variables: list, output_file: str) -> None:
    """
    Save residuals in a text file in a specific format.

    Args:
        residuals (dict): A dictionary containing the residuals for different variables.
                          The keys are the variable names and the values are lists of residuals.
        variables (list): A list of variable names for which the residuals are provided.
        output_file (str): The path to the output file where the residuals will be saved.

    Returns:
        None. The function saves the residuals in the specified output file.
    """
    if not os.path.isfile(output_file):
        raise FileNotFoundError(f"Output file '{output_file}' does not exist.")
    if not os.access(output_file, os.W_OK):
        raise PermissionError(f"Output file '{output_file}' is not writable.")
    
    with open(output_file, 'w') as file:
        header = f'Iteration {" ".join(variables)}\n'
        file.write(header)
        
        lines = [f"{i} {' '.join(str(residuals[variable][i]) for variable in variables)}\n" for i, value in enumerate(residuals[variables[0]])]
        file.writelines(lines)


# File path to your interFoamlog
log_file_path = '/home/joaofn/OpenFOAM/joaofn-10/run/H3_OF10_Turbulent/interFoamLog'

# Extracting residuals for different variables
variables = ['Ux', 'Uy', 'p'] #'omega', 'k', 'p']
residuals = {}
for variable in variables:
    residuals[variable] = extract_residuals(log_file_path, variable)
    output_file = 'residuals.txt'
    save_residuals(residuals, variables, output_file)
    
# Plotting
fig, axs = plt.subplots(len(variables), 1, figsize=(10, 6*(len(variables))))

for i, variable in enumerate(variables):
    ax = axs[i]
    ax.plot(residuals[variable], label=variable)
    ax.set_yscale('log')  # Use 'linear' for linear scale
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Residual')
    ax.set_title(f'Residuals for {variable}')
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()