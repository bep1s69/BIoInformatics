import numpy as np
import matplotlib.pyplot as plt
import os

def visualize_and_save(matrix, trajectory, output_folder="assignment_results"):
    """
    Plots the matrix as a heatmap and the vector evolution as a line graph.
    Saves the figures to the specified folder.
    """
    # 1. Create the Output Folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created folder: {output_folder}")

    # --- Plot 1: The Matrix Heatmap ---
    plt.figure(figsize=(8, 6))
    # imshow displays data as an image, suitable for matrices
    plt.imshow(matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Value')
    plt.title("Matrix Visualization (Heatmap)")
    plt.xlabel("Target State (j)")
    plt.ylabel("Source State (i)")
    
    # Save the matrix plot
    matrix_path = os.path.join(output_folder, "matrix_heatmap.png")
    plt.savefig(matrix_path)
    plt.close() # Close to free up memory
    print(f"Saved matrix heatmap to: {matrix_path}")

    # --- Plot 2: The Prediction Trajectory ---
    # Convert trajectory list to a 2D array for easier plotting (Steps x States)
    traj_array = np.array(trajectory)
    steps = range(len(trajectory))

    plt.figure(figsize=(10, 6))
    # Plot each component of the vector as a separate line
    for i in range(traj_array.shape[1]):
        plt.plot(steps, traj_array[:, i], marker='o', label=f'State {i+1}')

    plt.title("Vector Prediction over 5 Steps")
    plt.xlabel("Step")
    plt.ylabel("Value")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)

    # Save the trajectory plot
    traj_path = os.path.join(output_folder, "prediction_trajectory.png")
    plt.savefig(traj_path)
    plt.close()
    print(f"Saved trajectory plot to: {traj_path}")


def run_prediction(matrix, initial_vector, steps=5):
    """
    Performs the prediction and calls the visualization function.
    """
    M = np.array(matrix)
    v = np.array(initial_vector)

    # Validation
    rows, cols = M.shape
    vec_len = v.shape[0]
    if rows != cols or cols != vec_len:
        raise ValueError(f"Dimension mismatch: Matrix {rows}x{cols}, Vector {vec_len}")

    print(f"--- Starting Prediction ---")
    
    # Store history of the vector
    trajectory = [v]
    current_v = v

    # Iteration loop
    for i in range(1, steps + 1):
        current_v = np.dot(M, current_v)
        trajectory.append(current_v)
        print(f"Step {i}: {np.round(current_v, 4)}")

    # Call the visualization function
    visualize_and_save(M, trajectory)

    return trajectory

# --- Example Usage ---

# Example Matrix (3x3)
example_matrix = [
    [0.9, 0.2, 0.1],
    [0.1, 0.6, 0.1],
    [0.0, 0.2, 0.8]
]

# Example Initial Vector
example_vector = [50, 20, 30] 

# Run the software
run_prediction(example_matrix, example_vector, steps=5)