import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

def generate_gabor(sigma=1, spatial_freq=0.2, phase=0, orientation=0, size=100, aspect_ratio=1.0):
    """
    Generate a Gabor patch stimulus.
    """
    # Create a grid of (x, y) coordinates
    x = np.linspace(-size / 2, size / 2, size)
    y = np.linspace(-size / 2, size / 2, size)
    X, Y = np.meshgrid(x, y)

    # Rotate the grid by the orientation angle
    theta = np.deg2rad(orientation)
    X_rot = X * np.cos(theta) + Y * np.sin(theta)
    Y_rot = -X * np.sin(theta) + Y * np.cos(theta)

    # Create the Gaussian envelope
    gaussian = np.exp(-(X_rot ** 2 + (Y_rot ** 2) * aspect_ratio ** 2) / (2 * sigma ** 2))

    # Create the sinusoidal wave
    sinusoid = np.cos(2 * np.pi * spatial_freq * X_rot + np.deg2rad(phase))

    # Combine them to get the Gabor patch
    gabor = gaussian * sinusoid

    return gabor

def plot_gabor_with_transparent_background(gabor, filename='gabor_transparent.pdf'):
    """
    Plot and save a black and white Gabor patch with a transparent background to a PDF file.

    Parameters:
    - gabor: 2D array, the Gabor patch to plot.
    - filename: string, the name of the output PDF file.
    """
    size = gabor.shape[0]

    # Normalize the Gabor patch to range between -1 and 1
    gabor_normalized = (gabor - np.min(gabor)) / (np.max(gabor) - np.min(gabor))

    # Create an RGBA image where the Gabor patch will be black-and-white
    gabor_rgba = np.zeros((size, size, 4))

    # Set RGB channels: 0 is black, 1 is white
    gabor_rgba[..., 0] = gabor_normalized  # Red channel (grayscale)
    gabor_rgba[..., 1] = gabor_normalized  # Green channel (grayscale)
    gabor_rgba[..., 2] = gabor_normalized  # Blue channel (grayscale)

    # Set alpha channel: 1 is opaque (Gabor), 0 is transparent (background)
    gabor_rgba[..., 3] = (gabor_normalized > 0) * 1.0  # Opaque Gabor, transparent background

    # Create a figure and save it to PDF with transparent background
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(gabor_rgba, extent=(-1, 1, -1, 1), interpolation='nearest')
    ax.axis('off')

    # Save the figure to PDF with a transparent background
    pdf = PdfPages(filename)
    pdf.savefig(fig, transparent=True)
    pdf.close()
    plt.close()

os.chdir('../../../Presentations/Posters/Clipart')

# Example usage
gabor = generate_gabor(sigma=50, spatial_freq=0.012, phase=0, orientation=0, size=256, aspect_ratio=1)
plot_gabor_with_transparent_background(gabor, filename='gabor_transparent.pdf')
