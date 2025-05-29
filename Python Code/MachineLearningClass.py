import numpy as np

##### Part a #####
# Load in training data and labels
face_data_dict = np.load("face_emotion_data.npz")
feat = face_data_dict["X"]
labels = face_data_dict["y"]
n, p = feat.shape

# Solve the least-squares solution. weights is the array of weight coefficients
weights = np.linalg.inv(feat.T@feat)@feat.T@labels

print(f"Part_4a._Found_weights:\n{weights}")

