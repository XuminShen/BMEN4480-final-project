from keras import Model
from sklearn.preprocessing import MinMaxScaler
import numpy as np
from sklearn.model_selection import GridSearchCV, train_test_split
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from pandas import read_csv
import os
import warnings
import tensorflow as tf
from keras import layers
from keras import losses
from keras import Sequential
from keras.utils.np_utils import to_categorical
import keras
from keras.optimizers import Adam
import random
from collections import Counter

file=pd.read_csv("top8antigens.csv")
file2=pd.read_csv("1600_informations.csv")

df = pd.DataFrame({'x': file['cdr3.beta'], 'v':file['v.beta'] ,'d':file['d.beta'],'j':file['j.beta'],'y': file['antigen.epitope']})
df2 = pd.DataFrame({'x': file2['cdr3.beta'], 'v':file2['v.beta'] ,'d':file2['d.beta'],'j':file2['j.beta'],'y': file2['antigen.epitope']})

letters = set(''.join(df['x']))
print(letters)
num_letters = len(letters)

# Print the result
print('The number of unique letters in the data is:', num_letters)
cdr3_len =[]
for i in df['x']:
  cdr3_len.append(len(i))
freq_dict = Counter(cdr3_len)
#
# # Create a histogram of the data
# plt.hist(cdr3_len, bins=100)
#
# # Set the title and axis labels
# plt.title('Frequency of cdr3.beta sequence length')
# plt.xlabel('cdr3.beta sequence length')
# plt.ylabel('Frequency')

# # Add the frequency on top of each bar
# for bin, freq in freq_dict.items():
#     plt.text(bin, freq, str(freq), ha='center', va='bottom')

# Show the plot
# plt.show()


# ONE-HOT
data = df['x']
data2 = df2['x']

# Define a dictionary to map characters to integers
char_to_int = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}

# Define the maximum length of the encoded sequences
max_len = max([len(s) for s in data])
max_len2 = max([len(s) for s in data2])

# Initialize an empty matrix to store the encoded sequences
encoded_data = np.zeros((len(data2), max_len, len(char_to_int)), dtype=np.int32)

# Encode each sequence using one-hot encoding and store in the matrix
for i, s in enumerate(data2):
    for j, c in enumerate(s):
        encoded_data[i, j, char_to_int[c]] = 1

# Print the encoded data
print(encoded_data[0].shape)
encoded_data_swapped = np.transpose(encoded_data, (0, 2, 1))

# Print the shape of the swapped encoded data
print(encoded_data_swapped.shape)

# Save the swapped encoded data to a file
#np.save('cdr3_onehot_swapped.npy', encoded_data_swapped)
np.save('cdr3_onehot_1600.npy',encoded_data)

unique_v = set(df['v'])
string_to_index = {string: i for i, string in enumerate(unique_v)}
num_strings = len(unique_v)
num_samples = len(df2['v'])

# Initialize the one-hot encoded matrix
one_hot_v = np.zeros((num_samples, num_strings))

# Fill in the matrix with ones where the corresponding string is present
for i, string in enumerate(df2['v']):
    if string in string_to_index:
        index = string_to_index[string]
        one_hot_v[i, index] = 1
one_hot_v[one_hot_v == 0] = np.nan
one_hot_v[np.isnan(one_hot_v)] = 0

unique_d = set(df['d'])
string_to_index = {string: i for i, string in enumerate(unique_d)}
num_strings = len(unique_d)
num_samples = len(df2['d'])


# Initialize the one-hot encoded matrix
one_hot_j = np.zeros((num_samples, num_strings))

# Fill in the matrix with ones where the corresponding string is present
for i, string in enumerate(df2['d']):
    if string in string_to_index:
        index = string_to_index[string]
        one_hot_j[i, index] = 1
one_hot_j[one_hot_j == 0] = np.nan
one_hot_j[np.isnan(one_hot_j)] = 0

combined_matrix = np.hstack((one_hot_v, one_hot_j))
print(combined_matrix.shape)
#combined_matrix = np.hstack((one_hot_v, one_hot_j))
np.save('vdj_onehot_1600.npy', combined_matrix)
