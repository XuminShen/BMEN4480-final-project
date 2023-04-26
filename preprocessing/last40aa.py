import pandas as pd
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# Set the path to the folder containing the CSV files
path = "C:/Users/xumin/Desktop/TCR/1600"

# Get a list of all CSV files in the folder
all_files = glob.glob(path + "/*.csv")

# Create an empty list to store the number of columns for each file
num_columns_list = []

# Loop through each CSV file, read its header row to get the number of columns, and append it to the list
files = []
for filename in all_files:
    # with open(filename, 'r') as f:
    #     num_columns = len(f.readline().split(','))
    file_i = np.asarray(pd.read_csv(filename,index_col=0))
    num_columns_list.append(file_i.shape[0])
    files.append(file_i)

# Create a histogram of the number of columns
# plt.hist(num_columns_list, bins=range(min(num_columns_list), max(num_columns_list) + 2, 1))
# plt.xlabel("Number of columns")
# plt.ylabel("Frequency")
# plt.title("Distribution of number of beta chains")
# plt.show()

# indexs = np.arange(len(num_columns_list))
matrixs_last_40 = np.zeros([len(num_columns_list),40,40])
for i in np.arange(len(num_columns_list)):
    matrixs_i = files[i]
    matrixs_last_40[i]  = matrixs_i[-40:,-40:]

np.save('matrixs_last_40.npy',matrixs_last_40)
print('end')




# informations_np = np.asarray(informations)
# informations_np_pd = pd.DataFrame(informations_np[indexs,0,:])
# informations_np_pd.to_csv(os.path.join(directory_path, file_name[1]+'_informations.csv'))
# np.save(os.path.join(directory_path, file_name[1]+'_matrix.npy'),matrixs_114)

