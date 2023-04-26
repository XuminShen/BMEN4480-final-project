import os
import  pandas as pd
import numpy as np
import matplotlib.pyplot as plt
directory_path = "C:/Users/xumin/Desktop/TCR/"
# file_name = ["Michael_align_data.csv","Michaelmatrix"]
file_name = ["1600_align_data.csv","1600"]

full_path_1 = os.path.join(directory_path, file_name[0])
full_path_2 = os.path.join(directory_path, file_name[1])

data = pd.read_csv(full_path_1)
matrixs = []
informations = []
for root, dirs, files in os.walk(full_path_2):
    for file in files:
        # path
        file_path = os.path.join(root, file)

        # return file name
        split_file_name = file.split('_')
        try:
            data1 = data[data['cdr3_chains']==split_file_name[0]]
            data2 = data1[data1['v_gene'] == split_file_name[2]]
            data3 = data2[data2['cdr3'] == split_file_name[1]]
            data4 = data3[data3['v_gene_vdj'] == split_file_name[3]]

            if data4.shape[0]==1:
                matrixs.append(pd.read_csv(file_path))
                informations.append(data4)
        except:
            continue

print(len(matrixs))
print(len(informations))

shapes = []
for i in np.arange(len(matrixs)):
    shapes.append(np.shape(matrixs[i])[0])

plt.hist(shapes,bins = np.max(shapes)-np.min(shapes))
plt.show()

indexs = np.arange(len(shapes))[np.asarray(shapes)<=114]
matrixs_114 = np.zeros([len(indexs),114,114])
for i in np.arange(len(indexs)):
    matrixs_i = np.asarray(matrixs[indexs[i]])
    matrixs_114[i,np.int32(57-shapes[indexs[i]]/2):np.int32(57+shapes[indexs[i]]/2),np.int32(57-shapes[indexs[i]]/2):np.int32(57+shapes[indexs[i]]/2)]  = matrixs_i[:,1:]

informations_np = np.asarray(informations)
informations_np_pd = pd.DataFrame(informations_np[indexs,0,:])
informations_np_pd.to_csv(os.path.join(directory_path, file_name[1]+'_informations.csv'))
np.save(os.path.join(directory_path, file_name[1]+'_matrix.npy'),matrixs_114)
