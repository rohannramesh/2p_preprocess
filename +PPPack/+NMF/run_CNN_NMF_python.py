# to suppress warnings and errors
import sys,os
sys.stdout = open(os.devnull, 'w')  # redirect the real STDOUT
sys.stderr = open(os.devnull, 'w')  # redirect the real STDOUT
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# import necessary modules
from keras.preprocessing.image import ImageDataGenerator
from keras.models import load_model
import numpy as np
import scipy.io as sio

CNN_path1 = str(sys.argv[1]) #"PCAICA_CNN.h5"
# CNN_path1 = "/Users/rohanramesh/Documents/Harvard/Andermann/Code/NeuralNetworks/PCAICA_im_CNN_ROI_classify_train_OA27_OA26_OA67_45runs_wtc_NMF.h5"
# print(CNN_path1)
Folder_path = str(sys.argv[2]) #"ROI_images"
# Folder_path = "/Users/rohanramesh/Documents/Harvard/Andermann/Code/NeuralNetworks/Fitting/ROI_images"
# print(Folder_path)

# load in the tc distributions

tmp_load_var = sio.loadmat(str(sys.argv[3]))
# tmp_load_var = sio.loadmat("/Users/rohanramesh/Documents/Harvard/Andermann/Code/NeuralNetworks/Fitting/ROI_images/tc_hist.mat")
All_histograms = {}
All_histograms['ROIs_test'] = tmp_load_var['all_hists']


# to build a generator with the tc
def generatorWithNonFolder(generator1, other_input):
    while True:
            X1i = generator1.next()
            curr_idx = ((generator1.batch_index-1) * generator1.batch_size) % generator1.n
            idx = generator1.index_array[curr_idx:curr_idx + generator1.batch_size]
            X2i = other_input[idx,:]
            nkeep = X1i[0].shape[0]
            yield [X1i[0], X2i[0:nkeep,:]], X1i[1]



# lets load in the CNN
model = load_model(CNN_path1)
# model.summary()
gen = ImageDataGenerator(featurewise_center=True, featurewise_std_normalization=True)
newvalid_path = Folder_path
im_size = [20, 20]

test_batch = gen.flow_from_directory(newvalid_path, target_size=(im_size[0], im_size[1]),
             batch_size = np.shape(All_histograms['ROIs_test'])[0],
             color_mode='grayscale', shuffle=False, class_mode='categorical')

# combines the folder and the tc
ComboHistNewTest = generatorWithNonFolder(test_batch, All_histograms['ROIs_test'])
classification_im = model.predict_generator(ComboHistNewTest, 1)
# print(classification_im)
# save
savepathnewnpy = newvalid_path + '/Predictions'
np.save(savepathnewnpy, classification_im)