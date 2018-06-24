# to suppress warnings and errors
import sys,os
sys.stdout = open(os.devnull, 'w')  # redirect the real STDOUT
sys.stderr = open(os.devnull, 'w')  # redirect the real STDOUT
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# import necessary modules
from keras.preprocessing.image import ImageDataGenerator
from keras.models import load_model
import numpy as np

CNN_path1 = str(sys.argv[1]) #"PCAICA_CNN.h5"
# print(CNN_path1)
Folder_path = str(sys.argv[2]) #"ROI_images"
# print(Folder_path)

# lets load in the CNN
model = load_model(CNN_path1)
# model.summary()
gen = ImageDataGenerator(featurewise_center=True, featurewise_std_normalization=True)
newvalid_path = Folder_path
im_size = [20, 20]
test_batch = gen.flow_from_directory(newvalid_path, target_size=(im_size[0], im_size[1]),
                                     batch_size=10000,
                                     color_mode='grayscale', shuffle=False)
classification_im = model.predict_generator(test_batch, 1)
# print classification_im
# save
savepathnewnpy = newvalid_path + '/Predictions'
np.save(savepathnewnpy, classification_im)