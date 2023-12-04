# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 10:10:57 2023

@author: Usuario
"""
import os

import rawpy
import imageio
import napari
import skimage 
from skimage import io
import numpy as np
import matplotlib.pyplot as plt

from skimage import measure, color
from skimage.measure import label
from skimage.segmentation import clear_border
from skimage.filters import threshold_otsu
from skimage import exposure
from skimage.morphology import binary_closing, binary_erosion, binary_dilation

import cv2

def euclidean_distance(cx,cy,x,y):
    return np.sqrt(((x-cx)**2.0)+(y-cy)**2.0)


# Image reading
image_name = 'CRW_6417.crw'
image_folder = os.path.join(os.getcwd(),'Platephoto2Nov2023')
image_path = os.path.join(image_folder,image_name)

with rawpy.imread(image_path) as raw:
    img = raw.postprocess()

io.imshow(img)
io.show()

# We see what channel is the best-looking channel
r = img[:,:,0]
g = img[:,:,1]
b = img[:,:,2]

fig, ax = plt.subplots(3,figsize=(15,6))
ax[0].set_ylabel('Red channel')
ax[0].imshow(r)
ax[1].set_ylabel('Green channel')
ax[1].imshow(g)
ax[2].set_ylabel('Blue channel')
ax[2].imshow(b)

# We generates the mask that selects only the content inside the petri dish
thresh = 110
mask_r = r > thresh

m,n = np.shape(r)

first_row = np.zeros(n)
last_row = np.zeros(n)

first_row = mask_r[0,:]
last_row = mask_r[-1,:]

first_row_indexes = np.where(first_row==True)[0]
last_row_indexes = np.where(last_row==True)[0]

cutmin = np.max([first_row_indexes[0],last_row_indexes[0]])
cutmax = np.min([first_row_indexes[-1],last_row_indexes[-1]])

r = r[:,cutmin:cutmax]


io.imshow_collection([r,mask_r],cmap='gray')
io.show()


# Viewing of the petri dish selected
petri_dish = r

# Image pre processing
enhanced_image = exposure.equalize_adapthist(petri_dish)
print(np.max(enhanced_image))
print(np.min(enhanced_image))
thresh = (enhanced_image<0.65)

io.imshow_collection([r,thresh])
io.show()

erosion = binary_erosion(thresh,footprint= np.ones((15,15)))

io.imshow_collection([thresh,erosion])
io.show()

dilated = binary_dilation(erosion,footprint=np.ones((20,20)))

io.imshow_collection([thresh,dilated])

label_image = label(dilated)

plt.imshow(color.label2rgb(label_image, bg_label=0))
plt.show()

binary_image = (label_image * 255).astype(np.uint8)

# Apply Hough Circle Transform on the labeled image
detected_circles = cv2.HoughCircles(
    r, cv2.HOUGH_GRADIENT, dp=1, minDist=20, param1=50, param2=30, minRadius=5, maxRadius=20
)

label_image_with_circles = r.copy()
print(detected_circles)
# Draw the detected circles
if detected_circles is not None:
    detected_circles = np.uint16(np.around(detected_circles))
    for i in detected_circles[0, :]:
        center = (i[0], i[1])
        radius = i[2]
        cv2.circle(label_image_with_circles, center, radius, (0, 255, 0), 2)

# Display the original labeled image and the image with detected circles
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.imshow(r)
plt.title("Original Labeled Image")

plt.subplot(1, 2, 2)
plt.imshow(label_image_with_circles)
plt.title("Labeled Image with Detected Circles")
plt.show()


