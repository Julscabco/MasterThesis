# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:35:59 2023

@author: Usuario
"""

import napari
from skimage.data import cells3d


# create a `Viewer` and `Image` layer here
viewer, image_layer = napari.imshow(cells3d())

# print shape of image datas
print(image_layer.data.shape)

# start the event loop and show the viewer
napari.run()
