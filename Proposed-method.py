#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.image as mpimg
import pydicom
import time


PATH = "D:\\下載\\gLymph test-20200429T121458Z-001\\gLymph test\\S5010 T2\\"



def read_file(direc, file, normalization):
    file = pydicom.read_file(direc+file)
    img = file.pixel_array
    if normalization == 1:
        data = img/(np.max(img)-np.min(img))  #normalization
    else:
        data = img
    return img, data



def initial_curve(center, width, img_size_row, img_size_col):
    initial_map = np.zeros((img_size_row,img_size_col))
    for i in range(len(initial_map)):
        for j in range(len(initial_map[0])):
            initial_map[i][j] = np.sqrt((i-center[0])**2+(j-center[1])**2)-width
    return initial_map


def force(image,s_gradient, mean1, mean2):
    edge_constrain = 1/(1+s_gradient)
    f1 = (image-mean1)**2
    f2 = (image-mean2)**2
    eps = 1e-10
    return (1/2+f1*edge_constrain)/((1+(f1+f2)*edge_constrain+2*eps)), (1/2+f2*edge_constrain)/((1+(f1+f2)*edge_constrain+2*eps)) 


def delta(phi, eps):
    return 1/pi*(eps/(phi**2+eps**2)) #分析最大值


def gradient_square(img):
    return (np.gradient(img)[1])**2+(np.gradient(img)[0])**2


def K(phi):
    phi_x = np.gradient(phi)[1]
    phi_y = np.gradient(phi)[0]
    phi_xx = np.gradient(phi_x)[1]
    phi_xy = np.gradient(phi_x)[0]
    phi_yy = np.gradient(phi_y)[0]
    
    return (phi_xx * phi_y**2 + phi_yy * phi_x**2 - 2*phi_xy * phi_y*phi_x)/((phi_x**2 + phi_y**2)**(3/2)+1e-10) 




def narrow(phi):
    phi[abs(phi)<1e-3] = 0
    return phi



def show_contour(img, phi):
    plt.clf()
    plt.imshow(img,cmap = plt.cm.bone)
    plt.contour(phi,[0],colors='r',linewidths = 0.5) 
    plt.show()


def main(data, phi, gamma1, gamma2, eta, gs_img, c1, c2, eps):
    f1, f2 = force(data, gs_img , c1, c2)
    curvature = K(phi)
    phi_t = phi - delta(phi,eps)*(gamma1*f1-gamma2*f2-eta*curvature)
    c1_t = np.mean(data[phi_t>0])
    c2_t = np.mean(data[phi_t<0])
    energy.append(functional(gamma1, gamma2, eta, f1, f2, phi_t, eps))
    return c1_t, c2_t, phi_t



def functional(gamma1, gamma2, eta, f1, f2, phi, eps):
     return np.sum((gamma1*f1*((phi>0).astype('int')))+gamma2*(f2*((phi<0).astype('int')))+eta*delta(phi, eps)*gradient_square(phi)**0.5)
    



center = np.array([100,250])
width = 30

# img, data = read_file(PATH, 'IM100', 1)

# from skimage import io
# img = io.imread(PATH+'egle.png', as_gray=True)
# data = img/(np.max(img)-np.min(img))
# plt.imshow(img)
# plt.show()

img_size_row, img_size_col= len(img), len(img[0])
phi = narrow(initial_curve(center, width, img_size_row, img_size_col))
phi = phi/(np.max(phi)-np.min(phi))
gs_img = gradient_square(data)

c1 = np.mean(data[phi>0]) 
c2 = np.mean(data[phi<0]) 
start = time.time()
energy = []
for t in range(0,10001):
    c1_img = np.mean(img[phi>0]) 
    c2_img = np.mean(img[phi>0]) 
    c1_t, c2_t, phi_t = main(data, phi, 0.1, 0.1, 0.001, gs_img, c1, c2, .5)
    
    c1_img_t = np.mean(img[phi>0]) 
    c2_img_t = np.mean(img[phi>0])
    error = np.sum(np.abs(phi-phi_t)**2)**(1/2)
    if t % 30 == 0:
        show_contour(img, phi)
        print(error)
# #         print(abs(c2-c2_t))
        
    if t >0:
        if error<1:
            break
#     print(np.sum(np.abs(phi-phi_t)**2)**(1/2))
    phi = phi_t        
    c1 = c1_t
    c2 = c2_t
end = time.time()
print(t)
print(end-start)

plt.plot(energy)




from matplotlib.pyplot import figure


plt.clf()
plt.contour(phi,[0],colors='r',linewidths = 0.5)
plt.axis( 'image')
plt.gca().invert_yaxis()
plt.show()



c_1 = np.mean(img*((phi>=0).astype(int)))









c_2 = np.mean(img*((phi<=0).astype(int)))





restore = np.zeros((img_size_row, img_size_col))
restore[(phi>=0)] = -c_1
restore[(phi<0)] = -c_2
plt.clf()
plt.imshow(restore, cmap = plt.cm.bone)
# plt.imsave( 'proposed_without_constrain.png',restore, dpi = 500)
plt.show()



plt.clf()
plt.imshow(img, cmap = plt.cm.bone)
plt.show()



plt.clf()
plt.imshow(-1*img*((phi<0).astype(int)), cmap = plt.cm.bone)
# plt.imsave( 'proposed_without_constrain_WM.png',img*((phi<0).astype(int)),cmap = plt.cm.bone,  dpi = 500)
plt.show()





edge = ((np.gradient(restore)[0])**2+(np.gradient(restore)[1])**2)
# print(np.mean(edge[edge>0]), np.var(edge))
plt.clf()
plt.imshow(edge.astype(int), cmap = plt.cm.bone)
# plt.imsave('edge_map.png',edge*(edge>150).astype(int), cmap = plt.cm.bone)
plt.show()
