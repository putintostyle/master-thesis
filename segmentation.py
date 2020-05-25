import numpy as np
import pydicom
import time
import matplotlib.pyplot as plt

class segmentation():
    def __init__(self, image, contour_center, contour_radius, eps=.5, stab_eps = 1e-10):
        self.image = image
        self.data = image/(np.max(image)-np.min(image))
        
        self.row = len(image)
        self.column = len(image[0])

        self.contour_center = contour_center
        self.contour_radius = contour_radius
        self.init_contour()
        self.update_gradient()

        
    def init_contour(self):
        center = self.contour_center
        radius = self.contour_radius
        self.contour = np.array([[round(((i-center[1])**2+(j-center[0])**2)**0.5-radius,2) for i in range(self.row)] for j in range(self.column)])
        
    def data_gradient(self):
        self.grad_x, self.grad_y = self.gradients(self.data)
		self.data_gradient_magnitude = np.power(self.grad_x, 2) + np.power(self.grad_y, 2)
    
    def kappa(self):
        phi = self.contour
        phi_x = np.gradient(phi)[1]
        phi_y = np.gradient(phi)[0]
        phi_xx = np.gradient(phi_x)[1]
        phi_xy = np.gradient(phi_x)[0]
        phi_yy = np.gradient(phi_y)[0]
        
        self.curvature = (phi_xx * phi_y**2 + phi_yy * phi_x**2 - 2*phi_xy * phi_y*phi_x + stab_eps)/((phi_x**2 + phi_y**2)**(3/2)+stab_eps) 
    
    def narrow(self):
        self.contour = self.contour[np.abs(self.contour)<=1e-3] = 0

    def gradients(self):
        return np.gradient(self.data)[1], np.gradient(self.data)[0]
    
    def mean(self):
        return np.mean(self.data[self.contour>0]), np.mean(self.data[self.contour<0])
    def force(self):
        
        edge_constrain = 1/(1+self.data_gradient_magnitude)
        c1, c2 = self.mean()
        f1 = (self.data-c1) #>0
        f2 = (self.data-c2) #<0
        =
        self.force1 = (1/2+f1*edge_constrain)/((1+(f1+f2)*edge_constrain+2*stab_eps))
        self.force2 = (1/2+f2*edge_constrain)/((1+(f1+f2)*edge_constrain+2*stab_eps))
    
    def delta_function(self):
        self.delta = 1/pi*(eps/(self.contour**2+eps**2)) #分析最大值
    def
        
