import numpy as np
import pydicom
import time
import matplotlib.pyplot as plt

class segmentation():
    def __init__(self, image, contour_center, contour_radius=30, gamma1=1, gamma2 = 1, eta = 1e-3, iteration=200, eps=.5, stab_eps = 1e-10):
        self.iteration = iteration
        self.delta_eps = eps
        self.numerical_eps = stab_eps
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.eta = eta

        self.image = image
        self.data = image/(np.max(image)-np.min(image))
        
        self.row = len(image)
        self.column = len(image[0])

        self.contour_center = contour_center
        self.contour_radius = contour_radius
        self.init_contour()
        self.gradients()
        self.update_gradient()
        self.data_gradient()
        
    def init_contour(self):
        center = self.contour_center
        radius = self.contour_radius
        self.contour = np.array([[round(((i-center[1])**2+(j-center[0])**2)**0.5-radius,2) for i in range(self.row)] for j in range(self.column)])
        
    def data_gradient(self):
		self.data_gradient_magnitude = self.grad_x**2 + self.grad_y**2
    
    def kappa(self):
        # phi = self.contour
        phi_x = np.gradient(self.contour)[1]
        phi_y = np.gradient(self.contour)[0]
        phi_xx = np.gradient(phi_x)[1]
        phi_xy = np.gradient(phi_x)[0]
        phi_yy = np.gradient(phi_y)[0]
        
        self.curvature = (phi_xx * phi_y**2 + phi_yy * phi_x**2 - 2*phi_xy * phi_y*phi_x + self.numerical_eps)/((phi_x**2 + phi_y**2)**(3/2)+self.numerical_eps) 
    
    
    def gradients(self):
        self.grad_x = np.gradient(self.data)[1]
        self.grad_y = np.gradient(self.data)[0]
    
    
    def force(self):
        
        edge_constrain = 1/(1+self.data_gradient_magnitude)
        # c1, c2 = self.means()
        f1 = (self.data-self.c1)**2 #>0
        f2 = (self.data-self.c2)**2 #<0
        
        self.force1 = (1/2+f1*edge_constrain)/((1+(f1+f2)*edge_constrain+2*stab_eps))
        self.force2 = (1/2+f2*edge_constrain)/((1+(f1+f2)*edge_constrain+2*stab_eps))
    
    def delta_function(self):
        self.delta = 1/pi*(eps/(self.contour**2+eps**2)) #分析最大值
    
    def minimizer(self):
        # curvature = self.curvature
        self.contour = self.contour[np.abs(self.contour)<=1e-3] = 0
        self.contour = self.contour/(np.max(self.contour)-np.min(self.contour))
        self.c1 = np.mean(self.data[self.contour>0])
        self.c2 = np.mean(self.data[self.contour<0])

        for t in range(self.iteration):
            phi_t = self.contour - self.delta * (self.gamma1*self.force1 - self.gamma2*self.force2 - self.eta*self.curvature)
            
