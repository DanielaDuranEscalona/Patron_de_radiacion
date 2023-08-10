# -*- coding: utf-8 -*-



import numpy as np
import matplotlib.pyplot as plt


#%%% Creación de un patron de radiación con distribución de brillo gaussiana en el plato.

def RadiationPatern_gaussiana(radio,x,dx, lamda, mu, sigma):
    
    #Definición de variables locales
    y=x
    dy=dx
    x0=x/(2*dx)
    y0=x0
    Nx=int(x/dx)
    Ny=int(y/dy)
    
    #creacion de la matriz
    matrix=np.full((Nx,Ny),0,dtype=float)
    
    
    #creación del patron en la matriz considerando una distribución gaussiana centrada en el plato
    for i in range(104):
        for j in range(104):
            if ((i-x0/dx)**2 + (j-y0/dy)**2)**0.5 <= radio/dx:
                matrix[i,j]= np.exp( -(((i-x0)**2+(j-y0)**2)/(2 * sigma**2)))
            else:
                matrix[i,j]=0
                    
    iluminacion_dB=10*np.log10(matrix)
    

    #Patron de radiación, el patron de radiación corresponde a la transformada de fourier de la distribucion de brillo en la apertura
    f=np.fft.fft2(matrix)
    fs=abs(np.fft.fftshift(f)) #
    f_amp=fs/np.amax(fs)
    f_pwd=f_amp**2
    
    #Patron de radiacion Normalizado
    normalizado= 10*np.log10(f_pwd) 

    
    #Cambio de ejes del par de fourier 
    std=np.std(normalizado)
    FWHM=2.355*std
    
    phi_xmax= 1/(dx/lamda) #phi maximo
    phi_ymax= 1/(dy/lamda)
    delta_phix=phi_xmax/(Nx) #delta phi
#    delta_phiy=phi_ymax/(Ny-0.5)
    
    eje=np.arange(-phi_xmax/2,phi_xmax/2,delta_phix)*(180/np.pi)*60

    print(FWHM)
    
    plt.figure(1)
    plt.imshow(matrix,extent=[-x/2,x/2,-x/2,x/2])
    #plt.plot(fs)
    plt.colorbar()
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Distribucion de brillo en al apertura')
    plt.show()

    
    plt.figure(2)
    plt.imshow(iluminacion_dB,extent=[-x/2,x/2,-x/2,x/2])
    #plt.plot(fs)
    plt.colorbar()
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Distribucion de brillo en al apertura en dB')
    plt.show()

    
    plt.figure(3)
    plt.imshow(normalizado,extent=[-phi_xmax/2*(180/np.pi)*60,phi_xmax/2*(180/np.pi)*60,-phi_ymax/2*(180/np.pi)*60,phi_ymax/2*(180/np.pi)*60])
    #plt.plot(fs)
    plt.colorbar()
    plt.xlabel('$r \phi_x [arcmin]$')
    plt.ylabel('$r \phi_x [arcmin]$')
    plt.title('Patron de radiación')
    plt.show()
    
    
    plt.figure(4)
    #plt.imshow(matrix)
    plt.plot(eje,normalizado[int(Nx/2),:])
    plt.axhline(-FWHM)
    plt.axhline(-3)
    #plt.colorbar()
    plt.xlabel('$r \phi_x [arcmin]$')
    plt.ylabel('Amplitud [dB]')
    plt.title('Seccion transversal del patron de radiación')
    plt.show()

        
#%%% Creación de un patron de radiación con distribución de brillo uniforme al rededor del plato.

def RadiationPatern_uniforme(radio,x,dx, lamda):
    
    y=x
    dy=dx
    x0=x/2
    y0=x0
    Nx=int(x/dx)
    Ny=int(y/dy)
    matrix=np.full((Nx,Ny),0,dtype=float)
    
    for i in range(Nx):
        for j in range(Ny):
            radi=((i-x0/dx)**2 + (j-y0/dy)**2)**0.5
            if radi <= radio/dx:
                matrix[i,j]=1
            else:
                matrix[i,j]=0
                    
                    
    iluminacion_dB=np.log10(matrix)
    

    #Patron de radiación, el patron de radiación corresponde a la transformada de fourier de la distribucion de brillo en la apertura
    f=np.fft.fft2(matrix)
    fs=abs(np.fft.fftshift(f)) #Patron de radiacion asi na mas
    f_amp=fs/Nx
    f_pwd=f_amp**2
    
    #Patron de radiacion Normalizado
    normalizado= 10*np.log10(f_pwd/np.amax(f_pwd)) 
    
    
    #Cambio de ejes del par de fourier 
    
    phi_xmax= 1/(dx/lamda) #phi maximo
    phi_ymax= 1/(dy/lamda)
    delta_phix=phi_xmax/(Nx-0.5) #delta phi
#    delta_phiy=phi_ymax/(Ny-0.5)
    
    eje=np.arange(-phi_xmax/2,phi_xmax/2,delta_phix)*(180/np.pi)*60

    
    
    plt.figure(1)
    plt.imshow(matrix,extent=[-x/2,x/2,-x/2,x/2])
    #plt.plot(fs)
    plt.colorbar()
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Distribucion de brillo en al apertura')
    
    
    plt.figure(2)
    plt.imshow(iluminacion_dB,extent=[-x/2,x/2,-x/2,x/2])
    #plt.plot(fs)
    plt.colorbar()
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Distribucion de brillo en al apertura en dB')
    
    
    plt.figure(3)
    plt.imshow(normalizado,extent=[-phi_xmax/2*(180/np.pi)*60,phi_xmax/2*(180/np.pi)*60,-phi_ymax/2*(180/np.pi)*60,phi_ymax/2*(180/np.pi)*60])
    #plt.plot(fs)
    plt.colorbar()
    plt.xlabel('$r \phi_x [arcmin]$')
    plt.ylabel('$r \phi_x [arcmin]$')
    plt.title('Patron de radiación')
    
    
    plt.figure(4)
    #plt.imshow(matrix)
    plt.plot(eje,normalizado[int(Nx/2),:])
    plt.axhline(-3)
    #plt.colorbar()
    plt.xlabel('$r \phi_x [arcmin]$')
    plt.ylabel('Amplitud [dB]')
    plt.title('Seccion transversal del patron de radiación')
    plt.show()


#%%% Main code
if __name__=='__main__':

    #valores fijos del telescopio y de la grilla
    x=21    #ancho de la grilla en metros un poco mas amplio que el diametro
    dx=0.01 #espaciamiento dx
    radio=0.52#Radio del telescopio en metros
    
    
    mu=0  #Fijo para poder tener el centro de la gausiana en el centro del telescopio 
    sigma=70000
    lamda=0.001 #longitud de onda en metros (1mm=0.001m)
    
    sigma1=2300
    lamda1=0.003 #longitud de onda en metros (3mm=0.003)
    
    
    sigma2=100
    lamda2=0.000035
    
    #Patron de radiación para \lambda = 1mm
    
#    RadiationPatern_gaussiana(radio,x,dx, lamda, mu, sigma)
#    RadiationPatern_gaussiana(radio,x,dx, lamda1, mu, sigma1)
#    RadiationPatern_gaussiana(radio,x,dx, lamda2, mu, sigma2)
    
    RadiationPatern_uniforme(radio,x,dx, lamda)
    
    
