


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numpy.lib.scimath import sqrt as csqrt
from scipy.special import fresnel


def TMM(period, section_length, ref_mod, delta, alpha):
    
    T11 = np.cosh(alpha*section_length)-1j*(delta/alpha)*np.sinh(alpha*section_length)
    T22=np.cosh(alpha*section_length)+1j*(delta/alpha)*np.sinh(alpha*section_length)
    T12=-1j*(kapaac/alpha)*np.sinh(alpha*section_length)
    T21=1j*(kapaac/alpha)*np.sinh(alpha*section_length)

    return np.array([[T11,T12],[T21,T22]])





def Fres_Diff(write_wavelength, slitwidth, x_pos, screen_pos):
    alpha_1 = (x_pos + slitwidth/2)*(np.sqrt(2/(write_wavelength*screen_pos)))
    alpha_2 = (x_pos - slitwidth/2)*(np.sqrt(2/(write_wavelength*screen_pos)))
    s1, c1 = fresnel(alpha_1) #Fresnel integrals
    s2, c2 = fresnel(alpha_2) #Fresnel integrals
    Intensity = 0.5*((c2-c1)**2+(s2-s1)**2)
    return Intensity

# constants - all dimensional untis in m

r_core = 5e-6 # radius of the core
lambda_cut = 700e-9 # wavelength in mm
n_core = 1.445 # refractive idnex of the core
n_clad = np.sqrt(-(2.405**2*lambda_cut**2)/(4*np.pi**2*r_core**2)+n_core**2) # claculate cladding oindex given cut off wavelength and core radius and refractive index
d_n_max = 1e-6 # magnitude of index modulation


Bragg_wavelength = 750e-9
period = Bragg_wavelength/(2*n_core)
L = 2e-3
no_sections = 100
section_length = L/no_sections

lambda_mid = 750e-9 
lambda_resolution = 10e-12
lambda_range = 10e-9
lambda_points = np.int(lambda_range/lambda_resolution)


wavelength = np.linspace(lambda_mid-lambda_range/2,lambda_mid+lambda_range/2,lambda_points)

ko = 2*np.pi/wavelength

#calculating effective index of fibre mode

vf=((2*np.pi/wavelength)*r_core*np.sqrt(n_core**2-n_clad**2))

u = (1+np.sqrt(2))*vf/(1+(4+vf**4)**0.25)
bfib = 1-u**2/vf**2

n_eff = np.sqrt(bfib*(n_core**2-n_clad**2)+n_clad**2)  
R=np.zeros((lambda_points))
dx = np.int(section_length/period)*period
Length = dx*no_sections
section_length = L/no_sections

## flat top
#delta_n  = d_n_max*np.ones(no_sections)
#x = np.linspace(-Length/2,Length/2,no_sections)



## Guassian
#x = np.linspace(-Length/2,Length/2,no_sections)
#delta_n = d_n_max*np.exp(-(2*x/(0.4*Length))**2)

## Raised Cosine
#x = np.linspace(-Length/2,Length/2,no_sections)
#delta_n = 0.5*d_n_max*(1+np.cos(np.pi*(x/(0.1*Length))))


# Fresnel single slit-width and screen posi
screen_pos = 0.45
write_wavelength = 266e-9
slitwidth = 0.5e-3
x = np.linspace(-Length/2,Length/2,no_sections)
intensity = Fres_Diff(write_wavelength, slitwidth, x, screen_pos)
delta_n = d_n_max*intensity/np.max(intensity)



for i in range(0,lambda_points):
    MAT = np.eye(2)
    


    for m in range(no_sections-1,0,-1):
        kapadc = 4*np.pi*delta_n[m]/wavelength[i]
        kapaac = kapadc/2
        delta = kapadc+0.5*(2*ko[i]*n_eff[i]-2*np.pi/period)
        alpha = csqrt(kapaac**2-delta**2)

        T = TMM(period, section_length, delta_n[m], delta, alpha)
        
        MAT =np.dot(MAT,T)
        
    R[i] = np.abs(MAT[0,1]/MAT[0,0]) 

   
   

#Fresnel variable slit-width fixed screen distance
screen_pos = 350e-3
write_wavelength = 266e-9
x = np.linspace(-Length/2,Length/2,no_sections)
slitwidth = np.linspace(0.1e-3,1e-3,100) # slitwidth 

Spectrum1 = np.zeros((len(wavelength),len(slitwidth))) # set up array to receive the intensity data


for i in range(1,len(slitwidth)-1):
    intensity = Fres_Diff(write_wavelength, slitwidth[i], x, screen_pos)
    delta_n = d_n_max*intensity/np.max(intensity)
    for p in range (0,lambda_points):
        MAT = np.eye(2)
        for m in range(no_sections-1,0,-1):
            kapadc = 4*np.pi*delta_n[m]/wavelength[p]
            kapaac = kapadc/2
            delta = kapadc+0.5*(2*ko[p]*n_eff[p]-2*np.pi/period)
            alpha = csqrt(kapaac**2-delta**2)

            T = TMM(period, section_length, delta_n[m], delta, alpha)
            
            MAT =np.dot(MAT,T)
        
        Spectrum1[p,i] = np.abs(MAT[0,1]/MAT[0,0]) 






#Fresnel fixed slit-width variable screen distance

screen_pos = np.linspace(0,1,1000)
write_wavelength = 266e-9
x = np.linspace(-Length/2,Length/2,no_sections)
slitwidth = 0.75e-3 

Spectrum2 = np.zeros((len(wavelength),len(screen_pos))) # set up array to receive the intensity data


for i in range(1,len(screen_pos)-1):
    intensity = Fres_Diff(write_wavelength, slitwidth, x, screen_pos[i])
    delta_n = d_n_max*intensity/np.max(intensity)
    for p in range (0,lambda_points):
        MAT = np.eye(2)
        for m in range(no_sections-1,0,-1):
            kapadc = 4*np.pi*delta_n[m]/wavelength[p]
            kapaac = kapadc/2
            delta = kapadc+0.5*(2*ko[p]*n_eff[p]-2*np.pi/period)
            alpha = csqrt(kapaac**2-delta**2)

            T = TMM(period, section_length, delta_n[m], delta, alpha)
            
            MAT =np.dot(MAT,T)
        
        Spectrum2[p,i] = np.abs(MAT[0,1]/MAT[0,0]) 



plt.close('all')
plt.figure(2)
plt.plot(wavelength,R)

plt.figure(3)
plt.plot(x,delta_n)


#plot range
d_min = np.min(slitwidth)
d_max = np.max(slitwidth)
wave_max = np.max(wavelength)
wave_min= np.min (wavelength)

plt.figure(4)
plt.imshow(Spectrum1, cmap=matplotlib.cm.gray, interpolation='nearest', aspect='auto', origin='lower' )
plt.xlabel('Slitwidth (mm)')
plt.ylabel('Wavelength (nm)')
plt.title("Bragg spectrum as function of slitwidth" + "\n" + "assuming writing wavelength 266 nm and viewing distance " +  str(screen_pos) + "mm")


plt.figure(5)
plt.imshow(Spectrum2, cmap=matplotlib.cm.gray, interpolation='nearest', aspect='auto', origin='lower' )
plt.xlabel('Slit to fibre separation(m)')
plt.ylabel('Wavelength (nm)')
plt.title("Bragg spectrum as function of ditance from slit to fibre" + "\n" + "assuming writing wavelength 266 nm and slit width " +  str(slitwidth*1e3) + "mm")
plt.show()


