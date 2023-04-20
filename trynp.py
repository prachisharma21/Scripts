import numpy as np 
import matplotlib.pyplot as plt

def band_structure_densityplot(datafile = None, smearing= 0.5,Efermi=None,Energyboundary = 5,dE = 200, savefig_name = "density_plot"):
    """ A function to plot the density plot for band structure data 
        Parameters: 
        datafile = DFT output file from quantum espresso 
        smearing = smaering parameter for the guassian to give density 
        Efermi = Fermi energy
        Energyboundary = The energy width of the plot 
        dE = energy step for the points.   
    """
    Emin = -Energyboundary + Efermi 
    Emax =  Energyboundary + Efermi
    energy=np.linspace(Emin,Emax,dE) 

    data = datafile

    # select only the data within the energy bounds 
    data = data[(data[:,1]>=Emin-max(smearing*10,0.1))*(data[:,1]<=Emax+max(smearing*10,0.1))]
    
    # unique k-points in the data
    nk = np.unique(data[:, 0])
    density=np.zeros((len(nk),dE),dtype=float)
    for k,E,w in data[:,:3]:
        # To find the k-point along the path 
        # index of the minimum in the list 
        ik=np.argmin(abs(k-nk))
        density[ik,:]+=w*np.exp( -(energy-E)**2/(2*smearing**2))
    
    
    k1=np.unique(data[:, 0])/0.5093340253675882
    E1=energy-Efermi
    k1,E1=np.meshgrid(k1,E1)
    plt.pcolormesh(k1,E1,density.T, cmap='magma', shading='auto')
    #plt.colorbar()
    plt.ylabel('$\epsilon - \epsilon_F$ (eV)',fontsize=20)
    plt.axvline(1.2581, linewidth=0.5, color='y', alpha=1)
    plt.axvline(2.5162, linewidth=0.5, color='y', alpha=1)
    plt.axvline(3.8095, linewidth=0.5, color='y', alpha=1)
    plt.axvline(4.3284, linewidth=0.5, color='y', alpha=1)
    plt.axhline(0, linewidth=0.5, color='y', alpha=0.8)
    plt.xticks(ticks= [0, 1.2581, 2.5162, 3.8095, 4.3284], 
           labels=['$\Gamma$','L', 'F','$\Gamma$','T'],fontsize=20)
    plt.title(r'$PdCo(Fe)O_2$ relaxed structure, Fe defect')
    plt.savefig(savefig_name+".png")

filename = input("Enter the file name to be plotted: ")
result_no_defect = np.loadtxt(filename)

band_structure_densityplot(datafile=result_no_defect,smearing=0.05,Efermi=14.532,Energyboundary=5,dE=500,savefig_name="no_defect")