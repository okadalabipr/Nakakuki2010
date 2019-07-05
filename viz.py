import numpy as np
import matplotlib.pyplot as plt

def plot_func(sim):
    plt.figure(figsize=(20,8))
    plt.rcParams['font.size'] = 16
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['lines.linewidth'] = 2.5
    plt.rcParams['lines.markersize'] = 10

    plt.subplots_adjust(wspace=0.5, hspace=0.5)

    plt.subplot(2,4,1)
    plt.plot(sim.t,sim.PMEK_cyt[:,0],'b',label='EGF')
    plt.plot(sim.t,sim.PMEK_cyt[:,1],'r',label='HRG')
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated MEK\n(cytoplasm)')
    plt.legend(loc='upper right',fontsize=12,frameon=False)

    plt.subplot(2,4,2)
    plt.plot(sim.t,sim.PERK_cyt[:,0]/np.max(sim.PERK_cyt[:,1]),'b')
    plt.plot(sim.t,sim.PERK_cyt[:,1]/np.max(sim.PERK_cyt[:,1]),'r')
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated ERK\n(cytoplasm)')

    plt.subplot(2,4,3)
    plt.plot(sim.t,sim.PRSK_wcl[:,0]/np.max(sim.PRSK_wcl[:,1]),'b')
    plt.plot(sim.t,sim.PRSK_wcl[:,1]/np.max(sim.PRSK_wcl[:,1]),'r')
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated RSK\n(whole cell)')

    plt.subplot(2,4,4)
    plt.plot(sim.t,sim.PCREB_wcl[:,0]/np.max(sim.PCREB_wcl[:,1]),'b')
    plt.plot(sim.t,sim.PCREB_wcl[:,1]/np.max(sim.PCREB_wcl[:,1]),'r')
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated CREB\n(whole cell)')

    plt.subplot(2,4,5)
    plt.plot(sim.t,sim.DUSPmRNA[:,0]/np.max(sim.DUSPmRNA[:,1]),'b')
    plt.plot(sim.t,sim.DUSPmRNA[:,1]/np.max(sim.DUSPmRNA[:,1]),'r')
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel(r'$\it{dusp}$'+' mRNA\nexpression')

    plt.subplot(2,4,6)
    plt.plot(sim.t,sim.cFosmRNA[:,0]/np.max(sim.cFosmRNA[:,1]),'b')
    plt.plot(sim.t,sim.cFosmRNA[:,1]/np.max(sim.cFosmRNA[:,1]),'r')
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel(r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA\nexpression')

    plt.subplot(2,4,7)
    plt.plot(sim.t,sim.cFosPro[:,0]/np.max(sim.cFosPro[:,1]),'b')
    plt.plot(sim.t,sim.cFosPro[:,1]/np.max(sim.cFosPro[:,1]),'r')
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('c-Fos Protein\nexpression')

    plt.subplot(2,4,8)
    plt.plot(sim.t,sim.PcFos[:,0]/np.max(sim.PcFos[:,1]),'b')
    plt.plot(sim.t,sim.PcFos[:,1]/np.max(sim.PcFos[:,1]),'r')
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated c-Fos\nProtein expression')

    plt.show()