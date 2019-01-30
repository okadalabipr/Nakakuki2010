#%%
%run -i ../model/setParamConst.py
%run -i ../model/setVarEnum.py
%run -i ../model/diffeq.py
%run -i ../model/initialValues.py
%run -i ../runSim.py
%matplotlib inline
%run -i ../plot.py
plt.savefig('../Nakakukiet_al_2016.png',bbox_inches='tight')