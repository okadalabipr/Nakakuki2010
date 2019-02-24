#%%
%run -i model/f_parameter.py
%run -i model/f_variable.py
%run -i model/differential_equation.py
%run -i model/initial_condition.py
%run -i simulation.py
%run -i plot_func.py
plt.savefig('./cfosmodel.png',bbox_inches='tight')