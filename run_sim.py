def include(pyfile):
    with open(pyfile,'r',encoding='utf-8') as f:
        script = f.read()
        exec(script,globals())


include('model/f_parameter.py')
include('model/f_variable.py')
include('model/differential_equation.py')
include('model/initial_condition.py')

include('simulation.py')
include('plot_func.py')

# plt.savefig('./cfosmodel.png',bboc_inches='tight')
plt.show()