from model.param_const import f_params
from model.initial_condition import initial_values
from simulation import Simulation
from viz import plot_func

def run_simulation():
    x = f_params()
    y0 = initial_values()

    sim = Simulation(x,y0)
    sim.numerical_integration(x,y0)

    plot_func(sim)


if __name__ == "__main__":
    run_simulation()