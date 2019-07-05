from simulation import Simulation
from viz import plot_func

def run_simulation():
    sim = Simulation()
    sim.numerical_integration()

    plot_func(sim)


if __name__ == "__main__":
    run_simulation()