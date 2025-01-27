import rebound
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass



@dataclass(frozen=True)
class Constants:
    E = 0.
    GM = 1.0
if __name__ == "__main__": CONST = Constants()


class Simulation:
    _round = 1
    _dt = 0.1
    def __init__(self):
        pass

    def simulate(self, particle : 'P', integrator, round : float = _round, dt : float = _dt, **kwargs):
        """Simulate a particle for multiple steps around a central star with a given integrator

        Args:
            particle (P): particle to simulate
            integrator (_type_): integrator can be explicit euler, rungeKutta, ...
            round (float, optional): approx. number of rounds around the center star. Defaults to _round.
            dt (float, optional): steps-size. lower gives a better precision. Defaults to _dt.
        """
        many = int(np.ceil(round/dt * 2*np.pi / (1-particle.e)**(3/2)))
        ps = []
        ps.append(particle)
        for i in range(many):
            ps.append(integrator(ps[-1], dt))

        if 'mylabel' in kwargs:
            kwargs['label'] = f"{kwargs['mylabel']}_{f'{dt:0.2f}'[1:]}" + ("" if particle.e==0 else f"_e={particle.e:.2f}")
            del kwargs['mylabel']

        self.Plots(ps, **kwargs)

    def addPlot(self, x, y, **kwargs):
        qs = list(zip(x,y)) # dunno why list() is necessary
        self.Plots(qs, **kwargs)

    def addPlotApocenter(self, e : float, **kwargs):
        qs = [[1,0], [-(1+e)/(1-e), 0]]
        self.Plots(qs, **kwargs)

    def nextPlot(self, name : str):
        self.Plots._subplotsTitle.append(name)

    def plot(self, showPerfect = True):
        def getSubplot(i : int):
            row = i // ncols
            col = i % ncols
            return axs[row,col] if nrows>1 else axs[col]

        many = self.Plots._plots[0]._subs
        nrows = 1 if many<=3 else 2
        ncols = int(np.ceil(many/nrows))
        fig, axs = plt.subplots(nrows, ncols, figsize=(5*ncols,5*nrows))
        lnspace = np.array(np.linspace(0,2*np.pi))

        for plot in self.Plots._plots:
            subplot = getSubplot(plot.subplot)
            subplot.plot(plot.x, plot.y, **plot.kwargs)

        for i, title in enumerate(self.Plots._subplotsTitle):
            subplot = getSubplot(i)
            if showPerfect: subplot.plot(np.sin(lnspace), np.cos(lnspace), 'k', label='perfect')
            subplot.set_title(title)
            subplot.plot(0,0,'ko')
            subplot.legend()
            subplot.axis('equal')


        plt.show()


    class Plots:
        _plots = []
        _subplotsTitle = []
        def __init__(self, ps, **kwargs):
            self.subplot = self._subs-1
            self.x = [p[0] for p in ps]
            self.y = [p[1] for p in ps]
            self.kwargs = kwargs

            self._plots.append(self)

        @property
        def _subs(self):
            return len(self._subplotsTitle)



class Particle(np.ndarray):
    def __new__(cls, x: float, y: float, v: float = 0, w: float = 0):
        data = np.array([x, y, v, w], dtype=float).view(cls)
        return data
    def r(self) -> 'Particle':
        """returns the first 2 elements (r)"""
        return Particle(self[0], self[1])
    def v(self) -> 'Particle':
        """return the last 2 elements (v)"""
        return Particle(self[2], self[3])
    

    @property
    def f(self) -> 'Particle':
        def gmFunc(q : float) -> float:
            return -CONST.GM * q / np.sqrt(self[0]**2 + self[1]**2)**3             ## **3 needed?
        return Particle(self[2], self[3], gmFunc(self[0]), gmFunc(self[1]))
    
    @property
    def a(self) -> 'Particle':
        def gmFunc(q : float) -> float:
            return -CONST.GM * q / (np.sqrt(self[0]**2 + self[1]**2))**3        ## **3 wrong?
        return Particle(gmFunc(self[0]), gmFunc(self[1]), 0,0)

    def explicitEuler(self, dt : float) -> 'Particle':
        return self + self.f * dt
    
    def rungeKutta2(self, dt : float) -> 'Particle':
        k1 = self.f
        k2 = self.explicitEuler(dt).f
        return self + (k1+k2)/2 * dt
    
    def rungeKutta4(self, dt : float) -> 'Particle':
        k1 = self.f
        k2 = (self+k1*dt/2).f
        k3 = (self+k2*dt/2).f
        k4 = (self+k3*dt).f
        return self + (k1/6 + k2/3 + k3/3 + k4/6)*dt
    
    def leapFrog(self, dt : float) -> 'Particle':

        vnhalf = self.v() + self.a*dt/2
        rnplus1 = self.r() + vnhalf*dt
        vnplus1 = vnhalf + rnplus1.a*dt/2

        return Particle(rnplus1[0], rnplus1[1], vnplus1[0], vnplus1[1])
    
    def semiImplicitEuler(self, dt : float) -> 'Particle':
        vnplus1 = self.v() + self.r().a*dt
        rnplus1 = self.r() + vnplus1*dt # unclear if this is correct (first thought of vnplus1.a*dt, but dependent on Ansatz for f(v_n+1))

        return Particle(rnplus1[0], rnplus1[1], vnplus1[0], vnplus1[1])

class P(Particle):
    def __new__(cls, e : float = CONST.E):
        data = np.array([1,0,0,np.sqrt(1+e)], dtype=float).view(cls)
        data.e = e
        return data


if __name__ == "__main__":
    mpl.rcParams['font.family'] = 'monospace'
    mpl.rcParams['font.monospace'] = 'Courier New'

    assert 0 <= CONST.E < 1
    s = Simulation()

    """ USAGE
    - plot different particles with different integrators
    - **kwargs: label, color, marker, linestyle, linewidth, markersize, alpha
    """


    s.nextPlot(name="Compare EE, RK2, RK4")
    s.simulate(P(0), Particle.explicitEuler,   1, 0.25, mylabel=' ee',  color='green', alpha=.2)
    s.simulate(P(0), Particle.explicitEuler, 2.4, 0.05, mylabel=' ee',  color='green')
    s.simulate(P(0), Particle.rungeKutta2,     1, 0.50, mylabel='rk2', color='firebrick', alpha=.3, marker='+')
    s.simulate(P(0), Particle.rungeKutta2,     1, 0.05, mylabel='rk2', color='firebrick', marker='+')
    s.simulate(P(0), Particle.rungeKutta4,     1, 0.50, mylabel='rk4', color='blue', alpha=.3)
    s.simulate(P(0), Particle.rungeKutta4,     1, 0.05, mylabel='rk4', color='blue')


    s.nextPlot(name="Runge Kutta after 100 laps")
    s.simulate(P(0), Particle.rungeKutta2,   100,  0.1, mylabel='rk2', color='firebrick')
    s.simulate(P(0), Particle.rungeKutta4,   100,  0.1, mylabel='rk4', color='blue')


    s.nextPlot(name="Semi Implicit Euler & Leap Frog after 100 laps")
    s.simulate(P(0), Particle.semiImplicitEuler, 100, 0.1, mylabel='sie', color='purple')
    s.simulate(P(0), Particle.leapFrog,          100, 0.1, mylabel=' lf', color='orange')


    s.nextPlot(name="RK4 with different e")
    s.simulate(P(0.00), Particle.rungeKutta4, 1, 0.01, mylabel='rk4', color='orchid')
    s.simulate(P(0.50), Particle.rungeKutta4, 1, 0.01, mylabel='rk4', color='red')
    s.simulate(P(0.90), Particle.rungeKutta4, 1, 0.01, mylabel='rk4', color='green')
    s.addPlot([0,2],[0,2])
    s.addPlotApocenter(.9, color='green', marker='^')
    s.addPlotApocenter(.5, color='red', marker='^')

    s.nextPlot(name="RK4 base, comp. RK4, RK2, EE with e=0.5")
    s.simulate(P(0.50), Particle.rungeKutta4,      1  , 0.01,   label='rk4-base', color='black', marker='+', markevery=30)
    s.simulate(P(0.50), Particle.rungeKutta4,      1  , 0.40, mylabel='rk4',      color='blue')
    s.simulate(P(0.50), Particle.rungeKutta2,      1.8, 0.40, mylabel='rk2',      color='firebrick')
    s.simulate(P(0.50), Particle.explicitEuler,    1.8, 0.04, mylabel=' ee',      color='green')
    s.addPlotApocenter(.5, color='black', marker='^')


    s.nextPlot(name="RK4 base, comp. SIE, LF with e=0.5")
    s.simulate(P(0.50), Particle.rungeKutta4,      1  , 0.01,   label='rk4-base', color='black', marker='+', markevery=30)
    s.simulate(P(0.50), Particle.semiImplicitEuler,1.8, 0.40, mylabel='sie',      color='purple')
    s.simulate(P(0.50), Particle.leapFrog,         1.8, 0.40, mylabel=' lf',      color='orange')
    s.addPlotApocenter(.5, color='black', marker='^')


    s.plot(showPerfect=False)