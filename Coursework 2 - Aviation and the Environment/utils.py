import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class Aircraft:
    def __init__(self, specs):
        # Sea-level atmospheric conditions
        self.T_sl = 288.15
        self.p_sl = 101325
        self.rho_sl = 1.225
        self.a_sl = 340.3

        # Current flight condition values
        self.h = -1
        self.s = 0.0 # Air distance travelled
        self.T = -1
        self.T0 = -1
        self.p = -1
        self.p0 = -1
        self.rho = -1
        self.rho0 = -1
        self.a = -1

        self.V = -1
        self.M = -1
        self.Cl = -1
        self.beta = -1
        self.EAS = -1
        self.EASOpt = -1

        # Aircraft weights
        self.W_E = specs['W_E']
        self.W_MP = specs['W_MP']
        self.W_MTO = specs['W_MTO']
        self.W_P = specs['W_P']
        self.W_F = specs['W_F']
        self.W_FB = 0.0 # Weight of fuel burnt

        self.W = self.W_E + self.W_P + self.W_F
        self.W_TO = self.W

        self.D = -1

        # Other aircraft specs
        self.maxPass = specs['maxPass']
        self.MPRange = specs['MPRange']
        self.W_F_MP = specs['W_F_MP'] # Fuel capacity at maximum payload - WHY DOES PAYLOAD AFFECT FUEL CAPACITY????
        self.SWing = specs['SWing']

        # Engine specs at cruise
        self.rEngine = specs['rEgnine']
        self.thetEngine = specs['thetEngine']
        self.effComp = specs['effComp']
        self.effTurb = specs['effTurb']
        self.FPREngine = specs['FPREngine']
        self.effFan = specs['effFan']

        # Engine parameters
        self.LCVFuel = specs['LCVFuel']
        self.H = -1 # Range parameter

        # Overall efficiency
        self.effProp = -1 # To be calculated later on
        self.effCycle = -1
        self.effTransfer = specs['effTransfer']
        self.effOverall = -1

        # Pollution parameters
        self.W_FBTimestep = -1
        self.W_CO2 = 0.0
        self.W_NOX = 0.0
        
        # Model parameters
        self.timestep = specs['timestep']
        self.time = 0.0

        self.k = specs['k'] # Breguet range take-off fuel burn offset
        self.c1 = specs['c1'] # weight correlation
        self.c2 = specs['c2']
        self.K1 = specs['K1'] # Parabolic drag model parameters
        self.K2 = specs['K2']

        self.vRat = specs['vRat'] # EAS ratio (design point)

        self.betaOpt = 2*(self.K1*self.K2)**0.5

    def checkVariables(self):
        if self.time == 0.0:
            # Initialisation checks that only have to happen the first iteration

            # Engine values
            assert self.thetEngine <= 6
            assert self.rEngine <= 50
            assert self.effComp <= 0.92
            assert self.effTurb <= 0.92
            assert self.effFan <= 0.94
            assert self.FPREngine >= 1.35

            # Weight constraints
            # assert self.W_F <= self.W_F_MP
            assert self.W_P <= self.W_MP
            assert self.W <= self.W_MTO

        if self.W_F < 0:
            raise Exception(f'Zero fuel remaining at air distance: {self.s}')
        if self.M > 0.85:
            raise Exception(f'Mach number exceeding transonic drag threshold of 0.85')
        
        return
    
    #########################################################################################################################################################
    #                                                                   AIRCRAFT
    #########################################################################################################################################################

    def ISACondUpdate(self):
        # International Standard Atmosphere (slide 21)
        if self.h <= 11:
            self.T = self.T_sl - 6.5*self.h
            TRat = self.T/self.T_sl

            self.p = TRat**5.256 * self.p_sl
            self.rho = TRat**4.256 * self.rho_sl
            self.a = TRat**0.5 * self.a_sl
        elif self.h <= 20:
            self.T = 216.65

            tropRat = np.exp(-0.1577*(self.h-11))
            self.p = tropRat * 22631 # Tropopause values
            self.rho = tropRat * 0.364
            self.a = (self.T/self.T_sl)**0.5 * 340.3
        else:
            raise Exception('Incorrect altitude')
        return

    def speedsUpdate(self):
        # Optimal equivalent air speed for given aricraft weight (slide 33)
        self.EASOpt = (self.W*9.81/0.5/self.rho_sl/self.SWing)**0.5 * (self.K2/self.K1)**0.25
        self.EAS = self.vRat * self.EASOpt

        # Update true air speed (slide 33)
        self.V = self.EAS / (self.rho/self.rho_sl)**0.5
        
        # AIR Mach number from velocity and temperature
        self.M = self.V / self.a
        return
    
    def speedsOvrd(self, V, M):
        if V==-1:
            self.M = M
            self.V = M*self.a
        elif M==-1:
            self.V = V
            self.M = V/self.a
        else:
            self.V = V
            self.M = M
        
        self.EASOpt = (self.W/0.5/self.rho_sl/self.SWing)**0.5 * (self.K2/self.K1)**0.25
        self.EAS = self.V * (self.rho/self.rho_sl)**0.5
        self.vRat = self.EAS/self.EASOpt
        return

    def ClUpdate(self):
        self.Cl = self.W * 9.81 / (0.5 * self.rho * self.V**2 * self.SWing)
        return
    
    def DBetaUpdate(self):
        # Beta (1/(L/D)) and drag for given Cl using the parabolic drag model (slide 31)

        self.beta = self.K1/self.Cl + self.K2*self.Cl
        # self.beta = 0.5*self.betaOpt*(self.vRat**2 + 1/self.vRat**2)
        
        self.D = self.W * self.beta
        return

    #########################################################################################################################################################
    #                                                                   ENGINE
    #########################################################################################################################################################
    
    def stagnationValuesUpdate(self):
        # Calculate stagnation temperature, pressure and density for engine computations
        insideBrackets = 1 + 0.2*self.M**2

        self.T0 = self.T / insideBrackets**-1
        self.p0 = self.p / insideBrackets**(-1.4/0.4)
        self.rho0 = self.rho / insideBrackets**(-1/0.4)
        return

    def engineEfficienciesUpdate(self):
        # Calculate the overall engine efficiency given the current flight conditions (slide 39)

        # Cycle efficiency
        gamRat = 0.4/1.4
        self.effCycle = (self.thetEngine * self.effTurb * (1-1/self.rEngine**gamRat) - (self.rEngine**gamRat - 1)/self.effComp) / (self.thetEngine - 1 - (self.rEngine**gamRat - 1)/self.effComp)
        
        # Propulsion efficiency``
        Mjet = (5*((self.FPREngine * self.p0/self.p)**gamRat - 1))**0.5
        TjRat = (1 + 0.2*self.M**2) / (1 + 0.2*Mjet**2) * self.FPREngine**(gamRat/self.effFan)
        self.effProp = 2 / (1 + Mjet * TjRat**0.5 / self.M)

        # Overall efficiency
        self.effOverall = self.effProp * self.effCycle * self.effTransfer

        return
    
    def rangeParamUpdate(self):
        self.H = self.LCVFuel * self.effOverall / 9.81 / self.beta
        return
    
    def sFuelburnUpdate(self):
        # Calculate the distance travelled this timestep using the speed
        # Calculate the change in fuel weight during this timestep and update all aircraft weights acoordingly

        if self.time == 0.0:
            self.W_FBTimestep = 0.015 * self.W_TO
        else:
            sTimestep = self.V * self.timestep
            self.s += sTimestep
            self.W_FBTimestep = self.W * (1-np.exp(-sTimestep/self.H))

        self.W_F -= self.W_FBTimestep
        self.W -= self.W_FBTimestep
        self.W_FB += self.W_FBTimestep
        return
    
    def pollutionUpdate(self):
        # NOX emissions
        T03 = 1 + (self.rEngine**(0.4/1.4) - 1)/self.effComp
        EI_NOX = 0.011445 * np.exp(0.00676593 * T03) * self.T0
        self.W_NOX += EI_NOX * self.W_FBTimestep * (15.1*2)

        # CO2 emissions
        self.W_CO2 += 3088 * self.W_FBTimestep

        # Water vapour - contrails

        return
    
    #########################################################################################################################################################
    #                                                                   SIMULATION
    #########################################################################################################################################################

    def updateAllFlightValues(self, h, ft=True, MOvrd=-1, VOvrd=-1):
        # Update altitude before everything else
        if ft == True:
            # Altitude in ft
            self.h = 0.0003048*h
        else:
            # Altitude in km
            self.h = h
        
        self.checkVariables()

        # Aircraft
        self.ISACondUpdate()

        # if self.time == 0.0:
        self.sFuelburnUpdate()

        if MOvrd!=-1 or VOvrd!=-1:
            self.speedsOvrd(VOvrd, MOvrd)
        else:
            self.speedsUpdate()

        self.ClUpdate()
        self.DBetaUpdate()

        # Engine
        self.stagnationValuesUpdate()
        self.engineEfficienciesUpdate()
        self.rangeParamUpdate()
        self.pollutionUpdate()

        self.time += self.timestep
        return


if __name__ == '__main__':
    
    notesData = '-13.779527559055055, 2508.960573476652; 930.1181102362204, 9677.419354838668; 1860.236220472441, 17275.985663082407; 2321.850393700787, 21003.584229390646; 2790.354330708661, 24587.81362007164; 3727.36220472441, 32329.749103942624; 3861.712598425197, 33548.38709677418; 0, 7956.98924731178; 926.6732283464567, 15985.663082437255; 1853.3464566929133, 24157.706093189918; 2328.740157480315, 28458.78136200717; 2790.354330708661, 32616.4874551971; 3727.36220472441, 41290.32258064514; 4664.370078740157, 50179.21146953403; 5594.488188976376, 58996.415770609296; 6204.232283464566, 65089.605734767'
    notesSFB, notesFB = zip(*[(float(datapoint.split(', ')[0]), float(datapoint.split(', ')[1])) for datapoint in notesData.split('; ')])

    timestep = 1 # s

    # Analysis variables
    vRat = 1.016 # Ratio of EAS to optimal EAS
    initCruiseAlt = 9.5
    W_P = 24_000
    W_F = 90_000

    aircraftSpecs = {'W_E':106_000, 'W_MP':40_000, 'W_MTO':220_000, 'W_P':W_P, 'W_F':W_F, # Weight specs
                     'maxPass':240, 'MPRange':12_000, 'W_F_MP':74_000, 'SWing':315, # Other aircraft specs
                     'rEgnine':45, 'thetEngine':6, 'effComp':0.9, 'effTurb':0.9, 'FPREngine':1.45, 'effFan':0.92, 'effTransfer':0.9, 'LCVFuel':42.7e6, # Engine
                     'timestep':timestep, 'k':0.015, 'c1':0.3, 'c2':1.0, 'K1':0.0125, 'K2':0.0446, 'vRat':vRat} # Model parameters
    

    aircraft = Aircraft(aircraftSpecs)


    # L/D against Cl
    # CLs = []
    # LDs = []
    # MRange = np.arange(0.6,1,0.01)
    # for M in MRange:
    #     aircraft.updateAllFlightValues(initCruiseAlt, ft=False)
    #     aircraft.speedsOvrd(V=-1, M=M)
    #     aircraft.ClUpdate()
    #     aircraft.DBetaUpdate()
    #     CLs.append(aircraft.Cl)
    #     LDs.append(1/aircraft.beta)

    # sns.set_palette('dark')

    # plt.subplot(1,2,1)
    # plt.plot(CLs, LDs)
    # plt.xlabel('C_l')
    # plt.ylabel('L/D')
    # plt.grid()

    # plt.subplot(1,2,2)
    # plt.plot(MRange, CLs)
    # plt.xlabel('Distance travelled, s / km')
    # plt.ylabel('Fuel burnt / kg')
    # plt.grid()
    # plt.show()

    # Fuel burn against distance travelled
    altitudes = initCruiseAlt * np.ones(10)

    distances = []
    fuelBurn = []

    while aircraft.W_F > 0:
        aircraft.updateAllFlightValues(31000, ft=True, MOvrd=-1) #0.85)
        distances.append(aircraft.s/1000)
        fuelBurn.append(aircraft.W_FB)

    plt.plot(distances, fuelBurn, label='fuel burn method 1')
    plt.hlines(W_F, 0, 12000, linestyles='--', color='purple')
    plt.scatter(notesSFB, notesFB, s=3, color='darkgreen')
    plt.yticks(np.arange(0,max(fuelBurn) + 5000,5000))
    plt.grid()
    plt.legend()
    plt.show()
    print('lol')