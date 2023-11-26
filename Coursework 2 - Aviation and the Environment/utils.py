import numpy as np
import matplotlib.pyplot as plt


class Aircraft:
    def __init__(self, specs):
        # Sea-level atmospheric conditions
        self.T_sl = 288.15
        self.p_sl = 101.325
        self.rho_sl = 1.225
        self.a_sl = 340.3

        # Aircraft weights
        self.W_E = specs['W_E']
        self.W_MP = specs['W_MP']
        self.W_MTO = specs['W_MTO']

        # Other aircraft specs
        self.maxPass = specs['maxPass']
        self.MPRange = specs['MPRange']
        self.W_F_MP = specs['W_F_MP']
        self.SWing = specs['SWing']

        # Engine specs at cruise
        self.rEngine = specs['rEgnine']
        self.thetEngine = specs['thetEngine']
        self.effComp = specs['effComp']
        self.effTurb = specs['effTurb']
        self.FPREngine = specs['FPREngine']
        self.effFan = specs['effFan']
        self.effTransfer = specs['effTransfer']
        
        # Empirical relations parameters
        self.k = specs['k'] # Breguet range take-off fuel burn offset
        self.c1 = specs['c1'] # weight correlation
        self.c2 = specs['c2']
        self.K1 = specs['K1'] # Parabolic drag model parameters
        self.K2 = specs['K2']

        # Variables to change for analysis
        self.vCruise = specs['vCruise']
        self.MCruise = specs['MCruise']
        self.initCruiseAlt = specs['initCruiseAlt']
        self.LDCruise = specs['LDCruise']

    def ISACond(self, h):
        # International Standard Atmosphere (slide 21)

        if h <= 11:
            T = self.T_sl - 6.5*h
            TRat = T/self.T_sl

            p = TRat**5.256 * self.p_sl
            rho = TRat**4.256 * self.rho_sl
            a = TRat**0.5 * self.a_sl

        elif h<=20:
            T = 216.65

            tropRat = np.exp(-0.1577*(h-11))
            p = tropRat * 22.631 # Tropopause values
            rho = tropRat * 0.364
            a = (T/self.T_sl)**0.5 * 340.3
        
        return T, p, rho, a
    
    def betaParabolic(self, Cl=-1, opt=False):
        # Beta (1/(L/D)) for given Cl using the parabolic drag model (slide 31)

        if opt==False:
            if Cl==-1:
                raise Exception('No Cl input given')
            return self.K1/Cl + self.K2*Cl
        
        return 2*(self.K1*self.K2)**0.5

    def EAS(self, rho=-1, opt=False, W=-1):
        # Equivalent air speed (slide 33)
        if opt==False:
            if rho==-1:
                raise Exception('No density input given')
            return (rho/self.rho_sl)
        
        if W==-1:
            raise Exception('No weight/wing surface area input given')
        
        return (W/0.5/self.rho_sl/self.SWing)**0.5 * (self.K2/self.K1)**0.25
    
    def betaSpeedRatio(self, rho, W):
        # Beta given aircraft weight and wing surface area
        Bstar = self.betaParabolic(opt=True)
        Ve = self.EAS(rho)
        Vestar = self.EAS(opt=True, W=W)
        vRat = Ve/Vestar

        return 0.5*Bstar*(vRat**2 + 1/vRat**2)


if __name__ == '__main__':
    
    ## Plotting temperature, pressure, density and sonic speed ratios over different altitudes
    # hRange = np.arange(0,20,0.01)

    # TArray, pArray, rhoArray, aArray = np.zeros(len(hRange)), np.zeros(len(hRange)), np.zeros(len(hRange)), np.zeros(len(hRange))

    # for i, h in enumerate(hRange):
    #     TRat, pRat, rhoRat, aRat = aircraft.ISACond(h)

    #     TArray[i], pArray[i], rhoArray[i], aArray[i] = TRat/288.15, pRat/101.325, rhoRat/1.225, aRat/340.3

    # plt.plot(hRange, TArray)
    # plt.plot(hRange, pArray)
    # plt.plot(hRange, rhoArray)
    # plt.plot(hRange, aArray)

    # plt.show()

    aircraftSpecs = {'W_E':106_000, 'W_MP':40_000, 'W_MTO':220_000, 'maxPass':240, 'MPRange':12_000, 'W_F_MP':74_000, 'SWing':315, 'rEgnine':45, 
    'thetEngine':6, 'effComp':0.9, 'effTurb':0.9, 'FPREngine':1.45, 'effFan':0.92, 'effTransfer':0.9, 'k':0.015, 'c1':0.3, 'c2':1.0, 'K1':0.0125, 'K2':0.0446, 
    'vCruise':256, 'MCruise':0.85, 'initCruiseAlt':9.449, 'LDCruise':21}

    aircraft = Aircraft(aircraftSpecs)
    