import numpy as np

R_d = 287.0 # dry ideal gas constant

def read_input(fpath,defaults={},prefix=''):
    """Read parameters with specified prefix from input file"""
    prob_parm = {key:val for key,val in defaults.items()}
    with open(fpath,'r') as f:
        for line in f:
            line = line.strip().split('#')[0] # ignore comments
            if line == '':
                continue
            if line.startswith(prefix):
                #print('parse:',line)
                key, defn = line.split('=')
                key = key[len(prefix):].rstrip()
                value = []
                for val in defn.split():
                    try:
                        val = float(val)
                    except ValueError:
                        val = val.strip().strip('"').strip("'")
                    finally:
                        value.append(val)
                if len(value) == 1:
                    value = value[0]
                prob_parm[key] = value
    return prob_parm


class IsentropicVortex(object):
    """Generate 2D analytical solution"""
    def __init__(self,
                 M_inf=0.0, alpha=0.0, # freestream def
                 p_inf=1e5, T_inf=300.0, gamma=1.4, # fluid params
                 beta=0.05, sigma=1.0, R=2.0, # vortex params
                ):
        """Initialize vortex, can directly use parameters parsed from
        the input file, e.g., `read_input(fpath,prefix='prob.')`
        """
        self.M_inf = M_inf
        self.alpha = alpha
        self.gamma = gamma
        self.beta = beta
        self.sigma = sigma
        self.R = R
        self.p_inf = p_inf
        self.T_inf = T_inf
        print('ref pressure:',self.p_inf)
        print('ref temperature:',self.T_inf)
        # calculate additional reference values
        self.rho_inf = p_inf / (R_d*T_inf)
        print('ref density:',self.rho_inf)
        self.a_inf = np.sqrt(gamma*R_d*T_inf)
        print('ref speed of sound:',self.a_inf)
        
    def evaluate(self,
                 left_edge,right_edge,dims,
                 xc_norm=0.5,yc_norm=0.5):
        """Evaluate on specified grid"""
        Nx,Ny = dims[0],dims[1]
        # handle edge inputs as YTArray
        try:
            left_edge = left_edge.value
        except AttributeError:
            pass
        try:
            right_edge = right_edge.value
        except AttributeError:
            pass
        xv = xc_norm * (right_edge[0]-left_edge[0])
        yv = yc_norm * (right_edge[1]-left_edge[1])
        print(f'vortex at ({xv:g},{yv:g})')
        
        # mesh nodes
        x1n = np.linspace(left_edge[0],right_edge[0],Nx+1)
        y1n = np.linspace(left_edge[1],right_edge[1],Ny+1)
        xxn,yyn = np.meshgrid(x1n,y1n,indexing='ij')
        
        # mesh cell centers
        x1 = (x1n[1:] + x1n[:-1]) / 2
        y1 = (y1n[1:] + y1n[:-1]) / 2
        xx,yy = np.meshgrid(x1,y1,indexing='ij')
        
        # mesh x- and y-face centers
        xx_xf,yy_xf = np.meshgrid(x1n,y1,indexing='ij')
        xx_yf,yy_yf = np.meshgrid(x1,y1n,indexing='ij')
        
        # calculate fields @ cell centers
        dx = (xx - xv) / self.R
        dy = (yy - yv) / self.R
        r2 = dx**2 + dy**2;
        Omg = self.beta * np.exp(-r2/(2*self.sigma**2))
        deltaT = -(self.gamma-1.0)/(2.0*self.sigma**2) * Omg**2
        inv_gm1 = 1.0 / (self.gamma - 1.0)
        rho = (1.0 + deltaT)**inv_gm1
        p = 1.0/self.gamma * rho**self.gamma
        
        # calculate x-face velocities
        dx = (xx_xf - xv) / self.R
        dy = (yy_xf - yv) / self.R
        r2 = dx**2 + dy**2;
        Omg_xf = self.beta * np.exp(-r2/(2*self.sigma**2))
        u = self.M_inf*np.cos(self.alpha) - dy*Omg_xf
        u = (u[:-1,:] + u[1:,:]) / 2
        
        # calculate y-face velocities
        dx = (xx_yf - xv) / self.R
        dy = (yy_yf - yv) / self.R
        r2 = dx**2 + dy**2;
        Omg_yf = self.beta * np.exp(-r2/(2*self.sigma**2))
        v = self.M_inf*np.sin(self.alpha) + dx*Omg_yf
        v = (v[:,:-1] + v[:,1:]) / 2
        
        # dimensionalize w/ characteristic values
        rho *= self.rho_inf
        u *= self.a_inf
        v *= self.a_inf
        p *= self.rho_inf*self.a_inf**2 # note: this is _not_ p_inf!
        T = self.T_inf * (1.0 + deltaT)
        assert np.allclose(p, rho*R_d*T)
        
        return xxn,yyn,rho,u,v,p,T
        
