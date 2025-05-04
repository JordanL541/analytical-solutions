import sympy as sy

gdot, ux, uy, theta, rho, tau, x, y, x0, y0, vix, viy = sy.symbols("gdot ux uy theta rho tau x y x0 y0 vix viy") #Import variables for symbolic manipulation.

ux = gdot * (y * sy.cos(theta)*sy.cos(theta) - x * sy.sin(theta)*sy.cos(theta))  #Define the velocity profile.
uy = gdot * (y * sy.sin(theta)*sy.cos(theta) - x * sy.sin(theta)*sy.sin(theta))

vu = vix*ux + viy*uy  #Define variables to make code more readable.
uu = ux*ux + uy*uy
vixux3vix = vix - ux + 3*vix*vu
viyuy3viy = viy - uy + 3*viy*vu


f0 = (2/3. - 0.5 * vix**2)*(2/3. - 0.5 * viy**2) * rho * (1 - 1.5 * uu + 3*vu + 4.5 * vu * vu) #Define f0, f1, f2 as given in the paper.

f1 = (2/3. - 0.5 * vix**2)*(2/3. - 0.5 * viy**2) * rho * (1.5 * gdot * tau * sy.sin(2 * theta)*(vix * vixux3vix - viy*viyuy3viy)
                     +3 * gdot * tau * vix * viyuy3viy * sy.sin(theta) * sy.sin(theta)
                     -3 * gdot * tau * viy * vixux3vix * sy.cos(theta) * sy.cos(theta)
                     )

f2 = 3*(2/3. - 0.5 * vix**2)*(2/3. - 0.5 * viy**2)*gdot**2*tau*(2*tau-1)*rho*(sy.sin(theta)**2*sy.cos(theta)**2 * (vix**4+viy**4-6*(vix*viy)**2)
        + (sy.sin(theta))**3 * sy.cos(theta) * (3*vix*viy*(vix**2-viy**2) + vix*viy)
        + sy.sin(theta) *(sy.cos(theta))**3 * (3*vix*viy*(viy**2-vix**2)+vix*viy)
        + (sy.sin(theta))**4/2 * (3 * (vix*viy)**2 - vix**2)
        + (sy.cos(theta))**4/2 * (3*(vix*viy)**2 - viy**2)                
      )


f = f0 + f1 + f2 #Define f symbolically

fcenter = f.subs([(x, x0), (y,y0)]) #Define f at a certain lattice point (x0, y0).

fneighbor = f.subs([(x, x0+vix), (y, y0+viy)]) #Define this as in the Lattice Boltzmann Equation.

myterm = (fneighbor - fcenter) * tau + fcenter - f0.subs([(x,x0), (y,y0)])  #Rewriting of the Lattice Boltzmann Equation. Want to show this = 0.
print(sy.trigsimp(sy.expand(myterm))) #We see there is a prefactor that is a polynomial in vix, viy. This turns out to be 0 as shown below. We use sy.trigsimp instead of sy.simplify because sy.trigsimp is well-defined and won't change in future updates.
print("\n")
myvariable =-0.75*tau*vix**6*viy**2 + 1.0*tau*vix**6 + 1.75*tau*vix**4*viy**2 - 2.33333333333333*tau*vix**4 - 0.75*tau*vix**2*viy**6 + 1.75*tau*vix**2*viy**4 - 2.0*tau*vix**2*viy**2 + 1.33333333333333*tau*vix**2 + 1.0*tau*viy**6 - 2.33333333333333*tau*viy**4 + 1.33333333333333*tau*viy**2 + 0.375*vix**6*viy**2 - 0.5*vix**6 - 0.875*vix**4*viy**2 + 1.16666666666667*vix**4 + 0.375*vix**2*viy**6 - 0.875*vix**2*viy**4 + 1.0*vix**2*viy**2 - 0.666666666666667*vix**2 - 0.5*viy**6 + 1.16666666666667*viy**4 - 0.666666666666667*viy**2
#print(myvariable)
print(myvariable.subs([(vix**6, vix**2), (vix**4, vix**2), (viy**6, viy**2), (viy**4, viy**2)])) #These simplifying substitutions show the prefactor is 0.



