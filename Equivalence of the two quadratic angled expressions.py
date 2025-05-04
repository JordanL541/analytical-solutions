import sympy as sy

tau, rho, viyprime, gdot, vixprime, vix, viy, theta, ux, uy, uxprime, uyprime, x, y = sy.symbols("tau rho viyprime gdot vixprime vix viy theta ux uy uxprime uyprime x y ") #Define symbols

f0flat = rho*(2/3. - 0.5 * vix**2)*(2/3. - 0.5 * viy**2)*(1-1.5*(uxprime*uxprime) + 3*(vixprime*uxprime) + 4.5 *(vixprime*uxprime)**2) #Define equil. dist. in rotated frame

f1flat = -3*tau*(2/3. - 0.5 * vix**2)*(2/3. - 0.5 * viy**2)*rho*viyprime*gdot*(vixprime + 3*uxprime*vixprime**2 - uxprime) #Define f1 in rotated frame

f2flat = (2/3. - 0.5 * vix**2)*(2/3. - 0.5 * viy**2)*tau*(tau-0.5)*rho*viyprime**2*gdot**2*3*(3*vixprime**2-1) #Define f2 in rotated frame

fflat =  f0flat+f1flat + f2flat #The solution in the rotated frame

fflatrotation = fflat.subs([(viyprime, -vix*sy.sin(theta) + viy*sy.cos(theta)), (vixprime, vix*sy.cos(theta) + viy*sy.sin(theta)), (uxprime, ux*sy.cos(theta)+uy*sy.sin(theta))]) #The solution in the standard frame now

vu = vix*ux + viy*uy #Define helpful variables
uu = ux*ux + uy*uy
vixux3vix = vix - ux + 3*vix*vu
viyuy3viy = viy - uy + 3*viy*vu

f0 = (2/3. - 0.5 * vix**2)*(2/3. - 0.5 * viy**2) * rho * (1 - 1.5 * uu + 3*vu + 4.5 * vu * vu) #Define f0, f1, f2 as given in the appendix

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

mydiff = f0+f1+f2 - fflatrotation #Look at the difference
mydiff2 = mydiff.subs([(ux,gdot * (y * sy.cos(theta)*sy.cos(theta) - x * sy.sin(theta)*sy.cos(theta))), (uy, gdot * (y * sy.sin(theta)*sy.cos(theta) - x * sy.sin(theta)*sy.sin(theta)))]) #Substitute in the flow profile
mydiff3 = sy.trigsimp(mydiff2) #Simplify
mydiff4 = mydiff3.subs([(vix**4, vix**2), (vix**6, vix**2), (viy**4, viy**2), (viy**6, viy**2)]) #Input simplifying relations
print(sy.simplify(mydiff4))





