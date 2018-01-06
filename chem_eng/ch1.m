/*
Edited by Hajime Fujita
Reference: "基礎式から学ぶ化学工学"（2017.11, 化学同人）
*/

/* Parameter and Constant */
V: Volume of tank
c_A0: Concentration of initial flow
F: Flow
T: Temperature
r_A: ?
p: density (kg/m.^3)
C_p: Heat capacity of flow

/* Equation */

V*diff(c_A0/t) = F*c_A0 - F*c_A + V*r_A /*Macro-balance of component A*/

V*p*C_p*diff(T, t) = F*p*C_p*T_0 - F*p*C_p*T + V*Q - A*q /*Macro-balance of energy*/

diff(h, t) = (1./A)*F - (1/(AR))*h

/*
*/
