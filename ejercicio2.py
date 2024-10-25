#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sympy as sp
from sympy import symbols, Matrix
from sympy import init_printing
from sympy import simplify
from sympy import symbols, expand, factor
from pytc2.general import to_latex

from schemdraw import Drawing

# Ahora importamos las funciones de PyTC2

from pytc2.remociones import remover_polo_dc, remover_polo_jw , remover_polo_jw2 ,remover_polo_infinito , remover_polo_dc2 , remover_polo_infinito2
from pytc2.dibujar import display, dibujar_tanque_serie, dibujar_puerto_entrada, dibujar_funcion_exc_abajo,  dibujar_elemento_serie, dibujar_elemento_derivacion,  dibujar_tanque_derivacion, dibujar_tanque_RC_serie,  dibujar_espacio_derivacion, Capacitor, Resistor, ResistorIEC
from pytc2.dibujar import dibujar_Pi, dibujar_Tee, dibujar_lattice

from pytc2.general import print_latex, print_subtitle, a_equal_b_latex_s
from IPython.display import display,  Markdown

from pytc2.cuadripolos import Z2Tabcd_s, Y2Tabcd_s, Tabcd2Z_s, Tabcd2Y_s
from pytc2.cuadripolos import calc_MAI_impedance_ij, calc_MAI_vtransf_ij_mn, calc_MAI_ztransf_ij_mn

import scipy.signal as sig
from pytc2.sistemas_lineales import analyze_sys



patron = "\n" + "/" * 75 + "\n" # 75 barras para ajustar el largo deseado
# Activar la impresión en formato LaTeX
init_printing()
# Definir la variable simbólica s
s = sp.symbols('s',complex=True)
R = sp.symbols('R')

# Definir la función de transferencia H(s)
numerador = (s**2 + 9)
denominador = s*(s**2 + 2)
Z21 = numerador / denominador

print_latex(a_equal_b_latex_s('Z21(s)', Z21))

##/////////////////////////////////////////////////////////////////////////////

numerador = (2*s**2 + 1)
denominador = s*(s**2 + 2)
Z22= numerador / denominador

print_latex(a_equal_b_latex_s('Z22(s)', Z22))
print(patron)

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

Y2, Yt1  = remover_polo_infinito(1/Z22, omega_zero=3)
C1 = Yt1/s
C1 = sp.nsimplify(C1)
Y2 = sp.nsimplify(sp.factor(Y2))
Yt3 = sp.nsimplify(sp.factor(Yt1))

print_latex(a_equal_b_latex_s('C1', C1))
print_latex(a_equal_b_latex_s('Yt1(s)', Yt1))
print_latex(a_equal_b_latex_s('Y2(s)', Y2))
print(patron)

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

Z4, Zt3 , L2, C2 = remover_polo_jw(1/Y2, omega=3, isImpedance=True)
L2 = sp.nsimplify(L2)
C2 = sp.nsimplify(C2)
Z4 = sp.nsimplify(sp.factor(Z4))
Zt3 = sp.nsimplify(sp.factor(Zt3))

print_latex(a_equal_b_latex_s('L2', L2))
print_latex(a_equal_b_latex_s('C2', C2))
print_latex(a_equal_b_latex_s('Zt3(s)', Zt3))
print_latex(a_equal_b_latex_s('Z4(s)', Z4))
print(patron)

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

Y6, Yt5  = remover_polo_infinito(1/Z4 , isSigma=False)
C3 = Yt5/s
C3 = sp.nsimplify(C3)
Y6 = sp.nsimplify(sp.factor(Y6))
Yt5 = sp.nsimplify(sp.factor(Yt5))

print_latex(a_equal_b_latex_s('C3', C3))
print_latex(a_equal_b_latex_s('Yt5(s)', Yt5))
print_latex(a_equal_b_latex_s('Y6(s)', Y6))

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////

# Dibujamos la red resultante:

d = dibujar_puerto_entrada(Drawing(unit=4),
                        voltage_lbl = ('+', '$V$', '-'), 
                        current_lbl = '$I$')


d = dibujar_elemento_derivacion(d, 'C', C3)

d = dibujar_tanque_serie(d, L2, C2 )

d = dibujar_elemento_derivacion(d, 'C', C1)

display(d)

##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
##/////////////////////////////////////////////////////////////////////////////
print("\n" + "/" * 75)
print("/" * 75 + "\n")

print('Calculos de verificacion')

Ya , Yb , Yc = sp.symbols('Ya , Yb , Yc', complex=True)


YY = sp.Matrix([[ Ya + Yb , -Yb ], [ -Yb , Yb + Yc ]])


print_latex(a_equal_b_latex_s('YY{1}', YY))

print(patron)

YY_sym = sp.simplify(YY.subs(Ya, s*C1 ) ) 
YY_sym = sp.simplify(YY_sym.subs(Yb, 1/(s*L2) + C2*s ))
YY_sym = sp.simplify(YY_sym.subs(Yc, s*C3 ) )

print_latex(a_equal_b_latex_s('YY_{T}', YY_sym) )

TT = Y2Tabcd_s (YY_sym)

ZZ = Tabcd2Z_s (TT)

Z11 = sp.factor(sp.simplify(sp.expand(ZZ[0])))
Z12 = sp.factor(sp.simplify(sp.expand(ZZ[1])))
Z21 = sp.factor(sp.simplify(sp.expand(ZZ[2])))
Z22 = sp.factor(sp.simplify(sp.expand(ZZ[3])))

print(patron)

print_latex(a_equal_b_latex_s('ZZ_{11}', Z11 ) )

print_latex(a_equal_b_latex_s('ZZ_{12}', Z12 ) )

print_latex(a_equal_b_latex_s('ZZ_{21}', Z21 ) )

print_latex(a_equal_b_latex_s('ZZ_{22}', Z22) )

print_latex(a_equal_b_latex_s('\\frac{V_2}{I_1}(s)', sp.factor(Z21/(1+Z22))))

H = sig.TransferFunction( [1, 0, 9], [1, 2, 2, 1] )

analyze_sys(H)