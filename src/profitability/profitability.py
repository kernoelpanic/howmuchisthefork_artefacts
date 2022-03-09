#!/usr/bin/env python3
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import itertools
import scipy.special
from decimal import *
import unittest as ut
from sympy.solvers import solve
from sympy import Symbol
from sympy import zoo
from scipy.constants import golden as phi

import sys
sys.version
#extract source to .py file later:
def Pr_success_join(p_B=0,p_E=0,p_A=None,p_V=0,p_I=0,p_m=0,z=1):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m ) 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    q = p_B + p_E + p_m
    p = p_A + p_V + p_I
    if (q >= 0.5):
        return 1.0
    return ( q/(1-q) )**(z+1)
#extract source to .py file later:
def Pr_success_abstain(p_B=0,p_E=0,p_A=None,p_V=0,p_I=0,p_m=0,z=1):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m ) 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    q = p_B + p_E 
    p = p_A + p_V + p_I + p_m
    if (q >= 0.5):
        return 1.0
    return ( q/(1-q) )**(z+1)
#extract to .py file later:
def EV_main(p_B=0,p_E=0,p_A=None,p_V=0,p_I=0,p_m=0,z=1):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m ) 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    
    P = Pr_success_abstain(p_B=p_B,
                           p_E=p_E,
                           p_A=p_A,
                           p_V=p_V,
                           p_I=p_I,
                           p_m=p_m,
                           z=z)
    return ( 
              ( 
                  (( 1 - P ) * p_m )
                  / ( p_A + p_V + p_I + p_m ) 
              ) 
           )
#extract to .py file later
def EV_fork(p_B=0,p_E=0,p_A=None,p_V=0,p_I=0,p_m=0,z=1,epsilon=0):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m ) 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    
    P = Pr_success_join(p_B=p_B,
                        p_E=p_E,
                        p_A=p_A,
                        p_V=p_V,
                        p_I=p_I,
                        p_m=p_m,
                        z=z)
    return (
                (
                   ( 
                      (P * p_m)/
                      ( p_B + p_E + p_m )
                   ) * ( epsilon + 1 ) 
                )
           )
#extract to .py file later
def e_min(p_B=0,
          p_E=0,
          p_A=None,
          p_V=0,
          p_I=0,
          p_m=0,
          z=1):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m )
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    
    P_join = Pr_success_join(p_B=p_B,
                         p_E=p_E,
                         p_A=p_A,
                         p_V=p_V,
                         p_I=p_I,
                         p_m=p_m,
                         z=z)
    
    P_abstain = Pr_success_abstain(p_B=p_B,
                           p_E=p_E,
                           p_A=p_A,
                           p_V=p_V,
                           p_I=p_I,
                           p_m=p_m,
                           z=z)
    """
    return (
             (  (1 - ( Pa))/
                (
                 (p_A + p_V + p_I + p_m)
                )
             )*
             (
                 (p_B + p_E + p_m)/
                 (
                     Pj
                 ) 
             ) - 1
           )
    """
    return (
             (  
                ((1 - P_abstain) * p_m)/
                (p_A + p_V + p_I + p_m)
             )/
             (
                 (P_join * p_m)/
                 (p_B + p_E + p_m)
             ) - 1
           )
#extract to .py file later
def EV_main_blocks(p_B=0,
                   p_E=0,
                   p_A=None,
                   p_V=0,
                   p_I=0,
                   p_m=0,
                   z=1,
                   n_main=0,
                   n_fork=0):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m ) 
    if p_A < 0:
        p_A = 0 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    assert z >= 0
    assert n_main >= 0
    
    P = Pr_success_abstain(p_B=p_B,
                           p_E=p_E,
                           p_A=p_A,
                           p_V=p_V,
                           p_I=p_I,
                           p_m=p_m,
                           z=z)
    return ( 
              ( 
                  (( 1 - P ) * p_m )
                  / ( p_A + p_V + p_I + p_m ) 
              ) + ( n_main * ( 1 - P ) ) + ( n_fork * P )
           )
#extract to .py file later
def EV_fork_blocks(p_B=0,
                    p_E=0,
                    p_A=None,
                    p_V=0,
                    p_I=0,
                    p_m=0,
                    z=1,
                    epsilon=0,
                    n_fork=0,
                    n_main=0):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m )
    if p_A < 0:
        p_A = 0 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    
    P = Pr_success_join(p_B=p_B,
                        p_E=p_E,
                        p_A=p_A,
                        p_V=p_V,
                        p_I=p_I,
                        p_m=p_m,
                        z=z)
    return (
                (
                   ( 
                      ( P  *  p_m ) / ( p_B + p_E + p_m ) 
                   ) * ( epsilon + 1 )
                ) + ( n_fork * P * ( epsilon + 1 ) ) + ( n_main * ( 1 - P ) )
            )
#extract to .py file later
def EV_fork_blocks_effort(p_B=0,
                          p_E=0,
                          p_A=None,
                          p_V=0,
                          p_I=0,
                          p_m=0,
                          z=1,
                          epsilon=0,
                          n_fork=0,
                          n_main=0):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m )
    if p_A < 0:
        p_A = 0 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    
    P = Pr_success_join(p_B=p_B,
                        p_E=p_E,
                        p_A=p_A,
                        p_V=p_V,
                        p_I=p_I,
                        p_m=p_m,
                        z=z)
    
    return (
                (
                   ( 
                      ( P  *  p_m ) / ( p_B + p_E + p_m ) 
                  #) * ( epsilon )
                   ) * ( 1 + epsilon )
                ) + ( n_fork * P * ( epsilon ) )
                  #+ (((1-P)*p_m)/(p_B + p_E + p_m))
                  #+ ((P_attack_chain_block * p_m)/(p_B + p_E + p_m))
                  + n_fork + (n_main * (1 - P))
            )
#extract to .py file later
def EV_main_blocks_effort(p_B=0,
                   p_E=0,
                   p_A=None,
                   p_V=0,
                   p_I=0,
                   p_m=0,
                   z=1,
                   n_main=0,
                   n_fork=0):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m ) 
    if p_A < 0:
        p_A = 0 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    assert z >= 0
    assert n_main >= 0
    
    P = Pr_success_abstain(p_B=p_B,
                           p_E=p_E,
                           p_A=p_A,
                           p_V=p_V,
                           p_I=p_I,
                           p_m=p_m,
                           z=z)
    return ( 
              ( 
                  (( 1 - P ) * p_m )
                  / ( p_A + p_V + p_I + p_m ) 
              ) + ( n_main * ( 1 - P ) ) + ( n_fork ) 
           )
#extract to .py file later:
def Pr_success(kr=1,
               kl=6,
               N=9,
               p_B=0,
               p_A=None,
               p_E=0,
               p_V=0,
               p_I=0,
               z=0, # This is the start state, also negative values work
               pprint=False):
    
    if p_A is None:
        p_A = 1 - (p_B + p_E  + p_V + p_I)
    #assert p_A > 0
    if p_A < 0.0:
        p_A = 0.0 
    assert math.isclose(p_B + p_A + p_E  + p_V + p_I, 1)
    if kl < z:
        # catch the case where the game is already lost
        return {"Pclearwin":0.0,"Pclearlose":1.0,"Pwin":0.0,"Plose":1.0,"Pdraw":0.0}
    
    p_B = Decimal(p_B)
    p_E = Decimal(p_E)
    p_V = Decimal(p_V)
    p_I = Decimal(p_I)
    p_A = Decimal(p_A)
    
    length = kr + kl + 1
    start = kl - z
    
    # start state 
    M_init = np.zeros(shape=(1,length))
    M_init[0][start] = 1 
    
    # transition matrix of markov chain
    M = np.zeros(shape=(length,length))
    
    # final states S_delta_kr and S_-delta_kl
    M[length-1][length-1] = 1 # stay in final attack win state
    M[0][0] = 1 # stay in final main win state 
    
    # All left states and middle state, equal length
    for s in range(1,kl+1):
        M[s][s-1] = p_A + p_V + p_I
        M[s][s+1] = p_B + p_E
   
    # All right states 
    for s in range(kl+1,kr + kl):
        M[s][s-1] = p_V
        M[s][s+1] = p_B + p_E + p_A + p_I

    # convert array to matrix
    M = np.asmatrix(M)
    M_init = np.asmatrix(M_init)
    
    P = (M_init * M**int(N))
    
    Pwin = 0
    for s in range( kl+1 ,length):
        Pwin += P.item(s)
    Plose = 0
    for s in range(0, kl ):
        Plose += P.item(s)
    Pdraw = P.item( kl )
    
    if pprint:
        print(M)
        print(M_init)
        print(P)
        #print("Pwin  = ",Pwin)
        #print("Pdraw = ",Pdraw)
        #print("Plose = ",Plose)    
        #print("Total = ",Pwin + Pdraw + Plose)
    
    assert math.isclose(Pwin + Pdraw + Plose,1)
    return {"Pclearwin":P.item(length-1),"Pclearlose":P.item(0),"Pwin":Pwin,"Plose":Plose,"Pdraw":Pdraw}
#extract to .py file later:
def EV_main_markov_blocks(p_B=0,
                          p_E=0,
                          p_A=None,
                          p_V=0,
                          p_I=0,
                          p_m=0,
                          z=1,
                          n_main=0,
                          kr=1,
                          kl=6,
                          N=6,
                          Pm_mode="V",
                          Prslt="Pwin",
                          n_fork=0,
                          epsilon=0):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m ) 
    if p_A < 0:
        p_A = 0 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    assert z >= 0
    assert n_main >= 0
    
    # either p_m is part of p_V, or of p_I
    assert (Pm_mode == "V" or Pm_mode == "I")
    if Pm_mode == "V":
        P = Pr_success(kr=kr, 
                       kl=kl, 
                       N=N,
                       z=z,
                       p_B=p_B,
                       p_A=p_A, 
                       p_I=p_I,
                       p_E=p_E,
                       p_V=p_V + p_m)[Prslt] # here it is part of p_V
    elif Pm_mode == "I":
        P = Pr_success(kr=kr, 
                       kl=kl, 
                       N=N,
                       z=z,
                       p_B=p_B,
                       p_A=p_A, 
                       p_I=p_I + p_m, # here it is part of p_I
                       p_E=p_E,
                       p_V=p_V)[Prslt] 

    assert N > 0 # because if N == 0 then the attack is over
    return ( 
              ( 
                  (( 1 - P ) * p_m )
                  / ( p_A + p_V + p_I + p_m ) 
              ) + ( n_main * ( 1 - P ) ) + ( P * n_fork * (epsilon+1) )
           )
#extract to .py file later:
def EV_fork_markov_blocks(p_B=0,
                          p_E=0,
                          p_A=None,
                          p_V=0,
                          p_I=0,
                          p_m=0,
                          z=1,
                          epsilon=0,
                          n_fork=0,
                          kr=1,
                          kl=6,
                          N=6,
                          Prslt="Pwin",
                          n_main=0):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m )
    if p_A < 0:
        p_A = 0 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    

    P = Pr_success(kr=kr, 
                       kl=kl, 
                       N=N,
                       z=z,
                       p_B=p_B,
                       p_A=p_A,
                       p_I=p_I,
                       p_E=p_E + p_m,
                       p_V=p_V)[Prslt]
        
    assert N > 0 # because if N == 0 then the attack is over 
    return (
                (
                   ( 
                      ( P  *  p_m ) / ( p_B + p_E + p_m ) 
                   ) * ( epsilon + 1 )
                ) + ( n_fork * P * ( epsilon + 1 ) ) + ( n_main * ( 1 - P))
            )
#extract to .py file later:
def EV_fork_markov_blocks_effort(p_B=0,
                                 p_E=0,
                                 p_A=None,
                                 p_V=0,
                                 p_I=0,
                                 p_m=0,
                                 z=1, # This is the current state in terms of blocks behind, could also be negative i.e., represent a lead
                                 epsilon=0,
                                 n_fork=0,
                                 kr=1,
                                 kl=6,
                                 N=6, # Remaining duration of attack 
                                 Prslt="Pwin",
                                 duration=None, # this is the remaining duration, often this will just be N
                                 n_main=0): 
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m )
    if p_A < 0:
        p_A = 0 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    
    P = Pr_success(kr=kr, 
                       kl=kl, 
                       N=N,
                       z=z,
                       p_B=p_B,
                       p_A=p_A,
                       p_I=p_I,
                       p_E=p_E + p_m,
                       p_V=p_V)[Prslt]

    
    # depricated version which assume 1.0 here:
    #if duration is None: 
    #    P_block = 1.0
    
    if duration is None:
        # If not explicit duration was given
        # it is assumed the duration is N
        duration = N
    
    # Depricate version which used time here:
    #P_block = Pr_find_block_in_t_fork(p_m=p_m,
    #                                  p_fork=(p_B+p_E+p_m),
    #                                  duration=duration,
    #                                  EX=10*60)

    # Depricate version with used M blocks here:
    #P_block = Pr_find_M_blocks_in_N(M=1,
    #                                N=N,
    #                                p_m=p_m,
    #                                p_fork=(p_B+p_E+p_m))
    
    assert N > 0 # the attack is technically over
    EV = (
                (
                   ( 
                     # ( P_block * P  *  p_m ) / ( p_B + p_E + p_m ) 
                       ( P  *  p_m ) / ( p_B + p_E + p_m )
                   ) * ( 1 + epsilon )
                ) #+ ( ( P_block * (1-P) * p_m )/ (p_B + p_E + p_m) ) 
                  + ( ( (1-P) * p_m )/ (p_B + p_E + p_m) )
                  + n_fork
                  + ( n_fork * P * epsilon )
                  + ( n_main * ( 1 - P))
             )
    return EV
#extract to .py file later:
def EV_main_markov_blocks_effort(p_B=0,
                          p_E=0,
                          p_A=None,
                          p_V=0,
                          p_I=0,
                          p_m=0,
                          z=1,
                          n_main=0,
                          kr=1,
                          kl=6,
                          N=6,
                          Pm_mode="V",
                          Prslt="Pwin",
                          epsilon=0,
                          n_fork=0):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m ) 
    if p_A < 0:
        p_A = 0 
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    assert z >= 0
    assert n_main >= 0
    
    # either p_m is part of p_V, or of p_I
    assert (Pm_mode == "V" or Pm_mode == "I")
    if Pm_mode == "V":
        P = Pr_success(kr=kr, 
                       kl=kl, 
                       N=N,
                       z=z,
                       p_B=p_B,
                       p_A=p_A, 
                       p_I=p_I,
                       p_E=p_E,
                       p_V=p_V + p_m)[Prslt] # here it is part of p_V
    elif Pm_mode == "I":
        P = Pr_success(kr=kr, 
                       kl=kl, 
                       N=N,
                       z=z,
                       p_B=p_B,
                       p_A=p_A, 
                       p_I=p_I + p_m, # here it is part of p_I
                       p_E=p_E,
                       p_V=p_V)[Prslt] 

    assert N > 0 # because if N == 0 then the attack is over
    return ( 
              ( 
                  (( 1 - P ) * p_m )
                  / ( p_A + p_V + p_I + p_m ) 
              ) + ( n_main * ( 1 - P ) ) + ( n_fork ) + ( n_fork * epsilon * P  )
           )
#extract to .py file later
def e_min_markov(p_B=0,
                 p_E=0,
                 p_A=None,
                 p_V=0,
                 p_I=0,
                 p_m=0,
                 z=1,
                 kr=1,
                 kl=6,
                 N=6,
                 Prslt="Pwin"):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m )
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    
    P_join = Pr_success(  kr=kr, 
                          kl=kl, 
                          N=N,
                          z=z,
                          p_B=p_B,
                          p_A=p_A,
                          p_I=p_I,
                          p_E=p_E + p_m,
                          p_V=p_V)[Prslt]   
    
    P_abstain = Pr_success(  kr=kr, 
                             kl=kl, 
                             N=N,
                             z=z,
                             p_B=p_B,
                             p_A=p_A,
                             p_I=p_I,
                             p_E=p_E,
                             p_V=p_V + p_m)[Prslt]
    return (
             (     
                (( 1-P_abstain ) * p_m)/
                (p_A + p_V + p_I + p_m)
             )/
             (
                 (P_join * p_m)/
                 (p_B + p_E + p_m)
             )-1
           )
#extract to .py file later
def get_sympy_expr_for_e_min_markov_blocks():
    e  = Symbol('e',  real=True) # x
    Pj = Symbol('Pj', real=True) # j 
    Pa = Symbol('Pa', real=True) # a
    nm = Symbol('nm', real=True) # m
    nf = Symbol('nf', real=True) # f
    pB = Symbol('pB', real=True) # B
    pA = Symbol('pA', real=True) # A 
    pE = Symbol('pE', real=True) # E
    pV = Symbol('pV', real=True) # V
    pI = Symbol('pI', real=True) # D
    pm = Symbol('pm', real=True) # p

    e_expr = solve(( 
                    (
                      ( 
                       ((1-Pa)*pm)/(pA+pV+pI+pm)
                      ) + (nm*(1-Pa))
                        + (Pa * nf * (e+1))       # corrected, included bribe
                        - (nf*(e+1)*Pj)
                        - (nm * (1-Pj))           # corrected
                     )/
                     (
                       (Pj*pm)/(pB+pE+pm)
                      )
                    )-(e+1), e)
    return e_expr
#extract to .py file later
def get_sympy_expr_for_e_min_markov_blocks_effort():
    e  = Symbol('e', real=True)
    Pj = Symbol('Pj', real=True) # Pr join 
    Pa = Symbol('Pa', real=True) # Pr abstain
    nm = Symbol('nm', real=True)
    nf = Symbol('nf', real=True)
    pB = Symbol('pB', real=True)
    pA = Symbol('pA', real=True)
    pE = Symbol('pE', real=True)
    pV = Symbol('pV', real=True)
    pI = Symbol('pI', real=True)
    pm = Symbol('pm', real=True)
    e_expr = solve(( 
                    (
                     ( 
                      ((1-Pa)*pm)/(pA+pV+pI+pm)
                     ) + (nm*(1-Pa))
                       - (nf*e*Pj)
                       - (nf)
                       - (((1-Pj)*pm)/(pB+pE+pm))
                       + (nf)                        # corrected
                       - (nm*(1-Pj))                 # corrected
                       + (nf * e * Pa)               # corrected, added bribe
                     )/
                     (
                      (Pj*pm)/(pB+pE+pm)
                     )
                    )-(1+e), e)
    return e_expr
#extract to .py file later:
def e_min_markov_sympy(p_B=0,
                       p_E=0,
                       p_A=None,
                       p_V=0,
                       p_I=0,
                       p_m=0,
                       z=1,
                       n_main=0,
                       n_fork=0,
                       kr=1,
                       kl=6,
                       N=6,
                       Prslt="Pwin",
                       expr=None,
                       dbg=False):
    if p_A is None:
        p_A = 1 - ( p_B + p_E + p_V + p_I + p_m )
    assert math.isclose(p_B + 
                        p_E + 
                        p_A + 
                        p_V + 
                        p_I +
                        p_m, 1)
    
    P_join = Pr_success(   kr=kr, 
                           kl=kl, 
                           N=N,
                           z=z,
                           p_B=p_B,
                           p_A=p_A,
                           p_I=p_I,
                           p_E=p_E + p_m,
                           p_V=p_V)[Prslt]   
    if dbg:
        print("Pr_join(success) = ",P_join)
    
    P_abstain = Pr_success(    kr=kr, 
                               kl=kl, 
                               N=N,
                               z=z,
                               p_B=p_B,
                               p_A=p_A,
                               p_I=p_I,
                               p_E=p_E,
                               p_V=p_V + p_m)[Prslt]
    if dbg:
        print("Pr_abstain(success) = ",P_abstain)
    
    e  = Symbol('e', real=True)
    Pj = Symbol('Pj', real=True)  
    Pa = Symbol('Pa', real=True)
    nm = Symbol('nm', real=True)
    nf = Symbol('nf', real=True)
    pB = Symbol('pB', real=True)
    pA = Symbol('pA', real=True)
    pE = Symbol('pE', real=True)
    pV = Symbol('pV', real=True)
    pI = Symbol('pI', real=True)
    pm = Symbol('pm', real=True)
    
    if dbg:
        print(expr)
    
    expr_0 = expr.subs(Pj,P_join).subs(Pa,P_abstain)
    expr_1 = expr_0.subs(nm,n_main).subs(nf,n_fork)
    expr_2 = expr_1.subs(pB,p_B).subs(pA,p_A)
    expr_3 = expr_2.subs(pE,p_E).subs(pV,p_V).subs(pI,p_I).subs(pm,p_m)
    
    if dbg:
        print(expr_3)
    if expr_3.has(zoo): # expr.has(oo, -oo, zoo, nan)
        e = float("inf")
    else:
        e = float(expr_3) # Return python float
    return e
