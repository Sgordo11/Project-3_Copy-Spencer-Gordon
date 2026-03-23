#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:34:19 2021

@author: kendrick
"""

import numpy as np

# compute unknown displacements 
def ComputeDisplacements(K, F, n_unknowns):
    # extract submatrix of unknowns
    K11 = K[0:n_unknowns,0:n_unknowns]
    F1 = F[0:n_unknowns]
    
    d = np.linalg.solve(K11,F1)
    
    return d

# postprocess the forces at known displacement nodes
def PostprocessReactions(K, d, F, n_unknowns, nodes):
    # These are computed net forces and do not
    # take into account external loads applied
    # at these nodes
    F = np.matmul(K[n_unknowns:,0:n_unknowns], d)
    
    # Postprocess the reactions
    for node in nodes:
        if node.xidx >= n_unknowns:
            node.AddReactionXForce(F[node.xidx-n_unknowns][0] - node.xforce_external)
        if node.yidx >= n_unknowns:
            node.AddReactionYForce(F[node.yidx-n_unknowns][0] - node.yforce_external)
        
    return F

# determine internal member loads
def ComputeMemberForces(bars):
    for bar in bars:
        A = bar.A
        E = bar.E
        L = bar.Length()
        lambdax,lambday = bar.LambdaTerms()
        Nx = bar.init_node.xdisp
        Ny = bar.init_node.ydisp
        Fx = bar.end_node.xdisp
        Fy = bar.end_node.ydisp

        loads = (A*E/L)*(-lambdax*Nx-lambday*Ny+lambdax*Fx+lambday*Fy)
        bar.axial_load = loads


    
# compute the normal stresses
def ComputeNormalStresses(bars):
    for bar in bars:
        F = bar.axial_load
        A = bar.A

        bar.normal_stress = F/A

# compute the critical buckling load of a member
def ComputeBucklingLoad(bars):
    for bar in bars:
        E = bar.E
        I = bar.It
        L = bar.Length()

        bar.buckling_load = (np.pi**2*E*I)/(L**2)

