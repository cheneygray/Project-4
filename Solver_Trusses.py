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
    # COMPLETE THIS FUNCTION
    # Compute member forces for all bars using equation 14-23 
    for bar in bars:
        E=bar.E
        A=bar.A
        L=bar.Length()
        c=bar.LambdaTerms()[0]
        s=bar.LambdaTerms()[1]
        
        n1=bar.init_node
        n2=bar.end_node
        
        u1=n1.xdisp
        u2=n2.xdisp
        v1=n1.ydisp
        v2=n2.ydisp
        
        u=np.array([u1,v1,u2,v2])
        T=np.array([-c,-s,c,s])
        P=(A*E/L)*np.dot(T,u)
        
        bar.axial_force=P
    return bars
    
# compute the normal stresses
def ComputeNormalStresses(bars):
    # COMPLETE THIS FUNCTION
    # Compute normal stress for all bars
    for bar in bars:
        sigma = bar.axial_force/bar.A
        bar.normal_stress = sigma
        bar.is_comuted = True
    pass

# compute the critical buckling load of a member
def ComputeBucklingLoad(bars):
    # COMPLETE THIS FUNCTION
    # Compute critical buckling load for all bars
    for bar in bars:
        E=bar.E
        I=bar.Iu
        L=bar.Length()*12
        K=1
        Pcr=(np.pi**2*E*I)/((K*L)**2)
        bar.buckling_load=Pcr
    return bars
