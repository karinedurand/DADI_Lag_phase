#!/usr/bin/env python
# -*- coding: utf-8 -*- 

"""
module modeledemo with demographic models of divergence
"""

import numpy
import dadi



def AM_migr(params, (n1,n2), pts):
    nu1, nu2, migr, Tam, Ts = params

    """
    Model with split, ancient migration

    nu1: Size of population invasive after split.
    nu2: Size of population native after split.
    m12: Migration (2*Na*m).
    Ts: The scale time between the ancient migration and the present.
    Tam: The scaled time between the split and the end of migration
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate m
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=migr, m21=migr)
    # We keep the population sizes after the split to nu1 and nu2 and set the migration rates to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0,m21=0)

    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

    
def AM_bottle_migr(params, (n1,n2), pts):
    nu1, nu2,b, migr, Tam, Ts = params

    """
    Model with split, ancient migration

    nu1: Size of population invasive after split.
    nu2: Size of population native after split.
    b:  the scaling factor for the reduction of the population size (nu1) during the bottleneck event
    m: Migration  (2*Na*m).
    Ts: The scale time between the ancient migration and the present.
    Tam: The scaled time between the split and the end of migration
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    # Ensure x is always greater than 1
    b = max(b, 1.01)  # This ensures that x is greater than 1
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2,  m12=migr, m21=migr)     
    # We set the invasive population botteleneck after the split < to nu1 and the migration rates to zero
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1/b, nu2, m12=0,m21=0)         
    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs   
    
def AM_bottleGrowth_migr(params, (n1,n2), pts):
    nu1, nu2, b, g1, migr, Tam,Tb, Ts = params

    """
    Model with split, ancient migration

    nu1: Size of population invasive after split.
    nu2: Size of population native after split.
    b:  the scaling factor for the reduction of the population size (nu1) during the bottleneck event
    g1: Population growth coefficient of invasive population after the bottleneck .
    m: Migration from pop 2 to pop 1 (2*Na*m).
    Tam: The scaled time between the split and the end of migration
    Tb: The scale time between the ancient migration and the bottleneck
    Ts: The scale time between the bottleneck and the present.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    # Ensure x is always greater than 1
    b = max(b, 1.01)  # This ensures that b is greater than 1
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and the migration rate to m
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2,  m12=migr, m21=migr)  
    # We set the invasive population botteleneck after the split nu1/b and the migration rates to zero
    phi = dadi.Integration.two_pops(phi, xx, Tb, nu1/b , nu2,  m12=0, m21=0)
    # We continue the population size change after bottleneck (until present)  in invasive population and keep the migration rates to zero with b1>1
    # Ensure x is always greater than 1
    g1 = max(g1, 1.01)  # This ensures that b1 is greater than 1
    gnu1_func = lambda t: nu1 * g1**(t/Ts)
    phi = dadi.Integration.two_pops(phi, xx, Ts, gnu1_func, nu2, m12=0,m21=0)
    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs   

def IM_migr(params, (n1,n2), pts):
    nu1, nu2, migr, Ts = params
    """
    Model with migration during the divergence.

    nu1: Size of population invasive after split.
    nu2: Size of population native after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the present (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=migr, m21=migr)
    # Finally, calculate the spectrum.
    # oriented
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    ### Sum the two spectra in proportion O
    return fs
    
def SI(params, (n1,n2), pts):
    nu1, nu2, Ts = params
    """
    Model with split and complete isolation.

    nu1: Size of population invasive after split.
    nu2: Size of population native after split.
    Ts: The scaled time between the split  (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # Finally, calculate the spectrum.
    ### Sum the two spectra in proportion O
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs

