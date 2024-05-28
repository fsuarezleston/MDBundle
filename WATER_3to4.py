#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 09:39:46 2024

@author: fabs

Fabián Suárez Lestón
Universidade de Santiago de Compostela
ORCID 0000-0003-0843-667X
https://github.com/fsuarezleston

"""

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# MODULES
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

import argparse
import MDAnalysis
import numpy as np

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ARGUMENTS
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

parser = argparse.ArgumentParser(
    prog = "WATER 3to4", 
    description = "Change the water model in a system from 3 point to 4 point",
    epilog = "NOTE: It is assumed that water molecules are written in the following order: \n OW - H1 - H2 (-VS)")

parser.add_argument( "-f", "--file", type= str, default = 'system.gro',
    help= ''' File with the system.\n
    Default: %(default)s ''' )

parser.add_argument( "-wr", "--water_resname", type = str, default = "SOL",
    help = ''' Name of the current water model residue.\n
    Default: %(default)s ''' )

parser.add_argument( "-nr", "--new_resname", type = str, default = "SOL",
    help = ''' Name of the new residue.\n
    Default: %(default)s ''' )
parser.add_argument( "-nn", "--new_names", type = str, nargs='+', default = [ "OH2", "H1", "H2", "MW" ],
    help = ''' Name of the atoms in the new model.\n
    Default: %(default)s ''' )
parser.add_argument( "-vd", "--vsite_distance", type = float, default = 0.1594, # Value for the OPC4 water model
    help = ''' Disance between the O atom and the vistual site.\n
    A positive value corresponds to the direction of the sum of OH vectors.
    Default: %(default)s Å''' )

parser.add_argument( "-o", "--output", type = str, default="system_out.gro",
    help= ''' System with the new water model.\n
    Default: %(default)s ''' )

args = parser.parse_args()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# CODE
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def DefVS( water ) -> np.ndarray:
    '''
    Get the position for a virtual site 

    Parameters
    ----------
    water : MDAnalysis.core.groups.AtomGroup
        A water molecule.

    Returns
    -------
    numpy.ndarray
        The position of the virtual site.

    '''
    OW, H1, H2 = water.atoms
    
    vector = (H1.position - OW.position) + (H2.position - OW.position)
    vector /= np.linalg.norm( vector )
    
    return vector * args.vsite_distance + OW.position
    

if __name__ == "__main__":
    
    # The original system
    u = MDAnalysis.Universe( args.file )
    
    # All atoms but the water
    NON_WATER = u.select_atoms( f"not resname {args.water_resname}" )
    
    # Water residues
    WATER = u.select_atoms( f"resname {args.water_resname}" ).residues
    
    OPC_pos = []
    for mol in WATER:
        OPC_pos.append( np.vstack( [ mol.atoms.positions, DefVS( mol ) ] ) )
    OPC_pos = np.vstack( OPC_pos )
    
        
    # Generate a new universe
    NAtoms = len(NON_WATER.atoms) + len(OPC_pos)
    NResidues = len(u.residues)
    ResIndices = np.hstack( [ NON_WATER.resindices, [ i for i in np.unique( WATER.atoms.resindices ) for _ in args.opc_names ] ] )
    SegIndices = [0] * NResidues
    
    U = MDAnalysis.Universe.empty( NAtoms,
                                   n_residues = NResidues,
                                   atom_resindex = ResIndices,
                                   residue_segindex = SegIndices,
                                   trajectory = True )
        
    # Add atributes
    U.add_TopologyAttr('name', np.hstack([ NON_WATER.names, args.opc_names*len(WATER) ]) )
    U.add_TopologyAttr('type', np.hstack([ NON_WATER.types, ["O","H","H",""]*len(WATER) ]) )
    U.add_TopologyAttr('resname', np.hstack([ NON_WATER.residues.resnames, [args.opc_resname]*len(WATER) ]) )
    U.add_TopologyAttr('resids', list(range(1,NResidues+1)) )
    U.atoms.positions = np.vstack( [ NON_WATER.positions, OPC_pos ] )
    U.dimensions = u.dimensions    
    
    # Save the system
    U.atoms.write( args.output )