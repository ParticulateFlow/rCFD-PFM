#ifndef RCFD_LAYERS
#define RCFD_LAYERS

#include "rCFD_macros.h"
#include "rCFD_globals.h"


/* (C)  2022
    Stefan Pirker
    Particulate Flow Modelling
    Johannes Kepler University, Linz, Austria
    www.particulate-flow.at
*/  
 
    void rCFD_map_parent(const short i_layer)
    {
#if RP_NODE     
        short debug_this_code = 1;
        
        if(upper_layer >= Solver_Dict.number_of_layers){
            
            Message("\nWARNING myid %d: Tried to map to non-existing parent of layer %d", i_layer);
            
            return;
        }
        
        int i_cell, i_data, i_phase;
        
        int upper_cell;
        
        if(debug_this_code){

            i_data = 0;
                
            i_phase = 0;
        
            loop_cells_of_upper_layer{              
                
                Topo.Cell[upper_layer].data[_i_data] = 0.0;
            }
            
            loop_cells{
                
                upper_cell = _C.parent_cell[i_cell];
                
                Topo.Cell[upper_layer].data[i_phase][upper_cell][i_data] += 1.0;        
            }   
        }
#endif      
    }  

    void rCFD_map_children(const short i_layer)
    {
#if RP_NODE     
        short debug_this_code = 1;
        
        if(i_layer < 1){
            
            Message("\nWARNING myid %d: Tried to map to non-existing children of layer %d", i_layer);
            
            return;
        }
        
        int i_cell, i_data, i_phase, i_child;
        
        int lower_cell;
        
        if(debug_this_code){

            i_data = 0;

            i_phase = 0;
            
            loop_cells_of_lower_layer{
                
                Topo.Cell[lower_layer].data[_i_data] = 0.0;
            }
            
            loop_cells{
                
                loop_children{
                    
                    lower_cell = _C.child_index[i_cell][i_child];
                    
                    Topo.Cell[lower_layer].data[i_phase][lower_cell][i_data] = _C.data[_i_data];
                }
            }   
        }
#endif      
    }   
    

#endif  