/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2017  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "lhrs.hpp"

#include <fstream>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/format.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/variant.hpp>

#include <core/utils/bitset_utils.hpp>
#include <core/utils/conversion_utils.hpp>
#include <core/utils/graph_utils.hpp>
#include <core/utils/range_utils.hpp>
#include <core/utils/temporary_filename.hpp>
#include <core/utils/terminal.hpp>
#include <core/utils/timer.hpp>
#include <classical/abc/gia/gia.hpp>
#include <classical/abc/gia/gia_utils.hpp>
#include <classical/abc/utils/abc_run_command.hpp>
#include <classical/functions/spectral_canonization.hpp>
#include <classical/io/read_blif.hpp>
#include <classical/optimization/exorcism_minimization.hpp>
#include <classical/utils/truth_table_utils.hpp>
#include <reversible/gate.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/functions/add_circuit.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <reversible/functions/circuit_from_string.hpp>
#include <reversible/functions/clear_circuit.hpp>
#include <reversible/io/print_circuit.hpp>
#include <reversible/synthesis/lhrs/legacy/ordering_lut_graph.hpp>
#include <reversible/synthesis/lhrs/legacy/stg_map_esop.hpp>
#include <reversible/synthesis/lhrs/legacy/stg_map_precomp.hpp>
#include <reversible/synthesis/lhrs/legacy/stg_partners.hpp>
#include <reversible/utils/circuit_utils.hpp>
#include <reversible/utils/costs.hpp>

namespace cirkit
{

namespace legacy
{

/******************************************************************************
 * Merge properties                                                           *
 ******************************************************************************/

properties::ptr merge_properties( const properties::ptr& p1, const properties::ptr& p2 )
{
  const auto p = std::make_shared<properties>();
  if ( p1 )
  {
    for ( const auto& kv : *p1 )
    {
      set( p, kv.first, kv.second );
    }
  }
  if ( p2 )
  {
    for ( const auto& kv : *p2 )
    {
      set( p, kv.first, kv.second );
    }
  }
  return p;
}

/******************************************************************************
 * Order heuristics                                                           *
 ******************************************************************************/

class lut_order_heuristic
{
public:
  enum step_type { pi, po, inv_po, zero_po, one_po, compute, uncompute };

  /* describes a single computation step */
  struct step
  {
    int                   node;          /* the node to synthesize */
    unsigned              target;        /* the target line for the result */
    step_type             type;          /* which step to perform */
    std::vector<unsigned> clean_ancilla; /* number of clean ancillae */
    std::vector<unsigned> line_map;      /* the mapping of lines to qubits */
  };
  unsigned qubits;
  unsigned lone_qubits;
  unsigned step_count;
  using step_vec = std::vector<step>;

public:
  explicit lut_order_heuristic( const gia_graph& gia, unsigned additional_ancilla = 0u )
    : _gia( gia ), _additional_ancilla( additional_ancilla )
  {
  }
  void set_pebble_stats(unsigned q, unsigned l, unsigned s)
  {
      qubits = q;
      lone_qubits = l;
      step_count = s;
  }
  unsigned get_qubits() const { return qubits;}
  unsigned get_lone_qubits() const { return lone_qubits;}
  unsigned get_step_count() const { return step_count;}

  virtual unsigned compute_steps() = 0;
  inline const step_vec& steps() const { return _steps; }

  inline unsigned& node_to_line( int index ) { return _node_to_line[index]; }
  inline unsigned node_to_line( int index ) const { return _node_to_line.find( index )->second; }
  inline unsigned& operator[]( int index ) { return _node_to_line[index]; }

  std::vector<unsigned> compute_line_map( int index ) const
  {
    std::vector<unsigned> line_map;

    _gia.foreach_lut_fanin( index, [this, &line_map]( int fanin ) {
        const auto it = _node_to_line.find( fanin );
        if ( it == _node_to_line.end() )
        {
          std::cout << "no line for node " << fanin << std::endl;
          assert( false );
        }
        line_map.push_back( node_to_line( fanin ) );
      } );

    const auto it = _node_to_line.find( index );
    if ( it == _node_to_line.end() )
    {
      std::cout << "no line for node " << index << std::endl;
      assert( false );
    }
    line_map.push_back( node_to_line( index ) );
    return line_map;
  }

  unsigned num_clean_ancilla()
  {
    return _constants.size();
  }
  void free_steps()
  {
      std::vector<step>().swap(_steps);
  }
protected:
  void add_default_input_steps()
  {
    _gia.foreach_input( [this]( int index, int e ) {
        const auto line = _next_free++;
        node_to_line( index ) = line;
        add_step( index, line, step_type::pi );
      } );
  }

  void add_default_output_steps()
  {
    _gia.foreach_output( [this]( int index, int e ) {
        const auto driver = abc::Gia_ObjFaninLit0p( _gia, abc::Gia_ManCo( _gia, e ) );
        if ( abc::Abc_Lit2Var( driver ) == 0 )
        {
          add_step( index, 0, abc::Abc_LitIsCompl( driver ) ? step_type::one_po : step_type::zero_po );
        }
        else
        {
          add_step( index, node_to_line( abc::Abc_Lit2Var( driver ) ), abc::Abc_LitIsCompl( driver ) ? step_type::inv_po : step_type::po );
        }
      } );
  }

  void add_step( int index, unsigned target, step_type type )
  {
    std::cout << _steps.capacity() << " " << _steps.size() << " " << type << "\n";
    if ( !_dry_run )
    {
      if(type == step_type::compute || type == step_type :: uncompute)
      {
        _steps.push_back( {index, target, type, _constants, compute_line_map(index)} );
      }
      else
      {
        _steps.push_back( {index, target, type, _constants, std::vector<unsigned>()} );
      }
    }
  }
 
  unsigned request_constant()
  {
    if ( !_constants.empty() )
    {
      const auto line = _constants.back();
      _constants.pop_back();
      return line;
    }

    return _next_free++;
  }

  void add_constants( unsigned max )
  {
    while ( _next_free < max )
    {
      _constants.insert( _constants.begin(), _next_free++ );
    }
  }

  void free_constant( unsigned line )
  {
    _constants.push_back( line );
  }

  inline const gia_graph& gia() const { return _gia; }
  inline unsigned additional_ancilla() const { return _additional_ancilla; }

  inline unsigned next_free() const { return _next_free; }

  void set_mem_point()
  {
    _constants_mem = _constants;
    _next_free_mem = _next_free;
  }

  void return_to_mem_point()
  {
    std::swap( _constants, _constants_mem );
    std::swap( _next_free, _next_free_mem );
  }

  void set_dry_run( bool dry_run )
  {
    _dry_run = dry_run;
  }

private:
  const gia_graph& _gia;
  unsigned _additional_ancilla;

  step_vec _steps;
  std::unordered_map<int, unsigned> _node_to_line;
  std::vector<unsigned> _constants;

  /* for memory */
  std::vector<unsigned> _constants_mem;
  unsigned _next_free_mem;

  bool _dry_run = false;

protected:
  unsigned _next_free = 0u;
};

std::ostream& operator<<( std::ostream& os, const lut_order_heuristic::step& s )
{
  switch ( s.type )
  {
  case lut_order_heuristic::pi:        os << "PI"; break;
  case lut_order_heuristic::po:        os << "PO"; break;
  case lut_order_heuristic::inv_po:    os << "PO'"; break;
  case lut_order_heuristic::zero_po:   os << "ZERO"; break;
  case lut_order_heuristic::one_po:    os << "ONE"; break;
  case lut_order_heuristic::compute:   os << "COMPUTE"; break;
  case lut_order_heuristic::uncompute: os << "UNCOMPUTE"; break;
  }
  os << boost::format( " %d ↦ %d" ) % s.node % s.target;
  return os;
}

class defer_lut_order_heuristic : public lut_order_heuristic
{
public:
  defer_lut_order_heuristic( const gia_graph& gia, unsigned additional_ancilla )
    : lut_order_heuristic( gia, additional_ancilla )
  {
    partners = find_stg_partners( gia );
  }

public:
  virtual unsigned compute_steps()
  {
    set_mem_point();
    set_dry_run( true );
    const auto next_free = compute_steps_int();
    set_dry_run( false );
    return_to_mem_point();

    return compute_steps_int( next_free + additional_ancilla() );
  }

private:
  unsigned compute_steps_int( unsigned add_frees = 0u )
  {
    gia().init_lut_refs();

    add_default_input_steps();

    if ( add_frees )
    {
      add_constants( add_frees );
    }

    adjust_indegrees();

    gia().foreach_lut( [this]( int index ) {
        const auto is_partners = partners.has_partners( index );

        if ( !is_partners || partners.partners( index ).front() == index )
        {
          const auto target = request_constant();
          (*this)[index] = target;

          if ( is_partners )
          {
            for ( auto p : partners.partners( index ) )
            {
              add_step( p, target, lut_order_heuristic::compute );
            }
          }
          else
          {
            add_step( index, target, lut_order_heuristic::compute );
          }

          /* start uncomputing */
          if ( gia().lut_ref_num( index ) == 0 )
          {
            visited.clear();
            decrease_children_indegrees( index );
            uncompute_children( index );
          }
        }
      } );

    add_default_output_steps();

    return next_free();
  }

  void adjust_indegrees()
  {
    gia().foreach_output( [this]( int index, int e ) {
        const auto driver = abc::Gia_ObjFaninId0p( gia(), abc::Gia_ManCo( gia(), e ) );
        output_luts.push_back( driver );
        gia().lut_ref_dec( driver );
      } );
  }

  void decrease_children_indegrees( int index )
  {
    gia().foreach_lut_fanin( index, [this]( int fanin ) {
        if ( gia().is_lut( fanin ) )
        {
          gia().lut_ref_dec( fanin );
        }
      } );
  }

  void uncompute_children( int index )
  {
    gia().foreach_lut_fanin( index, [this]( int fanin ) {
        if ( gia().is_lut( fanin ) && gia().lut_ref_num( fanin ) == 0 )
        {
          uncompute_node( fanin );
        }
      } );
  }

  void uncompute_node( int index )
  {
    if ( is_visited( index ) ) return;
    assert( gia().lut_ref_num( index ) == 0 );

    if ( !is_output_lut( index ) )
    {
      const auto target = (*this)[index];

      const auto is_partners = partners.has_partners( index );
      if ( !is_partners || partners.partners( index ).front() == index )
      {
        if ( is_partners )
        {
          for ( auto p : partners.partners( index ) )
          {
            add_step( p, target, lut_order_heuristic::uncompute );
          }
        }
        else
        {
          add_step( index, target, lut_order_heuristic::uncompute );
        }
        free_constant( target );
      }
    }

    visited.push_back( index );

    decrease_children_indegrees( index );
    uncompute_children( index );
  }

  void print_lut_refs()
  {
    std::cout << "[i] LUT refs:";
    gia().foreach_lut( [this]( int index ) {
        std::cout << boost::format( "  %d:%d" ) % index % gia().lut_ref_num( index );
      });
    std::cout << std::endl;
  }

private:
  inline bool is_visited( int index ) const
  {
    return std::find( visited.begin(), visited.end(), index ) != visited.end();
  }

  inline bool is_output_lut( int index ) const
  {
    return std::find( output_luts.begin(), output_luts.end(), index ) != output_luts.end();
  }

  std::vector<int> visited;
  std::vector<int> output_luts;

  stg_partners partners;
};


class rev_pebble_heuristic : public lut_order_heuristic
{
public:
  rev_pebble_heuristic( const gia_graph& gia, unsigned additional_ancilla )
    : lut_order_heuristic( gia, additional_ancilla )
  {
    partners = find_stg_partners( gia );
  }

public:
  virtual unsigned compute_steps()
  {
    /* pebble game begins */
    set_mem_point();
    set_dry_run( false );
    gia().init_lut_refs();
    std::cout << "Copying gia to Boost GML\n";
    std::vector<std::pair<int,int>> edges;
    gia().foreach_output( [this]( int index, int e ) {                                                                                              
        auto const driver = abc::Gia_ObjFaninId0p( gia(), abc::Gia_ManCo( gia(), e ) );
        if(gia().is_lut(driver) && driver != 0)
	    	output_LUT_set.insert( driver );
        // gia().lut_ref_dec( driver );
        // Remove after testing
        #ifdef ORDERING_DEBUG

        std::cout <<"X" << index << ","<< e << " driver:" << driver << "[" ;
        std::cout << gia().is_lut(driver) << "]\n";
        gia().foreach_lut_fanin( driver, [driver,this]( int fanin ) {
            std::cout <<  (int) driver << "<-" << fanin;
            if ( gia().is_lut( fanin ) )
            {
                gia().lut_ref_dec( fanin );
                std::cout <<"LUT";
            }
            else
                std::cout << "PI";
            std::cout << "\n";
        } );

        #endif
    } );

       std::set <int> addedLut;
       std::set <int> connectedLut; 
       std::vector <int> toProcessLut;

       for(auto lut:output_LUT_set)
       {
           toProcessLut.push_back(lut);
       }

       std::vector<int> zeroFanInLut;
       bool isZeroPresent = false;
       while(!toProcessLut.empty())
       {
           int lut = toProcessLut.back();
           toProcessLut.pop_back();
          
    
           if(addedLut.find(lut) != addedLut.end())
               continue;
           addedLut.insert(lut);  
           bool zeroFanIn = true; 
           if(lut == 0){
               std::cout << "Node zero found \n";
               isZeroPresent = true;
               // continue;
           } 

           gia().foreach_lut_fanin( lut, [&toProcessLut,&connectedLut, &edges,&zeroFanIn, addedLut,lut,this]( int fanin ) {
               
		        if(gia().is_lut(fanin) && fanin != 0)                
                {
                    edges.push_back(std::make_pair(fanin,lut));
                    #ifdef ORDERING_DEBUG
                    std::cout <<  (int) lut << "<-" << fanin << "[LUT]\n";
                    #endif
               
                    zeroFanIn = false;
                    connectedLut.insert(fanin);
                    connectedLut.insert(lut);
                
                    if (addedLut.find(fanin) == addedLut.end()) 
               	    { 
                        toProcessLut.push_back(fanin);
                    }
                }
                else
                {
                     #ifdef ORDERING_DEBUG
                    std::cout <<  (int) lut << "<-" << fanin << "[PI]\n";
                    #endif
               

                }
            } );
            
       }
    for(auto v: addedLut)
    {
        if(connectedLut.find(v) == connectedLut.end())
            zeroFanInLut.push_back(v);
    }
    std::vector<std::pair<int,int>> orderedEdges;
    std::map<int,int> vertexMap;
    int vertexLabel = 0;
    int source, dest;
       std::cout << "\n";
       for(auto edge : edges)
       {
           source = (int) edge.first;
           dest = (int) edge.second;
                  if(vertexMap.find(source) == vertexMap.end())
           {
               
               vertexMap[source] = vertexLabel;
               //std::cout << vertexLabel << "new assigned\n";
               vertexLabel++;
           }
           if(vertexMap.find(dest) == vertexMap.end())
           {
               vertexMap[dest] = vertexLabel;
               //std::cout << vertexLabel << "new assigned\n";
               vertexLabel++;
           }
           source = vertexMap[source];
           dest = vertexMap[dest];
           #ifdef ORDERING_DEBUG

           std::cout << edge.first << "[" << source << "] ->";
           std::cout << edge.second << "[" << dest << "]\n";
           #endif

    
           orderedEdges.push_back(std::make_pair(source,dest));
       }

    std::map<int,int> newVertexMap;
    //std::cout << "vertex map:\n";
    for(auto edge: vertexMap)
    {
        #ifdef ORDERING_DEBUG 
        std::cout << "vmap_f" << edge.first << "->" << edge.second << "\n";
        #endif 
        newVertexMap[edge.second] = edge.first; 

    }  

    std::vector<int> bglOutput;
    for(auto out: output_LUT_set)
    {
        if(vertexMap.find(out) == vertexMap.end())
        {
            vertexMap[out] = vertexLabel;
            std::cout << "new output vertex not present in BGL" << out << "\n"; 
            vertexLabel++;
        }
        else
            bglOutput.push_back(vertexMap[out]);
    }
    unsigned int reqQbits = 0;
    int available=0; 
    int loneLutQubits = 0;
    LutGraph *lG;
    lG = new LutGraph(orderedEdges,bglOutput, output_LUT_set.size());
 
    if(!edges.empty())
    {
        lG->pebble();
        reqQbits = lG->getQubitCount();    
        available = lG->getQubitsFree(); //starting qubits - no. of outputs
    }
    //delete lG;
    std::set<int> lutMapped;
    int stepCount = 0;
    for(auto const step : lG->getPebbleSteps())
    {
        int index = newVertexMap[std::get<1>(step)];
        lutMapped.insert(index); 
        stepCount++;
    }
    loneLutQubits = zeroFanInLut.size();
    
    //std::cout << "lone : " << loneLutQubits << "\n";
    loneLutQubits = available > loneLutQubits ? 0: loneLutQubits - available;
    std::ofstream debugStats ("debugStats.txt", std::ios::out|std::ios::app);
    debugStats << "req:" << reqQbits << " lone:" << loneLutQubits;
    debugStats << " total:" << reqQbits+loneLutQubits;
    debugStats << " steps:"<< stepCount <<"\n";
    debugStats.close();
    set_pebble_stats(reqQbits+loneLutQubits, loneLutQubits, stepCount);
    // NOTE : code for adding steps 
    /*
    add_default_input_steps();
    //  add the additional lines needed 
    add_constants(reqQbits+loneLutQubits+additional_ancilla());
    
    std::map<int,int> targetMap;
    std::map<int,int> my_qubit_their_qubit_map;

    int max_qubit = -1;
    for(auto lut : zeroFanInLut)
    {
       int target =  request_constant();
       max_qubit = max_qubit < target ? target : max_qubit;

        (*this)[lut] = target;
        add_step(lut, target, lut_order_heuristic::compute);   
        targetMap[lut] = target;
        #ifdef ORDERING_DEBUG
            std::cout << "zero fan in steps" << lut << "\n";
        #endif
    }
    for(uint i=1;i<= reqQbits;i++)
     {
     	my_qubit_their_qubit_map[i] = request_constant(); //my ith qbit maps to req_const()
     	max_qubit = max_qubit < my_qubit_their_qubit_map[i] ? my_qubit_their_qubit_map[i] : max_qubit;
        #ifdef ORDERING_DEBUG 
        std::cout << "virtual q:" << i << " mapped:" << my_qubit_their_qubit_map[i] << " "; 
        #endif 
     }
    for(auto const step : lG->getPebbleSteps())
    {
        int index = newVertexMap[std::get<1>(step)];
        
        if(index == 0)
        {
            std::cout << " Constant node 0 : not added to steps\n";
            continue;
        }
        
        //#ifdef ORDERING_DEBUG
            std::cout << "step : " << std::get<1>(step) << ":";
            std::cout << index  << "->" << std::get<2>(step) << " " << std::get<0>(step) << "\n " ;
        //#endif 

        if(std::get<0>(step) == 0)
        {
            int target =  my_qubit_their_qubit_map[std::get<2>(step)];//target is a new qubit
            (*this)[index] = target;
            add_step(index, target, lut_order_heuristic::compute);   
            targetMap[index] = target;         
        }
        else if(std::get<0>(step) == 1)
        {
            int target = targetMap[index];
            add_step( index, target, lut_order_heuristic::uncompute );
            free_constant( target );
        }
        else
        {
        	std::cout<<"INVALID STEP TYPE\n";
        	exit(1);
        }
    }
    delete lG;
    std::cout<<"MAX_QUBITS REQD"<<max_qubit+1<<"\n";
    
    add_default_output_steps();
    std::cout << "Completed adding steps\n";
  */ 
    int neededLines = next_free();
    return_to_mem_point();
    //complete_add_step();
    return  neededLines; 
    
  }

private:
  unsigned compute_steps_int( unsigned add_frees = 0u )
  {
    gia().init_lut_refs();

    add_default_input_steps();

    if ( add_frees )
    {
      add_constants( add_frees );
    }

    adjust_indegrees();

    gia().foreach_lut( [this]( int index ) {
        const auto is_partners = partners.has_partners( index );

        if ( !is_partners || partners.partners( index ).front() == index )
        {
          const auto target = request_constant();
          (*this)[index] = target;

          if ( is_partners )
          {
            for ( auto p : partners.partners( index ) )
            {
              add_step( p, target, lut_order_heuristic::compute );
            }
          }
          else
          {
            add_step( index, target, lut_order_heuristic::compute );
          }

          if ( gia().lut_ref_num( index ) == 0 )
          {
            visited.clear();
            decrease_children_indegrees( index );
            uncompute_children( index );
          }
        }
      } );

    add_default_output_steps();

    return next_free();
  }

  void adjust_indegrees()
  {
    gia().foreach_output( [this]( int index, int e ) {
        const auto driver = abc::Gia_ObjFaninId0p( gia(), abc::Gia_ManCo( gia(), e ) );
        output_luts.push_back( driver );
        gia().lut_ref_dec( driver );
      } );
  }

  void decrease_children_indegrees( int index )
  {
    gia().foreach_lut_fanin( index, [this]( int fanin ) {
        if ( gia().is_lut( fanin ) )
        {
          gia().lut_ref_dec( fanin );
        }
      } );
  }

  void uncompute_children( int index )
  {
    gia().foreach_lut_fanin( index, [this]( int fanin ) {
        if ( gia().is_lut( fanin ) && gia().lut_ref_num( fanin ) == 0 )
        {
          uncompute_node( fanin );
        }
      } );
  }

  void uncompute_node( int index )
  {
    if ( is_visited( index ) ) return;
    assert( gia().lut_ref_num( index ) == 0 );

    if ( !is_output_lut( index ) )
    {
      const auto target = (*this)[index];

      const auto is_partners = partners.has_partners( index );
      if ( !is_partners || partners.partners( index ).front() == index )
      {
        if ( is_partners )
        {
          for ( auto p : partners.partners( index ) )
          {
            add_step( p, target, lut_order_heuristic::uncompute );
          }
        }
        else
        {
          add_step( index, target, lut_order_heuristic::uncompute );
        }
        free_constant( target );
      }
    }

    visited.push_back( index );

    decrease_children_indegrees( index );
    uncompute_children( index );
  }

  void print_lut_refs()
  {
    std::cout << "[i] LUT refs:";
    gia().foreach_lut( [this]( int index ) {
        std::cout << boost::format( "  %d:%d" ) % index % gia().lut_ref_num( index );
      });
    std::cout << std::endl;
  }

private:
  inline bool is_visited( int index ) const
  {
    return std::find( visited.begin(), visited.end(), index ) != visited.end();
  }

  inline bool is_output_lut( int index ) const
  {
    return std::find( output_luts.begin(), output_luts.end(), index ) != output_luts.end();
  }

  std::vector<int> visited;
  std::vector<int> output_luts;
  std::set<int> output_LUT_set;
  stg_partners partners;
};

/******************************************************************************
 * Manager                                                                    *
 ******************************************************************************/

class lut_based_synthesis_manager
{
public:
  lut_based_synthesis_manager( circuit& circ, const gia_graph& gia, const lhrs_params& params, lhrs_stats& stats )
    : circ( circ ),
      gia( gia ),
      params( params ),
      stats( stats ),
      pbar( "[i] step %5d/%5d   dd = %5d   ld = %5d   cvr = %6.2f   esop = %6.2f   map = %6.2f   clsfy = %6.2f   total = %6.2f", params.progress )
  {
      if( params.assignment_strategy == cirkit::legacy::qubit_assignment_strategy::defer ) 
      {
        order_heuristic = std::make_shared<defer_lut_order_heuristic>( gia, params.additional_ancilla ) ;
      }
      else if ( params.assignment_strategy == cirkit::legacy::qubit_assignment_strategy::rev_pebble )
      {
        order_heuristic = std::make_shared<rev_pebble_heuristic>( gia, params.additional_ancilla );
      }

  }

  bool run()
  {
    clear_circuit( circ );

    const auto lines = order_heuristic->compute_steps();
    circ.set_lines( lines );

    std::vector<std::string> inputs( lines, "0" );
    std::vector<std::string> outputs( lines, "0" );
    std::vector<constant> constants( lines, false );
    std::vector<bool> garbage( lines, true );
    stats.qubits = order_heuristic->get_qubits();
    stats.lone_qubits = order_heuristic->get_lone_qubits();
    stats.step_count = order_heuristic->get_step_count();

    std::unordered_map<unsigned, lut_order_heuristic::step_type> orig_step_type; /* first step type of an output */
    /*
    auto step_index = 0u;
    pbar.keep_last();
    for ( const auto& step : order_heuristic->steps() )
    {
      if ( params.verbose )
      {
        std::cout << step << std::endl;
      }
      pbar( ++step_index, order_heuristic->steps().size(), stats.num_decomp_default, stats.num_decomp_lut, stats.map_esop_stats.cover_runtime, stats.map_esop_stats.exorcism_runtime, stats.map_luts_stats.mapping_runtime, stats.map_precomp_stats.class_runtime, stats.synthesis_runtime );
      increment_timer t( &stats.synthesis_runtime );

      switch ( step.type )
      {
      case lut_order_heuristic::pi:
        inputs[step.target] = outputs[step.target] = gia.input_name( abc::Gia_ManIdToCioId( gia, step.node ) );
        constants[step.target] = boost::none;
        orig_step_type[step.target] = lut_order_heuristic::po;
        break;

      case lut_order_heuristic::zero_po:
      case lut_order_heuristic::one_po:
        circ.set_lines( circ.lines() + 1 );
        inputs.push_back( step.type == lut_order_heuristic::zero_po ? "0" : "1" );
        constants.push_back( step.type == lut_order_heuristic::zero_po ? false : true );
        outputs.push_back( gia.output_name( abc::Gia_ManIdToCioId( gia, step.node ) ) );
        garbage.push_back( false );
        break;

      case lut_order_heuristic::po:
      case lut_order_heuristic::inv_po:
        if ( outputs[step.target] != "0" )
        {
          circ.set_lines( circ.lines() + 1 );
          inputs.push_back( "0" );
          constants.push_back( false );
          outputs.push_back( gia.output_name( abc::Gia_ManIdToCioId( gia, step.node ) ) );
          garbage.push_back( false );

          const auto pol = orig_step_type[step.target] == step.type ? true : false;
          if ( !params.onlylines )
          {
            append_cnot( circ, make_var( step.target, pol ), circ.lines() - 1 );
          }
        }
        else
        {
          outputs[step.target] = gia.output_name( abc::Gia_ManIdToCioId( gia, step.node ) );
          garbage[step.target] = false;

          if ( step.type == lut_order_heuristic::inv_po && !params.onlylines )
          {
            append_not( circ, step.target );
          }
          orig_step_type.insert( {step.target, step.type} );
        }
        break;

      case lut_order_heuristic::compute:
        if ( !params.onlylines )
        {
          synthesize_node( step.node, false, step.clean_ancilla, step.line_map );
        }
        break;

      case lut_order_heuristic::uncompute:
        if ( !params.onlylines )
        {
          synthesize_node( step.node, true, step.clean_ancilla , step.line_map );
        }
        break;
      }
    }
    */
    std::cout << "completed run()\n";
    /*circ.set_inputs( inputs );
    circ.set_outputs( outputs );
    circ.set_constants( constants );
    circ.set_garbage( garbage ); */
    //order_heuristic->free_steps();
    return true;
  }

private:
  boost::dynamic_bitset<> get_affected_lines( unsigned begin, unsigned end )
  {
    boost::dynamic_bitset<> mask( circ.lines() );

    while ( begin != end )
    {
      mask |= get_line_mask( circ[begin++], circ.lines() );
    }

    return mask;
  }

  inline void synthesize_node( int index, bool lookup, const std::vector<unsigned>& clean_ancilla, const std::vector<unsigned>& line_map)
  {
    /* track costs */
    const auto begin = circ.num_gates();
   // const auto line_map = order_heuristic->compute_line_map( index );

    switch ( params.mapping_strategy )
    {
    case lhrs_mapping_strategy::direct:
      synthesize_node_direct( index, lookup, line_map, clean_ancilla );
      break;
    case lhrs_mapping_strategy::lut_based_min_db:
    case lhrs_mapping_strategy::lut_based_best_fit:
      synthesize_node_lut_based( index, lookup, line_map, clean_ancilla );
      break;
    case lhrs_mapping_strategy::lut_based_pick_best:
      synthesize_node_pick_best( index, lookup, line_map, clean_ancilla );
      break;
    case lhrs_mapping_strategy::shannon:
      synthesize_node_shannon( index, lookup, line_map, clean_ancilla );
      break;
    }

    /* track costs */
    if ( params.count_costs )
    {
      const auto end = circ.num_gates();
      stats.gate_costs.push_back( costs( circ, begin, end, costs_by_gate_func( t_costs() ) ) );
      stats.line_maps.push_back( line_map );
      stats.affected_lines.push_back( get_index_vector( get_affected_lines( begin, end ) ) );
      stats.clean_ancillas.push_back( clean_ancilla );
    }
  }

  void synthesize_node_direct( int index, bool lookup, const std::vector<unsigned>& line_map, const std::vector<unsigned>& clean_ancilla )
  {
    const auto lut = gia.extract_lut( index );

    const auto sp = pbar.subprogress();
    stg_map_esop( circ, lut, line_map, params.map_esop_params, stats.map_esop_stats );
  }

  void synthesize_node_lut_based( int index, bool lookup, const std::vector<unsigned>& line_map, const std::vector<unsigned>& clean_ancilla )
  {
    const auto lut = gia.extract_lut( index );

    if ( params.max_func_size == 0u )
    {
      params.map_luts_params.max_cut_size = params.map_precomp_params.class_method == 0u ? 5 : 4;
    }
    else
    {
      params.map_luts_params.max_cut_size = params.max_func_size;
    }

    const auto sp = pbar.subprogress();
    stg_map_luts( circ, lut, line_map, clean_ancilla, params.map_luts_params, stats.map_luts_stats );
  }

  inline void append_circuit_fast( const circuit& src )
  {
    boost::push_back( circ.gates, src.gates );
  }

  inline circuit get_fast_circuit() const
  {
    circuit c;
    c._lines = circ.lines();
    return c;
  }

  void synthesize_node_pick_best( int index, bool lookup, const std::vector<unsigned>& line_map, const std::vector<unsigned>& clean_ancilla )
  {
    using candidate_t = std::pair<circuit, cost_t>;
    std::vector<candidate_t> candidates;

    const auto lut = gia.extract_lut( index );

    const auto sp = pbar.subprogress();

    const auto old_strategy = params.map_luts_params.strategy;
    for ( const auto& strategy : {stg_map_luts_params::mapping_strategy::mindb, stg_map_luts_params::mapping_strategy::bestfit} )
    {
      /* cut size 4 */
      params.map_luts_params.strategy = strategy;
      {
        auto lcirc = get_fast_circuit();
        params.map_luts_params.max_cut_size = 4;
        stg_map_luts( lcirc, lut, line_map, clean_ancilla, params.map_luts_params, stats.map_luts_stats );
        if ( lcirc.num_gates() )
        {
          candidates.push_back( {lcirc, costs( lcirc, costs_by_gate_func( t_costs() ) )} );
        }
      }

      /* cut size 5 */
      if ( params.map_precomp_params.class_method == 0u )
      {
        auto lcirc = get_fast_circuit();
        params.map_luts_params.max_cut_size = 5;
        stg_map_luts( lcirc, lut, line_map, clean_ancilla, params.map_luts_params, stats.map_luts_stats );
        if ( lcirc.num_gates() )
        {
          candidates.push_back( {lcirc, costs( lcirc, costs_by_gate_func( t_costs() ) )} );
        }
      }
    }

    params.map_luts_params.strategy = old_strategy;

    if ( !candidates.empty() )
    {
      const auto best_candidate = std::min_element( candidates.begin(), candidates.end(),
                                                    []( const candidate_t& c1, const candidate_t& c2 ) {
                                                      return c1.second < c2.second;
                                                    } );

      append_circuit_fast( best_candidate->first );
      ++stats.num_decomp_lut;
      return;
    }

    stg_map_esop( circ, lut, line_map, params.map_esop_params, stats.map_esop_stats );
  }

  void synthesize_node_shannon( int index, bool lookup, const std::vector<unsigned>& line_map, const std::vector<unsigned>& clean_ancilla )
  {
    const auto lut = gia.extract_lut( index );

    const auto sp = pbar.subprogress();
    stg_map_shannon( circ, lut, line_map, clean_ancilla, params.map_shannon_params, stats.map_shannon_stats );
  }


private:
  circuit& circ;
  const gia_graph& gia;

  const lhrs_params& params;
  lhrs_stats& stats;

  std::unordered_map<unsigned, circuit> computed_circuits;

  std::shared_ptr<lut_order_heuristic> order_heuristic;

  progress_line pbar;
};

/******************************************************************************
 * Public functions                                                           *
 ******************************************************************************/

bool lut_based_synthesis( circuit& circ, const gia_graph& gia, const lhrs_params& params, lhrs_stats& stats )
{
  /* timing */
  reference_timer t( &stats.runtime );

  lut_based_synthesis_manager mgr( circ, gia, params, stats );
  const auto result = mgr.run();

  return result;
}

}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
