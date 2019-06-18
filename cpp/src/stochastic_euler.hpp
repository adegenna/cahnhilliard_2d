/*
 libs/numeric/odeint/examples/stochastic_euler.hpp

 Copyright 2012 Karsten Ahnert
 Copyright 2012 Mario Mulansky

 Stochastic euler stepper example and Ornstein-Uhlenbeck process

 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#include <vector>
#include <iostream>
#include <boost/random.hpp>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

//[ stochastic_euler_class
 //template< size_t N >
class stochastic_euler
{
public:

    typedef unsigned short order_type;

    typedef boost::numeric::odeint::stepper_tag stepper_category;

    static order_type order( void ) { return 1; }

    template< class System >
    void do_step( System system , std::vector< double > &x , double t , double dt ) const
    {
        std::vector<double> det(x.size()) , stoch(x.size()) ;
        system.first( x , det , t );
        system.second( x , stoch );
        for( size_t i=0 ; i<x.size() ; ++i )
          x[i] += dt * det[i] + sqrt( dt ) * stoch[i];
    }
};
//]

//[ stochastic_euler_ornstein_uhlenbeck_def
//const static size_t N = 1;;

struct ornstein_det
{
  void operator()( const std::vector< double > &x , std::vector< double > &dxdt ) const
    {
        dxdt[0] = -x[0];
    }
};

struct ornstein_stoch
{
    boost::mt19937 &m_rng;
    boost::normal_distribution<> m_dist;

  ornstein_stoch( boost::mt19937 &rng , double sigma ) : m_rng( rng ) , m_dist( 0.0 , sigma ) { }

  void operator()( const std::vector< double > &x , std::vector< double > &dxdt )
    {
      for (int i=0; i<dxdt.size(); i++) {
        dxdt[i] = m_dist( m_rng );
      }
    }
};
//]

struct streaming_observer
{
    template< class State >
    void operator()( const State &x , double t ) const
    {
        std::cout << t << "\t" << x[0] << "\n";
    }
};


// int main( int argc , char **argv )
// {
//     using namespace std;
//     using namespace boost::numeric::odeint;

//     //[ ornstein_uhlenbeck_main
//     boost::mt19937 rng;
//     double dt = 0.1;
//     state_type x = {{ 1.0 }};
//     integrate_const( stochastic_euler< N >() , make_pair( ornstein_det() , ornstein_stoch( rng , 1.0 ) ),
//             x , 0.0 , 10.0 , dt , streaming_observer() );
//     //]
//     return 0;
// }
