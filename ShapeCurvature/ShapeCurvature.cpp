/*
Copyright (c) 2016, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include <Include/PreProcessor.h>

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>

#include <Eigen/Sparse>
#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif // EIGEN_USE_MKL_ALL

#include <Misha/CmdLineParser.h>
#include <Misha/Miscellany.h>
#include <Misha/Algebra.h>
#include <Misha/Ply.h>
#include <Misha/PlyVertexData.h>
#include <Misha/FEM.h>
#include <Misha/MultiThreading.h>

#include <Include/CurvatureMetric.h>
#include <Include/GradientDomain.h>

using namespace MishaK;
using namespace MishaK::CmdLineParser;
using namespace MishaK::Geometry;
using namespace MishaK::Array;
using namespace MishaK::MultiThreading;

CmdLineParameter< std::string >
	In( "in" ) ,
	Out( "out" );

CmdLineParameter< float >
	NormalSmoothingWeight( "nWeight" , 1e-4f );

CmdLineReadable
	Verbose( "verbose" );

CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&NormalSmoothingWeight ,
	&Verbose ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input geometry>\n" , In.name.c_str() );
	printf( "\t[--%s <output geometry>]\n" , Out.name.c_str() );
	printf( "\t[--%s <normal smoothing weight>=%g]\n" , NormalSmoothingWeight.name.c_str() , NormalSmoothingWeight.value );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

#ifdef EIGEN_USE_MKL_ALL
using Solver = Eigen::PardisoLLT< Eigen::SparseMatrix< double , Eigen::ColMajor , __int64 > >;
#else // !EIGEN_USE_MKL_ALL
using Solver = Eigen::SimplicialLLT< Eigen::SparseMatrix< double > >;
#endif // EIGEN_USE_MKL_ALL

enum{ VERTEX_POSITION , VERTEX_NORMAL , VERTEX_COLOR };

int main( int argc , char* argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	static const unsigned int K = 2;
	static const unsigned int Dim = K+1;
	using Real = double;
	using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , Dim > , VertexFactory::NormalFactory< float , Dim > , VertexFactory::RGBColorFactory< float > >;
	using Vertex = typename Factory::VertexType;

	Miscellany::Timer timer;
	{
		std::vector< SimplexIndex< K , int > > simplices;
		std::vector< Vertex > vertices;

		int file_type;

		Miscellany::PerformanceMeter pMeter( '.' );

		//////////////////////
		// Read in the data //
		pMeter.reset();
		{
			Factory factory;
			bool *readFlags = new bool[ factory.plyReadNum() ];
			file_type = PLY::ReadSimplices( In.value , factory , vertices , simplices , readFlags );

			bool hasNormals = factory.template plyValidReadProperties< VERTEX_NORMAL >( readFlags );

			// Create normals if they are not already there
			if( !hasNormals )
			{
				for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get< VERTEX_NORMAL >() = Point< float , Dim >();
				for( int i=0 ; i<simplices.size() ; i++ )
				{
					Simplex< float , Dim , K > s;
					for( unsigned int k=0 ; k<=K ; k++ ) s[k] = vertices[ simplices[i][k] ].template get< VERTEX_POSITION >();
					Point< float , Dim > n = s.normal();
					for( int k=0 ; k<=K ; k++ ) vertices[ simplices[i][k] ].template get< VERTEX_NORMAL >() += n;
				}
				for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get< VERTEX_NORMAL >() /= Point< float , Dim >::Length( vertices[i].template get< VERTEX_NORMAL >() );
			}
			delete[] readFlags;

			if( Verbose.set ) std::cout << "Source Vertices / Simplices: " << vertices.size() << " / " << simplices.size() << std::endl;
		}
		if( Verbose.set ) std::cout << pMeter( "Read mesh" ) << std::endl;
		// Read in the data //
		//////////////////////

		//////////////////////////
		// Set the Riemannian mesh
		pMeter.reset();
		Real originalArea = 0;
		FEM::RiemannianMesh< Real > mesh( GetPointer( simplices ) , simplices.size() );
		{
			// Create the embedded metric and normalize to have unit area
			mesh.template setMetricFromEmbedding< Dim >( [&]( unsigned int i ){ return vertices[i].template get< VERTEX_POSITION >(); } );
			originalArea = mesh.area();
			mesh.makeUnitArea();
		}
		if( Verbose.set ) std::cout << pMeter( "Set Riemannian mesh" ) << std::endl;
		// Set the Riemannian mesh
		//////////////////////////

		/////////////////////
		// Smooth the normals
		if( NormalSmoothingWeight.value>0 )
		{
			Solver solver;

			///////////////////////////////////////
			// System matrix symbolic factorization
			{
				pMeter.reset();
				solver.analyzePattern( mesh.template massMatrix< FEM::BASIS_0_WHITNEY , true >() );
				if( Verbose.set ) std::cout << pMeter( "Symbolic factorization" ) << std::endl;
			}
			// System matrix symbolic factorization
			///////////////////////////////////////

			std::vector< Point< Real , Dim > > normals( vertices.size() );
			for( unsigned int i=0 ; i<normals.size() ; i++ ) normals[i] = vertices[i].template get< VERTEX_NORMAL >();

			Miscellany::PerformanceMeter pMeter( '.' );
			normals = GradientDomain::ProcessVertexVertex< Solver , Point< Real , Dim > , Real >( solver , mesh , (Real)1. , NormalSmoothingWeight.value , [&]( unsigned int v ){ return normals[v]; } , []( unsigned int ){ return Point< Real , Dim >(); } );
			ThreadPool::ParallelFor
			(
				0 , normals.size() ,
				[&]( unsigned int , size_t i ){ normals[i] /= Point< Real , Dim >::Length( normals[i] ); }
			);
			for( unsigned int i=0 ; i<normals.size() ; i++ ) vertices[i].template get< VERTEX_NORMAL >() = normals[i];
			if( Verbose.set ) std::cout << pMeter( "Normal smoothing" ) << std::endl;
		}
		// Smooth the normals
		/////////////////////

		/////////////////////////
		// Compute the curvatures
		{
			Miscellany::PerformanceMeter pMeter( '.' );

			Real scl = (Real)(1./sqrt(originalArea) );
			std::vector< ProjectiveData< Real , Real > > _simplcexCurvatures( simplices.size() );
			ThreadPool::ParallelFor
				(
					0 , simplices.size() ,
					[&]( unsigned int , size_t i )
					{
						Simplex< Real , Dim , K > s;
						Point< Real , Dim > n[K+1];
						Real curvatures[K];
						for( unsigned int k=0 ; k<=K ; k++ ) s[k] = Point< Real , Dim >( vertices[ simplices[i][k] ].template get< VERTEX_POSITION >()  ) * scl , n[k] = Point< Real , Dim > ( vertices[ simplices[i][k] ].template get< VERTEX_NORMAL >() );
						CurvatureMetric::PrincipalCurvatures< Real , K >( s , n , curvatures );

						Real area = s.measure() , kappa = 0;
						for( unsigned int k=0 ; k<K ; k++ ) kappa += curvatures[k] * curvatures[k];
						kappa = (Real)sqrt( kappa );
						_simplcexCurvatures[i] = ProjectiveData< Real , Real >( kappa * area , area );
					}
				);

			if( Verbose.set ) std::cout << pMeter( "Computed curvatures" ) << std::endl;

			std::vector< ProjectiveData< Real , Real > > _vertexCurvatures( vertices.size() );
			for( unsigned int i=0 ; i<simplices.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ ) _vertexCurvatures[ simplices[i][k] ] += _simplcexCurvatures[i];

			Real sigma = 0 , min = std::numeric_limits< Real >::infinity() , max = -std::numeric_limits< Real >::infinity();
			{
				ProjectiveData< Real , Real > _variance;

				for( unsigned int i=0 ; i<_vertexCurvatures.size() ; i++ )
				{
					Real k = _vertexCurvatures[i].value();
					_variance += ProjectiveData< Real , Real >( k * k * _vertexCurvatures[i].weight , _vertexCurvatures[i].weight );
					min = std::min< Real >( min , k );
					max = std::max< Real >( max , k );
				}
				sigma = (Real)sqrt( _variance.value() );
			}
			if( Verbose.set ) std::cout << "Standard deviation: " << sigma << " : [ " << min << " , " << max << " ]" << std::endl;

			auto CurvatureToColor = [&]( Real kappa )
				{
					Point< Real , 3 > clr;

					Real v = kappa / ( 2 * sigma );
					v = std::min< Real >( (Real)1 , std::max< Real >( (Real)-1 , v ) );

					clr = Point< Real , 3 >( (Real)1 , (Real)1 , (Real)1 ) * (1-v);
					return clr * 255;
				};

			for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get< VERTEX_COLOR >() = CurvatureToColor( _vertexCurvatures[i].value() );
		}
		// Adjust the metric to take into account the curvature
		///////////////////////////////////////////////////////

		////////////////////
		// Output the result
		if( Out.set ) PLY::WriteSimplices( Out.value , Factory() , vertices , simplices , file_type );
		// Output the result
		////////////////////
	}
	if( Verbose.set ) std::cout << "Processed: " << timer() << std::endl;

	return EXIT_SUCCESS;
}
