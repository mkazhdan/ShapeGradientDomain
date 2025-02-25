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


enum struct SignalType
{
	POSITION ,
	NORMAL ,
	COLOR ,
	NORMAL_TO_POSITION ,
	COUNT
};
const std::string SignalNames[] = { std::string( "position" ) , std::string( "normal" ) , std::string( "color" ) , std::string( "normal to position" ) };

Misha::CmdLineParameter< std::string >
	In( "in" ) ,
	Out( "out" );

Misha::CmdLineParameter< float >
	CurvatureWeight( "kWeight" , 0.f ) ,
	NormalSmoothingWeight( "nWeight" , 1e-4f ) ,
	GradientWeight( "gWeight" , 1e-4f ) ,
	GradientScale( "gScale" , 1.f ) ,
	NormalProjectionWeight( "npWeight" , 1e2f );

Misha::CmdLineParameter< int >
	Signal( "signal" , static_cast< int >( SignalType::POSITION ) );

Misha::CmdLineReadable
	Anisotropic( "anisotropic" ) ,
	Verbose( "verbose" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&GradientWeight ,
	&GradientScale ,
	&NormalProjectionWeight ,
	&Anisotropic ,
	&Signal ,
	&CurvatureWeight ,
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
	printf( "\t[--%s <gradient weight>=%g]\n" , GradientWeight.name.c_str() , GradientWeight.value );
	printf( "\t[--%s <gradient scale>=%g]\n" , GradientScale.name.c_str() , GradientScale.value );
	printf( "\t[--%s <curvature weight>=%f]\n" , CurvatureWeight.name.c_str() , CurvatureWeight.value );
	printf( "\t[--%s <normal projection weight>=%g]\n" , NormalProjectionWeight.name.c_str() , NormalProjectionWeight.value );
	printf( "\t[--%s <signal type>=%d]\n" , Signal.name.c_str() , Signal.value );
	for( int i=0 ; i<static_cast< int >( SignalType::COUNT ) ; i++ ) printf( "\t\t%d] %s\n" , i , SignalNames[i].c_str() );
	printf( "\t[--%s]\n" , Anisotropic.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

#ifdef EIGEN_USE_MKL_ALL
using Solver = Eigen::PardisoLLT< Eigen::SparseMatrix< double , Eigen::ColMajor , __int64 > >;
#else // !EIGEN_USE_MKL_ALL
using Solver = Eigen::SimplicialLLT< Eigen::SparseMatrix< double > >;
#endif // EIGEN_USE_MKL_ALL

enum{ VERTEX_POSITION , VERTEX_NORMAL , VERTEX_COLOR };

using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::NormalFactory< float , 3 > , VertexFactory::RGBColorFactory< float > >;
using Vertex = typename Factory::VertexType;


template< bool HasNormals , bool HasColors > struct OutPLY;

template<>
struct OutPLY< true , true >
{
	static const unsigned int NormalIndex = 1 , ColorIndex = 1;
	using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::NormalFactory< float , 3 > , VertexFactory::RGBColorFactory< float > >;
};

template<>
struct OutPLY< true , false >
{
	static const unsigned int NormalIndex = 1 , ColorIndex = -1;
	using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::NormalFactory< float , 3 > >;
};

template<>
struct OutPLY< false , true >
{
	static const unsigned int NormalIndex = -1 , ColorIndex = 1;
	using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::RGBColorFactory< float > >;
};

template<>
struct OutPLY< false , false >
{
	static const unsigned int NormalIndex = -1 , ColorIndex = -1;
	using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > >;
};

template< bool HasNormals , bool HasColors >
void WriteMesh( std::string fileName , const std::vector< Vertex > &vertices , const std::vector< TriangleIndex > &triangles , int file_type )
{
	using OutFactory = typename OutPLY< HasNormals , HasColors >::Factory;
	using OutVertex = typename OutFactory::VertexType;
	static const unsigned int NormalIndex = OutPLY< HasNormals , HasColors >::NormalIndex;
	static const unsigned int ColorIndex = OutPLY< HasNormals , HasColors >::ColorIndex;

	std::vector< OutVertex > outVertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ )
	{
		outVertices[i].template get<0>() = vertices[i].template get< VERTEX_POSITION >();
		if constexpr( HasNormals ) outVertices[i].template get< NormalIndex >() = vertices[i].template get< VERTEX_NORMAL >();
		if constexpr( HasColors ) outVertices[i].template get< ColorIndex >() = vertices[i].template get< VERTEX_COLOR >();
	}
	PLY::WriteTriangles( fileName , OutFactory() , outVertices , triangles , file_type );
}

template< class Real >
int _main( void )
{
	std::vector< TriangleIndex > triangles;
	std::vector< Vertex > vertices;

	int file_type;
	bool hasNormals , hasColors;

	Solver solver;

	Miscellany::PerformanceMeter pMeter( '.' );

	//////////////////////
	// Read in the data //
	pMeter.reset();
	{
		Factory factory;
		bool *readFlags = new bool[ factory.plyReadNum() ];
		file_type = PLY::ReadTriangles( In.value , factory , vertices , triangles , readFlags );

		hasNormals = factory.template plyValidReadProperties< VERTEX_NORMAL >( readFlags );
		hasColors  = factory.template plyValidReadProperties< VERTEX_COLOR  >( readFlags );

		// Create normals if they are not already there
		if( !hasNormals )
		{
			for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<1>() = Point3D< float >();
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				Point3D< float > n = Point3D< float >::CrossProduct( vertices[ triangles[i][1] ].template get< VERTEX_POSITION >() - vertices[ triangles[i][0] ].template get< VERTEX_POSITION >() , vertices[ triangles[i][2] ].template get< VERTEX_POSITION >() - vertices[ triangles[i][0] ].template get< VERTEX_POSITION >() );
				for( int j=0 ; j<3 ; j++ ) vertices[ triangles[i][j] ].template get< VERTEX_NORMAL >() += n;
			}
			for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get< VERTEX_NORMAL >() /= Point< float , 3 >::Length( vertices[i].template get<1>() );
		}

		if( SignalType( Signal.value )==SignalType::NORMAL ) hasNormals = true;
		if( SignalType( Signal.value )==SignalType::COLOR && !hasColors ) ERROR_OUT( "No color channel specified" );
		delete[] readFlags;
		if( Verbose.set ) std::cout << "Source Vertices / Triangles: " << vertices.size() << " / " << triangles.size() << std::endl;
	}
	if( Verbose.set ) std::cout << pMeter( "Read mesh" ) << std::endl;
	// Read in the data //
	//////////////////////

	//////////////////////////
	// Set the Riemannian mesh
	pMeter.reset();
	Real originalArea = 0;
	FEM::RiemannianMesh< Real > mesh( GetPointer( triangles ) , triangles.size() );
	{
		// Create the embedded metric and normalize to have unit are
		mesh.template setMetricFromEmbedding< 3 >( [&]( unsigned int i ){ return vertices[i].template get<0>(); } );
		originalArea = mesh.area();
		mesh.makeUnitArea();
	}
	if( Verbose.set ) std::cout << pMeter( "Set Riemannian mesh" ) << std::endl;
	// Set the Riemannian mesh
	//////////////////////////

	///////////////////////////////////////
	// System matrix symbolic factorization
	pMeter.reset();
	solver.analyzePattern( mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY , true >() );
	if( Verbose.set ) std::cout << pMeter( "Symbolic factorization" ) << std::endl;
	// System matrix symbolic factorization
	///////////////////////////////////////

	///////////////////////////////////////////////////////
	// Adjust the metric to take into account the curvature
	if( CurvatureWeight.value>0 )
	{
		pMeter.reset();
		std::vector< Point< Real , 3 > > curvatureNormals( vertices.size() );
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) curvatureNormals[i] = vertices[i].template get<1>();

		///////////////////
		// Normal smoothing
		if( NormalSmoothingWeight.value>0 )
		{
			Miscellany::PerformanceMeter pMeter( '.' );
			curvatureNormals = GradientDomain::ProcessVertexVertex< Solver , Point3D< Real > , Real >( solver , mesh , (Real)1. , NormalSmoothingWeight.value , [&]( unsigned int v ){ return curvatureNormals[v]; } , []( unsigned int ){ return Point< Real , 3 >(); } );
			ThreadPool::ParallelFor
			(
				0 , curvatureNormals.size() ,
				[&]( unsigned int , size_t i ){ curvatureNormals[i] /= Point< Real , 3 >::Length( curvatureNormals[i] ); }
			);
			if( Verbose.set ) std::cout << pMeter( "Normal smoothing" ) << std::endl;
		}
		// Normal smoothing
		///////////////////

		{
			Miscellany::PerformanceMeter pMeter( '.' );

			// Input: Principal curvature values
			// Output: Positive entries of the diagonal matrix describing the scaling along the principal curvature directions
			//         Outputting the identity matrix reproduces the embedding metric
			auto PrincipalCurvatureFunctor = [&]( unsigned int , Point< Real , 2 > pCurvatures )
				{
					pCurvatures[0] *= pCurvatures[0];
					pCurvatures[1] *= pCurvatures[1];
					if( !Anisotropic.set ) pCurvatures[0] = pCurvatures[1] = ( pCurvatures[0] + pCurvatures[1] ) / (Real)2.;
					return Point< Real , 2 >( (Real)1. , (Real)1. ) + pCurvatures * CurvatureWeight.value;
				};

			Real s = (Real)(1./sqrt(originalArea) );
			CurvatureMetric::SetCurvatureMetric
			(
				mesh ,
				[&]( unsigned int idx ){ return Point< Real , 3 >( vertices[idx].template get< VERTEX_POSITION >() ) * s; } ,
				[&]( unsigned int idx ){ return Point< Real , 3 >( curvatureNormals[idx] ); } ,
				PrincipalCurvatureFunctor
			);
			if( Verbose.set ) std::cout << pMeter( "Metric update" ) << std::endl;
		}
		if( Verbose.set ) std::cout << pMeter( "Adjusted metric" ) << std::endl;
	}
	// Adjust the metric to take into account the curvature
	///////////////////////////////////////////////////////


	///////////////
	// Set the data
	std::vector< Point3D< Real > > signal( vertices.size() );
	switch( SignalType( Signal.value ) )
	{
	case SignalType::POSITION: for( int i=0 ; i<vertices.size() ; i++ ) signal[i] = Point3D< Real >( vertices[i].template get< VERTEX_POSITION >() ) ; break;
	case SignalType::NORMAL_TO_POSITION:
	case SignalType::NORMAL  : for( int i=0 ; i<vertices.size() ; i++ ) signal[i] = Point3D< Real >( vertices[i].template get< VERTEX_NORMAL   >() ) ; break;
	case SignalType::COLOR   : for( int i=0 ; i<vertices.size() ; i++ ) signal[i] = Point3D< Real >( vertices[i].template get< VERTEX_COLOR    >() ) ; break;
	default: ERROR_OUT( "Unrecognized signal type: " , Signal.value );
	}
	// Set the data
	///////////////

	////////////////////
	// Smooth the signal
	if( GradientWeight.value>0 && GradientScale.value!=1 )
	{
		pMeter.reset();
		{
			// Low frequencies described in terms of values at vertices
			auto LowFrequencyVertexValues = [&]( unsigned int v ){ return signal[v]; };

			// High frequencies described in terms of scaled values at vertices
			auto HighFrequencyVertexValues = [&]( unsigned int v ){ return signal[v]*(Real)GradientScale.value; };

			signal = GradientDomain::ProcessVertexVertex< Solver , Point3D< Real > , Real >( solver , mesh , (Real)1. , (Real)GradientWeight.value , LowFrequencyVertexValues , HighFrequencyVertexValues );
		}
		if( Verbose.set ) std::cout << pMeter( "Processed" ) << std::endl;
	}
	// Smooth the signal
	////////////////////

	//////////////////////////////////
	// Fit the geometry to the normals
	if( SignalType( Signal.value )==SignalType::NORMAL_TO_POSITION )
	{
		pMeter.reset();
		signal = GradientDomain::FitToNormals< Solver >( solver , mesh , (Real)1. , (Real)NormalProjectionWeight.value , [&]( unsigned int v ){ return vertices[v].template get< VERTEX_POSITION >(); } , [&]( unsigned int v ){ return signal[v]; } );
		if( Verbose.set ) std::cout << pMeter( "Fit geometry" ) << std::endl;
	}
	// Fit the geometry to the normals
	//////////////////////////////////

	////////////////////////////////////////////
	// Copy the processed signal into the vertices
	switch( SignalType( Signal.value ) )
	{
	case SignalType::NORMAL_TO_POSITION:
	case SignalType::POSITION: for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get< VERTEX_POSITION >() = Point3D< float >( signal[i] ) ; break;
	case SignalType::NORMAL  : for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get< VERTEX_NORMAL   >() = Point3D< float >( signal[i] ) ; break;
	case SignalType::COLOR   : for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get< VERTEX_COLOR    >() = Point3D< float >( signal[i] ) ; break;
	default: ERROR_OUT( "Unrecognized singal type: " , Signal.value );
	}
	// Copy the processed signal to the vertices
	////////////////////////////////////////////

	////////////////////
	// Output the result
	if( Out.set )
	{
		if( hasColors && hasNormals ) WriteMesh< true  , true  >( Out.value , vertices , triangles , file_type );
		else if( hasColors )          WriteMesh< false , true  >( Out.value , vertices , triangles , file_type );
		else if( hasNormals )         WriteMesh< true  , false >( Out.value , vertices , triangles , file_type );
		else                          WriteMesh< false , false >( Out.value , vertices , triangles , file_type );
	}
	// Output the result
	////////////////////
	return EXIT_SUCCESS;
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	Miscellany::Timer timer;
	_main< double >();
	if( Verbose.set ) std::cout << "Processed: " << timer() << std::endl;

	return EXIT_SUCCESS;
}
