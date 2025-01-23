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
#include <Misha/NormalSmoother.h>

#include <Include/CurvatureMetric.h>
#include <Include/GradientDomain.h>

Misha::CmdLineParameter< std::string >
	In( "in" ) ,
	Out( "out" );

Misha::CmdLineParameter< float >
	NormalSmoothingDiffusionTime( "nTime" , 1e-4f ) ,
	CurvatureWeight( "kWeight" , 0.f ) ,
	ValueWeight( "vWeight" , 1.f ) ,
	GradientWeight( "gWeight" , 1e-4f ) ,
	GradientScale( "gScale" , 1.f );

Misha::CmdLineParameter< int >
	NormalSmoothingIters( "nIters" , 0 );

Misha::CmdLineReadable
	UseEdgeDifferences( "edgeDifferences" ) ,
	Anisotropic( "anisotropic" ) ,
	UseColors( "useColors" ) ,
	Verbose( "verbose" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&ValueWeight ,
	&GradientWeight ,
	&GradientScale ,
	&Anisotropic ,
	&UseColors ,
	&CurvatureWeight ,
	&NormalSmoothingIters ,
	&NormalSmoothingDiffusionTime ,
	&UseEdgeDifferences ,
	&Verbose ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input geometry>\n" , In.name.c_str() );
	printf( "\t[--%s <output geometry>]\n" , Out.name.c_str() );
	printf( "\t[--%s <value weight>=%f]\n" , ValueWeight.name.c_str() , ValueWeight.value );
	printf( "\t[--%s <gradient weight>=%f]\n" , GradientWeight.name.c_str() , GradientWeight.value );
	printf( "\t[--%s <gradient scale>=%f]\n" , GradientScale.name.c_str() , GradientScale.value );
	printf( "\t[--%s <curvature weight>=%f]\n" , CurvatureWeight.name.c_str() , CurvatureWeight.value );
	printf( "\t[--%s <normal smoothing iterations>=%d]\n" , NormalSmoothingIters.name.c_str() , NormalSmoothingIters.value );
	printf( "\t[--%s <normal smoothing diffusion time>=%g]\n" , NormalSmoothingDiffusionTime.name.c_str() , NormalSmoothingDiffusionTime.value );
	printf( "\t[--%s]\n" , Anisotropic.name.c_str() );
	printf( "\t[--%s]\n" , UseEdgeDifferences.name.c_str() );
	printf( "\t[--%s]\n" , UseColors.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

#ifdef EIGEN_USE_MKL_ALL
using Solver = Eigen::PardisoLLT< Eigen::SparseMatrix< double , Eigen::ColMajor , __int64 > >;
#else // !EIGEN_USE_MKL_ALL
using Solver = Eigen::SimplicialLLT< Eigen::SparseMatrix< double > >;
#endif // EIGEN_USE_MKL_ALL

template< class Real > Real Clamp( Real v , Real min , Real max ){ return std::min< Real >( max , std::max< Real >( min , v ) ); }
template< class Real > Point3D< Real > ClampColor( Point3D< Real > c ){ return Point3D< Real >( Clamp( c[0] , (Real)0 , (Real)255 ) , Clamp( c[1] , (Real)0 , (Real)255 ) , Clamp( c[2] , (Real)0 , (Real)255 ) ); }

template< class Real >
int _main( void )
{
	using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::NormalFactory< float , 3 > , VertexFactory::RGBColorFactory< float > >;
	using Vertex = typename Factory::VertexType;

	std::vector< TriangleIndex > triangles;
	std::vector< Vertex > vertices;

	int file_type;
	bool hasNormals , hasColors;

	//////////////////////
	// Read in the data //
	{
		Factory factory;
		bool *readFlags = new bool[ factory.plyReadNum() ];
		PLY::ReadTriangles( In.value , factory , vertices , triangles , readFlags , file_type );

		hasNormals = factory.template plyValidReadProperties<1>( readFlags );
		hasColors  = factory.template plyValidReadProperties<2>( readFlags );

		if( !hasNormals )
		{
			for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<1>() = Point3D< float >();
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				Point3D< float > n = Point3D< float >::CrossProduct( vertices[ triangles[i][1] ].template get<0>() - vertices[ triangles[i][0] ].template get<0>() , vertices[ triangles[i][2] ].template get<0>() - vertices[ triangles[i][0] ].template get<0>() );
				for( int j=0 ; j<3 ; j++ ) vertices[ triangles[i][j] ].template get<1>() += n;
			}
			for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<1>() /= Point< float , 3 >::Length( vertices[i].template get<1>() );
		}
		if( !hasColors ) for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<2>() = ( Point3D< float >( 1.f , 1.f , 1.f ) + vertices[i].template get<1>() ) / 2.f;
		delete[] readFlags;
		if( Verbose.set ) std::cout << "Source Vertices / Triangles: " << vertices.size() << " / " << triangles.size() << std::endl;
	}
	// Read in the data //
	//////////////////////

	///////////////////
	// Normal smoothing
	if( NormalSmoothingIters.value>0 && NormalSmoothingDiffusionTime.value>0 )
	{
		Miscellany::Timer timer;
		std::vector< Point3D< Real > > v( vertices.size() ) , n( vertices.size() );
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) v[i] = vertices[i].template get<0>() , n[i] = vertices[i].template get<1>();
		NormalSmoother::Smooth< 2 , int , Solver >( v , n , triangles , NormalSmoothingIters.value , NormalSmoothingDiffusionTime.value );
		for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<1>() = n[i];
		if( Verbose.set ) std::cout << "\tSmoothed normals: " << timer() << std::endl;
	}
	// Normal smoothing
	///////////////////
	
	std::vector< Point3D< Real > > signal( vertices.size() );

	//////////////////
	// Set the data //
	if( UseColors.set ) for( int i=0 ; i<vertices.size() ; i++ ) signal[i] = Point3D< Real >( vertices[i].template get<2>() );
	else                for( int i=0 ; i<vertices.size() ; i++ ) signal[i] = Point3D< Real >( vertices[i].template get<0>() );
	// Set the data //
	//////////////////

	////////////////////
	// Smooth the signal
	{
		// Read in the mesh and normalize the scale
		FEM::RiemannianMesh< Real > mesh( GetPointer( triangles ) , triangles.size() );
		mesh.template setMetricFromEmbedding< 3 >( [&]( unsigned int i ){ return vertices[i].template get<0>(); } );
		Real area = mesh.area();
		mesh.makeUnitArea();

		////////////////////
		// Modify the metric
		if( CurvatureWeight.value>0 )
		{
			Miscellany::Timer timer;
			auto PrincipalCurvatureFunctor = [&]( Point< Real , 2 > pCurvatures )
				{
					pCurvatures[0] *= pCurvatures[0];
					pCurvatures[1] *= pCurvatures[1];
					if( !Anisotropic.set ) pCurvatures[0] = pCurvatures[1] = ( pCurvatures[0] + pCurvatures[1] ) / (Real)2.;
					return Point< Real , 2 >( (Real)1. , (Real)1. ) + pCurvatures * CurvatureWeight.value;
				};
			Real s = (Real)(1./sqrt(area) );
			CurvatureMetric::SetCurvatureMetric
			(
				mesh ,
				[&]( unsigned int idx ){ return Point< Real , 3 >( vertices[idx].template get<0>() ) * s; } ,
				[&]( unsigned int idx ){ return Point< Real , 3 >( vertices[idx].template get<1>() ); } ,
				PrincipalCurvatureFunctor
			);
			if( Verbose.set ) std::cout << "\tUpdated metric: " << timer() << std::endl;
		}
		// Modify the metric
		////////////////////

		///////////////////////////////
		// Compute and solve the system
		{
			Miscellany::Timer timer;

			auto LowFrequencyVertexValues = [&]( unsigned int v ){ return signal[v]; };
			if( UseEdgeDifferences.set )
			{
				auto HighFrequencyEdgeValues = [&]( unsigned int e )
					{
						int v1 , v2;
						mesh.edgeVertices( e , v1 , v2 );
						return ( signal[v2] - signal[v1] ) * GradientScale.value;
					};
				signal = GradientDomain::ProcessVertexEdge< Solver , Point3D< Real > , Real >( mesh , ValueWeight.value , GradientWeight.value , LowFrequencyVertexValues , HighFrequencyEdgeValues );
			}
			else
			{
				auto HighFrequencyVertexValues = [&]( unsigned int v ){ return signal[v]*GradientScale.value; };
				signal = GradientDomain::ProcessVertexVertex< Solver , Point3D< Real > , Real >( mesh , ValueWeight.value , GradientWeight.value , LowFrequencyVertexValues , HighFrequencyVertexValues );
			}
			if( Verbose.set ) std::cout << "\tSolved the system: " << timer() << std::endl;
		}
		// Compute and solve the system
		///////////////////////////////
	}
	// Smooth the signal
	////////////////////

	////////////////////////////////////////////
	// Copy the processed signal to the vertices
	if( UseColors.set ){ hasColors = true ; for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<2>() = ClampColor( Point3D< float >( signal[i] ) ); }
	else                                    for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<0>() = Point3D< float >( signal[i] );
	// Copy the processedsignal to the vertices
	////////////////////////////////////////////

	////////////////////
	// Output the result
	if( Out.set )
	{
		if( hasColors && hasNormals ) PLY::WriteTriangles( Out.value , Factory() , vertices , triangles , file_type );
		else if( hasColors )
		{
			using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::RGBColorFactory< float > >;
			using Vertex = typename Factory::VertexType;

			std::vector< Vertex > _vertices( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i].template get<0>() , _vertices[i].template get<1>() = vertices[i].template get<2>();
			PLY::WriteTriangles( Out.value , Factory() , _vertices , triangles , file_type );
		}
		else if( hasNormals )
		{
			using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::NormalFactory< float , 3 > >;
			using Vertex = typename Factory::VertexType;

			std::vector< Vertex > _vertices( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = vertices[i].template get<0>() , _vertices[i].template get<1>() = vertices[i].template get<1>();
			PLY::WriteTriangles( Out.value , Factory() , _vertices , triangles , file_type );
		}
		else
		{
			using Factory = VertexFactory::PositionFactory< float , 3 >;
			using Vertex = typename Factory::VertexType;

			std::vector< Vertex > _vertices( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[i] = vertices[i].template get<0>();
			PLY::WriteTriangles( Out.value , Factory() , _vertices , triangles , file_type );
		}
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
	return _main< double >();
}
