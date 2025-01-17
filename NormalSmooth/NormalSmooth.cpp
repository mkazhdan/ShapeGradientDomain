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
#include <Misha/CmdLineParser.h>
#include <Misha/Algebra.h>
#include <Misha/Ply.h>
#include <Misha/PlyVertexData.h>
#include <Misha/NormalSmoother.h>
#include <Misha/Miscellany.h>
#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif // EIGEN_USE_MKL_ALL
#include <Eigen/Sparse>

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< int > Iters( "iters" , 1 );
Misha::CmdLineParameter< float > DiffusionTimeStep( "dTime" , 1e-4f );
Misha::CmdLineReadable Verbose( "verbose" );

Misha::CmdLineReadable* params[] = { &In , &Out , &DiffusionTimeStep , &Verbose , &Iters , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name.c_str() );
	printf( "\t[--%s <output mesh>]\n" , Out.name.c_str() );
	printf( "\t[--%s <smoothing iterations>=%d]\n" , Iters.name.c_str() , Iters.value );
	printf( "\t[--%s <diffusion time step>=%f]\n" , DiffusionTimeStep.name.c_str() , DiffusionTimeStep.value );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

#ifdef EIGEN_USE_MKL_ALL
using Solver = Eigen::PardisoLLT< Eigen::SparseMatrix< double , Eigen::ColMajor , __int64 > >;
#else // !EIGEN_USE_MKL_ALL
using Solver = Eigen::SimplicialLLT< Eigen::SparseMatrix< double > >;
#endif // EIGEN_USE_MKL_ALL

template< class Real >
void _main( void )
{
	int file_type;
	std::vector< SimplexIndex< 2 > > triangles;
	std::vector< Point3D< Real > > vertices , normals;

	//////////////////////
	// Read in the data //
	{
		using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::NormalFactory< float , 3 > >;
		using Vertex = typename Factory::VertexType;
		Factory factory;

		std::vector< Vertex > _vertices;
		bool *vertexFlags = new bool[ factory.plyReadNum() ];
		PLY::ReadTriangles( In.value , factory , _vertices ,  triangles , vertexFlags , file_type );
		vertices.resize( _vertices.size() ) , normals.resize( _vertices.size() );
		for( int i=0 ; i<_vertices.size() ; i++ ) vertices[i] = Point3D< Real >( _vertices[i].template get<0>() );
		if( factory.template plyValidReadProperties< 1 >( vertexFlags ) ) for( int i=0 ; i<_vertices.size() ; i++ ) normals[i] = Point3D< Real >( _vertices[i].template get<1>() );
		else
		{
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				Point3D< Real > v[] = { Point3D< Real >( vertices[ triangles[i][0] ] ) , Point3D< Real >( vertices[ triangles[i][1] ] ) , Point3D< Real >( vertices[ triangles[i][2] ] ) };
				Point3D< Real > n = Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
				for( int j=0 ; j<3 ; j++ ) normals[ triangles[i][j] ] += n;
			}
			for( int i=0 ; i<normals.size() ; i++ ) normals[i] /= Point< Real , 3 >::Length( normals[i] );
		}
		delete[] vertexFlags;
		if( Verbose.set ) std::cout << "Source Vertices / Triangles: " << vertices.size() << " / " << triangles.size() << std::endl;
	}
	// Read in the data //
	//////////////////////

	/////////////////////
	// Smooth the normals
	{
		Miscellany::Timer timer;
		NormalSmoother::Smooth< 2 , Solver >( vertices , normals , triangles , Iters.value , DiffusionTimeStep.value );
		if( Verbose.set ) std::cout << "Smoothed normals: " << timer() << std::endl;
	}
	/////////////////////

	////////////////////
	// Output the result
	if( Out.set )
	{
		using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::NormalFactory< float , 3 > >;
		using Vertex = typename Factory::VertexType;

		std::vector< Vertex > _vertices( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = Point3D< float >( vertices[i] ) , _vertices[i].template get<1>() = Point3D< float >( normals[i] );
		PLY::WriteTriangles( Out.value , Factory() , _vertices , triangles , file_type );
	}
	// Output the result
	////////////////////
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	_main< double >();
	return EXIT_SUCCESS;
}
