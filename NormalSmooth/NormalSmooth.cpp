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

// Enable this to force testing of array access

#include <stdio.h>
#include <stdlib.h>
#ifdef NEW_CODE
#else // !NEW_CODE
#include <omp.h>
#endif // NEW_CODE
#include <algorithm>
#include <vector>
#include <Misha/CmdLineParser.h>
#include <Misha/Algebra.h>
#include <Misha/Ply.h>
#include <Misha/PlyVertexData.h>
#ifdef NEW_CODE
#include <Misha/NormalSmoother.h>
#include <Misha/Miscellany.h>
#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif // EIGEN_USE_MKL_ALL
#include <Eigen/Sparse>
#else // !NEW_CODE
#include <Misha/LinearSolvers.h>
#include <Misha/Timer.h>
#include <Misha/FEM.h>
#endif // NEW_CODE

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" ) , Sphere( "sphere" );
Misha::CmdLineParameter< int > Iters( "iters" , 1 );
Misha::CmdLineParameter< float > ValueWeight( "vWeight" , 1e4f ) , GradientWeight( "gWeight" , 1.f ) , GradientScale( "gScale" , 1.f );
Misha::CmdLineReadable Verbose( "verbose" );

Misha::CmdLineReadable* params[] = { &In , &Out , &ValueWeight , &GradientWeight , &Verbose , &Iters , NULL };

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name.c_str() );
	printf( "\t[--%s <output mesh>]\n" , Out.name.c_str() );
	printf( "\t[--%s <smoothing iterations>=%d]\n" , Iters.name.c_str() , Iters.value );
	printf( "\t[--%s <value weight>=%f]\n" , ValueWeight.name.c_str() , ValueWeight.value );
#ifndef FOR_RELEASE
	printf( "\t[--%s <gradient/smoothing weight>=%f]\n" , GradientWeight.name.c_str() , GradientWeight.value );
#endif // !FOR_RELEASE
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

#ifdef NEW_CODE
#ifdef EIGEN_USE_MKL_ALL
using Solver = Eigen::PardisoLLT< Eigen::SparseMatrix< double , Eigen::ColMajor , __int64 > >;
#else // !EIGEN_USE_MKL_ALL
using Solver = Eigen::SimplicialLLT< Eigen::SparseMatrix< double > >;
#endif // EIGEN_USE_MKL_ALL
#else // !NEW_CODE
template< class Real >
void WriteMesh( std::string fileName , const std::vector< TriangleIndex >& triangles , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& normals )
{
	using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::RGBColorFactory< float > >;
	using Vertex = typename Factory::VertexType;

	std::vector< Vertex > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = Point3D< float >( vertices[i] ) , _vertices[i].template get<1>() = Point3D< float >( normals[i] + Point3D< Real >( 1 , 1 , 1 ) ) * 128.f;
	PLY::WriteTriangles( fileName , Factory() ,  _vertices ,  triangles , PLY_BINARY_NATIVE );
}

template< class Real >
void WriteColorMesh( char* fileName , const std::vector< TriangleIndex >& triangles , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& colors )
{
	using Factory = VertexFactory::Factory< float , VertexFactory::PositionFactory< float , 3 > , VertexFactory::RGBColorFactory< float > >;
	using Vertex = typename Factory::VertexType;

	std::vector< Vertex > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].template get<0>() = Point3D< float >( vertices[i] ) , _vertices[i].template get<1>() = Point3D< float >( colors[i] );
	PlyWriteTriangles( fileName , Factory() , _vertices ,  triangles , PLY_BINARY_NATIVE );
}

template< class Real >
void WriteMesh( char* fileName , const std::vector< TriangleIndex >& triangles , const std::vector< Point3D< Real > >& vertices )
{
	using Factory = VertexFactory::PositionFactory< float , 3 >;
	using Vertex = typename Factory::VertexType;

	std::vector< Vertex > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i] = Point3D< float >( vertices[i] );
	PlyWriteTriangles( fileName , Factory() , _vertices ,  triangles , PLY_BINARY_NATIVE );
}
#endif // NEW_CODE

template< class Real >
void _main( void )
{
#ifdef NEW_CODE
#else // !NEW_CODE
	typedef EigenSolverCholeskyLLt< Real , typename SparseMatrix< Real , int >::RowIterator > Solver;
#endif // NEW_CODE
	int file_type;
#ifdef NEW_CODE
	std::vector< SimplexIndex< 2 > > triangles;
#else // !NEW_CODE
	std::vector< TriangleIndex > triangles;
#endif // NEW_CODE
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

#ifdef NEW_CODE
	/////////////////////
	// Smooth the normals
	{
		Miscellany::Timer timer;
		NormalSmoother::Smooth< 2 , Solver >( vertices , normals , triangles , Iters.value , 1./ValueWeight.value );
		if( Verbose.set ) std::cout << "Smoothed normals: " << timer() << std::endl;
	}
	/////////////////////
#else // !NEW_CODE
	/////////////////////
	// Smooth the normals
	{
		FEM::RiemannianMesh< Real > mesh( GetPointer( triangles ) , triangles.size() );
		mesh.setMetricFromEmbedding( GetPointer( vertices ) );
		mesh.makeUnitArea();
		SparseMatrix< Real , int > M , _M = mesh.template massMatrix< FEM::BASIS_0_WHITNEY >() , _S = mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY >();
		M.resize( 2*vertices.size() );
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) for( int ii=0 ; ii<2 ; ii++ )
		{
			M.SetRowSize( 2*i+ii , 2*_M.rowSizes[i] );
			for( int j=0 ; j<_M.rowSizes[i] ; j++ ) for( int jj=0 ; jj<2 ; jj++ ) M[2*i+ii][2*j+jj].N = _M[i][j].N*2+jj;
		}
		std::vector< Point3D< Real > > tangents( vertices.size()*2 );
		std::vector< Real > b( vertices.size()*2 ) , o( vertices.size()*2 );

		Solver solver( M , true );

		Real mWeight = (Real)ValueWeight.value , sWeight = (Real)GradientWeight.value;

		for( int iter=0 ; iter<Iters.value ; iter++ )
		{
			Timer t;

			// Set the tangent directions
#pragma omp parallel for
			for( int i=0 ; i<vertices.size() ; i++ )
			{
				Point3D< Real > v( 1 , 0 , 0 );
				if( fabs( Point3D< Real >::Dot( v , normals[i] ) )>0.99 ) v = Point3D< Real >( 0 , 1 , 0 );
				tangents[2*i+0] = Point3D< Real >::CrossProduct( normals[i] , v               ) ; tangents[2*i+0] /= Point< Real , 3 >::Length( tangents[2*i+0] );
				tangents[2*i+1] = Point3D< Real >::CrossProduct( normals[i] , tangents[2*i+0] ) ; tangents[2*i+1] /= Point< Real , 3 >::Length( tangents[2*i+1] );
			}

			// Solve for the tangent offsets minimizing the dirichlet energy:
			// E( o1 , o2 ) = || \sum o[i] * T[i] ||^2 + e * || \nabla( \sum n[i] + o[i] * T[i] ) ||^2
			//              = o^t * T^t * M * T * o + e * [ o^t * T^t * S * T * o + 2 * o^t * T^t * S * n + n^t * S * n ]
			// \nabla E = 0:
			// 0 = T^t * ( M + e * S ) * T * o + e * T^t * S * n
			{
				Timer t;
#pragma omp parallel for 
				for( int i=0 ; i<vertices.size() ; i++ ) for( int ii=0 ; ii<2 ; ii++ ) 
				{
					b[2*i+ii] = 0;
					for( int j=0 ; j<_M.rowSizes[i] ; j++ )
					{
						for( int jj=0 ; jj<2 ; jj++ ) M[2*i+ii][2*j+jj].Value = ( _M[i][j].Value*mWeight + _S[i][j].Value * sWeight ) * Point3D< Real >::Dot( tangents[2*i+ii] , tangents[ 2*_M[i][j].N+jj ] );
						b[2*i+ii] -= _S[i][j].Value * Point3D< Real >::Dot( normals[ _S[i][j].N ] , tangents[2*i+ii] ) * sWeight;
					}
				}
				if( Verbose.set ) printf( "\tSet system matrix[%d]: %.2f(s)\n" , iter+1 , t.elapsed() );
			}
			{
				Timer t;
				solver.update( M );
				solver.solve( GetPointer( b ) , GetPointer( o ) );
#pragma omp parallel for
				for( int i=0 ; i<vertices.size() ; i++ ) normals[i] += tangents[2*i+0] * o[2*i+0] + tangents[2*i+1] * o[2*i+1] , normals[i] /= Point< Real , 3 >::Length( normals[i] );
				if( Verbose.set ) printf( "\tSolved system[%d]: %.2f(s)\n" , iter+1 , t.elapsed() );
			}
		}
	}
	// Smooth the normals
	/////////////////////
#endif // NEW_CODE

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
