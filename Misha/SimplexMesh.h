/*
Copyright (c) 2022, Michael Kazhdan
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

#ifndef SIMPLEX_MESH_INCLUDED
#define SIMPLEX_MESH_INCLUDED
#include <iostream>
#include <vector>
#include <functional>
#include "Geometry.h"
#include "SimplexBasis.h"
#include "MultiThreading.h"

template< unsigned int ... > struct SimplexMesh;

template< unsigned int Dim >
struct SimplexMesh< Dim >
{
	struct Sample
	{
		unsigned int sIdx;
		Point< double , Dim+1 > bcCoordinates;
	};
};

template< unsigned int Dim , unsigned int Degree >
struct SimplexMesh< Dim , Degree >
{
	typedef MultiIndex< Degree , unsigned int , false > NodeMultiIndex;
	typedef MultiIndex< Dim , unsigned int , false > FaceMultiIndex;
	static const unsigned int NodesPerSimplex = SimplexElements< Dim , Degree >::NodeNum;

	SimplexMesh( void ){}

	template< unsigned int EmbeddingDimension , typename Index >
	using VertexFunctor = std::function< Point< double , EmbeddingDimension > ( Index ) >;
	template< typename Index >
	using MetricFunctor = std::function< SquareMatrix< double , Dim > ( Index ) >;

	template< unsigned int EmbeddingDimension , typename Index >
	static SimplexMesh Init( const std::vector< SimplexIndex< Dim , Index > > &simplices , VertexFunctor< EmbeddingDimension , Index > vFunction );

	template< typename Index >
	static SimplexMesh Init( const std::vector< SimplexIndex< Dim , Index > > &simplices , MetricFunctor< Index > gFunction );

	size_t simplices( void ) const { return _simplices.size(); }
	SimplexIndex< Dim , unsigned int > simplex( unsigned int idx ) const{ return _simplices[idx]; }
	SquareMatrix< double , Dim > metric( unsigned int idx ) const { return _g[idx]; }
	size_t nodes( void ) const { return _nodeMap.size(); }
	Eigen::SparseMatrix< double > mass( void ) const;
	Eigen::SparseMatrix< double > stiffness( void ) const;
	Eigen::SparseMatrix< double > bistiffness( void ) const;
	Eigen::SparseMatrix< double > system( std::function< SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum > ( const SquareMatrix< double , Dim > & ) > m2s ) const;
	Eigen::SparseMatrix< double > crossFaceGradientEnergy( void ) const;
	template< typename UseFaceFunctor >
	Eigen::SparseMatrix< double > crossFaceGradientEnergy( const UseFaceFunctor &useFaceFunctor ) const;
	template< unsigned int FaceDim >
	Eigen::SparseMatrix< double > _crossFaceGradientEnergy( void ) const;
	Eigen::SparseMatrix< double > evaluationMatrix( const std::vector< typename SimplexMesh< Dim >::Sample > &samples ) const;
	double evaluate( const Eigen::VectorXd &coefficients , typename SimplexMesh< Dim >::Sample sample ) const;
	Polynomial::Polynomial< Dim , Degree , double > evaluate( const Eigen::VectorXd &coefficients , unsigned int simplexIndex ) const;

	SimplexMesh( SimplexMesh &&sm ){ std::swap( _simplices , sm._simplices ) , std::swap( _g , sm._g ) , std::swap( _nodeMap , sm._nodeMap ); }
	SimplexMesh & operator = ( SimplexMesh &&sm ){ std::swap( _simplices , sm._simplices ) , std::swap( _g , sm._g ) , std::swap( _nodeMap , sm._nodeMap ) ; return *this; }

	NodeMultiIndex nodeMultiIndex( unsigned int s , unsigned int n ) const;
	unsigned int nodeIndex( const NodeMultiIndex &multiIndex ) const;
	unsigned int nodeIndex( unsigned int s , unsigned int n ) const { return _localToGlobalNodeIndex.size() ? _localToGlobalNodeIndex[ s*NodesPerSimplex+n ] : nodeIndex( nodeMultiIndex( s , n ) ); }
	typename NodeMultiIndex::map::const_iterator cbegin( void ) const { return _nodeMap.cbegin(); }
	typename NodeMultiIndex::map::const_iterator cend  ( void ) const { return _nodeMap.cend  (); }
	const typename NodeMultiIndex::map &nodeMap( void ) const { return _nodeMap; }

	double volume( void ) const;
	void makeUnitVolume( void );
	void setMetric( std::function< SquareMatrix< double , Dim > (unsigned int) > metricFunction );
	void hashLocalToGlobalNodeIndex( void );

protected:
	template< unsigned int EmbeddingDimension , typename Index >
	void _init( const std::vector< SimplexIndex< Dim , Index > > &simplices , VertexFunctor< EmbeddingDimension , Index > vFunction );

	template< typename Index >
	void _init( const std::vector< SimplexIndex< Dim , Index > > &simplices , MetricFunctor< Index > gFunction );

	std::vector< SimplexIndex< Dim , unsigned int > > _simplices;
	typename NodeMultiIndex::map _nodeMap;
	std::vector< SquareMatrix< double , Dim > > _g;
	std::vector< unsigned int > _localToGlobalNodeIndex;
};

#include "SimplexMesh.inl"
#endif // SIMPLEX_MESH_INCLUDED