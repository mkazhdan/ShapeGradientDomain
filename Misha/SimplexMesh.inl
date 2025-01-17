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

/////////////////
// SimplexMesh //
/////////////////
template< unsigned int Dim , unsigned int Degree >
template< unsigned int EmbeddingDimension , typename Index >
SimplexMesh< Dim , Degree > SimplexMesh< Dim , Degree >::Init( const std::vector< SimplexIndex< Dim , Index > > &simplices , VertexFunctor< EmbeddingDimension , Index > vFunction )
{
	SimplexMesh sm;
	sm._init( simplices , vFunction );
	return sm;
}

template< unsigned int Dim , unsigned int Degree >
template< typename Index >
SimplexMesh< Dim , Degree > SimplexMesh< Dim , Degree >::Init( const std::vector< SimplexIndex< Dim , Index> > &simplices , MetricFunctor< Index > gFunction )
{
	SimplexMesh sm;
	sm._init( simplices , gFunction );
	return sm;
}

template< unsigned int Dim , unsigned int Degree >
template< unsigned int EmbeddingDimension , typename Index >
void SimplexMesh< Dim , Degree >::_init( const std::vector< SimplexIndex< Dim , Index> > &simplices , VertexFunctor< EmbeddingDimension , Index > vFunction )
{
	std::function< SquareMatrix< double , Dim > ( unsigned int ) > gFunction = [&]( unsigned int s )
	{
		Simplex< double , EmbeddingDimension , Dim > simplex;
		for( unsigned int d=0 ; d<=Dim ; d++ ) simplex[d] = vFunction( simplices[s][d] );
		return RightSimplex< Dim >::Metric( simplex );
	};
	_init( simplices , gFunction );
}

template< unsigned int Dim , unsigned int Degree >
template< typename Index >
void SimplexMesh< Dim , Degree >::_init( const std::vector< SimplexIndex< Dim , Index > > &simplices , MetricFunctor< Index > gFunction )
{
	_simplices.resize( simplices.size() );
	_g.resize( _simplices.size() );
#ifdef USE_UNORDERED_SET_MAP
	_nodeMap.reserve( _simplices.size() * NodesPerSimplex );
#endif // USE_UNORDERED_SET_MAP
	for( unsigned int s=0 ; s<_simplices.size() ; s++ )
	{
		for( unsigned int d=0 ; d<=Dim ; d++ ) _simplices[s][d] = (unsigned int)simplices[s][d];
		_g[s] = gFunction(s);
		for( unsigned int n=0 ; n<SimplexElements< Dim , Degree >::NodeNum ; n++ ) _nodeMap[ nodeMultiIndex( s , n ) ] = 0;
	}
	unsigned int nodeCount = 0;
	for( auto iter=_nodeMap.begin() ; iter!=_nodeMap.end() ; iter++ ) iter->second = nodeCount++;
}

template< unsigned int Dim , unsigned int Degree >
typename SimplexMesh< Dim , Degree >::NodeMultiIndex SimplexMesh< Dim , Degree >::nodeMultiIndex( unsigned int s , unsigned int n ) const
{
	unsigned int indices[ Degree ];
	SimplexElements< Dim , Degree >::FactorNodeIndex( n , indices );
	for( unsigned int j=0 ; j<Degree ; j++ ) indices[j] = _simplices[s][ indices[j] ];
	return NodeMultiIndex( indices );
}

template< unsigned int Dim , unsigned int Degree >
unsigned int SimplexMesh< Dim , Degree >::nodeIndex( const NodeMultiIndex &mIdx ) const
{
	auto iter = _nodeMap.find( mIdx );
	if( iter==_nodeMap.end() ) ERROR_OUT( "Could not find node index: " , mIdx );
	return iter->second;
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > SimplexMesh< Dim , Degree >::mass( void ) const
{
	return system( []( const SquareMatrix< double , Dim > &g ){ return SimplexElements< Dim , Degree >::MassMatrix( g ); } );
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > SimplexMesh< Dim , Degree >::stiffness( void ) const
{
	return system( []( const SquareMatrix< double , Dim > &g ){ return SimplexElements< Dim , Degree >::GradientSquareNormMatrix( g ); } );
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > SimplexMesh< Dim , Degree >::bistiffness( void ) const
{
	return system( []( const SquareMatrix< double , Dim > &g ){ return SimplexElements< Dim , Degree >::HessianSquareNormMatrix( g ); } );
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > SimplexMesh< Dim , Degree >::system( std::function< SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum > ( const SquareMatrix< double , Dim > & ) > m2s ) const
{
	std::vector< Eigen::Triplet< double > > entries;
	entries.resize( _simplices.size() * NodesPerSimplex * NodesPerSimplex );
#pragma omp parallel for
	for( int s=0 ; s<(int)_simplices.size() ; s++ )
	{
		unsigned int indices[ SimplexElements< Dim , Degree >::NodeNum ];
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum ; i++ ) indices[i] = nodeIndex( s , i );

		SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum > M = m2s( _g[s] );
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum ; i++ )
			for( unsigned int j=0 ; j<SimplexElements< Dim , Degree >::NodeNum ; j++ )
				entries[ s*NodesPerSimplex*NodesPerSimplex + i * NodesPerSimplex + j ] = Eigen::Triplet< double >( indices[i] , indices[j] , M(i,j) );
	}
	Eigen::SparseMatrix< double > M( nodes() , nodes() );
	M.setFromTriplets( entries.begin() , entries.end() );
	return M;
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > SimplexMesh< Dim , Degree >::crossFaceGradientEnergy( void ) const
{
	return crossFaceGradientEnergy( []( FaceMultiIndex ){ return true; } );
}

template< unsigned int Dim , unsigned int Degree >
template< typename UseFaceFunctor >
Eigen::SparseMatrix< double > SimplexMesh< Dim , Degree >::crossFaceGradientEnergy( const UseFaceFunctor &useFaceFunctor ) const
{
	std::vector< Eigen::Triplet< double > > entries;
	// Ballpark, there are:
	//	(Dim+1) faces on each simplex
	//	two simplices incident on each face,
	//	NodesPerSimplex nodes associated with each simplex (some of which will be double-counted)
	entries.reserve( _simplices.size() * (Dim+1) * ( 2 * NodesPerSimplex ) * ( 2 * NodesPerSimplex ) );

	// Functionality for returning the multi-index associated with a face of a simplex
	auto GetFaceMultiIndex = [&]( unsigned int s , unsigned int f , Permutation< Dim > &p )
	{
		SimplexIndex< Dim-1 , unsigned int > faceIndex = RightSimplex< Dim >::Face( f );
		p = Permutation< Dim >( [&]( unsigned int i , unsigned int j ){ return _simplices[s][ faceIndex[i] ]<_simplices[s][ faceIndex[j] ]; } );
		for( unsigned int d=0 ; d<Dim ; d++ ) faceIndex[d] = _simplices[s][ faceIndex[d] ];
		return FaceMultiIndex( &faceIndex[0] );
	};

	// A mapping from simplex faces to indices
	typename FaceMultiIndex::map faceMap;
	{
		Permutation< Dim > p;
		for( unsigned int s=0 ; s<_simplices.size() ; s++ )	for( unsigned int f=0 ; f<=Dim ; f++ ) faceMap[ GetFaceMultiIndex( s , f , p ) ] = 0;

		unsigned int faceCount = 0;
		for( auto & [ mIdx , idx ] : faceMap ) idx = faceCount++;
	}

	// The per face gradient components of the node functions
#ifdef USE_UNORDERED_SET_MAP
	std::vector< std::unordered_map< unsigned int , Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > > > faceGradientComponents( faceMap.size() );
#else // !USE_UNORDERED_SET_MAP
	std::vector< std::map< unsigned int , Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > > > faceGradientComponents( faceMap.size() );
#endif // USE_UNORDERED_SET_MAP
	// The per face metrics
	std::vector< SquareMatrix< double , Dim-1 > > faceMetrics( faceMap.size() );

	for( unsigned int s=0 ; s<_simplices.size() ; s++ )
	{
		// Compute the gradient components and metrics for the faces
		Matrix< Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > _faceGradientComponents;
		Point< SquareMatrix< double , Dim-1 > , Dim+1 > _faceMetrics;

		_faceGradientComponents = SimplexElements< Dim , Degree >::FaceGradientOrthogonalComponents( _g[s] , &_simplices[s][0] );
		_faceMetrics = SimplexElements< Dim , Degree >::FaceMetrics( _g[s] , &_simplices[s][0] );

		unsigned int indices[ SimplexElements< Dim , Degree >::NodeNum ];
		for( unsigned int i=0 ; i<SimplexElements< Dim , Degree >::NodeNum ; i++ ) indices[i] = nodeIndex( s , i );

		for( unsigned int f=0 ; f<=Dim ; f++ )
		{
			Permutation< Dim > p;
			FaceMultiIndex fmi = GetFaceMultiIndex( s , f , p );
			if( useFaceFunctor(fmi) )
			{
				bool evenParity = ( p.parity()&1 )==0;
				unsigned int fi = faceMap[fmi];
				Matrix< double , Dim , Dim-1 > A = RightSimplex< Dim-1 >::AffineTransform( p );
				// Accumulate the metrices and (signed) gradient components
				faceMetrics[fi] += _faceMetrics[f]/2;
				for( unsigned int n=0 ; n<SimplexElements< Dim , Degree >::NodeNum ; n++ )
				{
					Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > &gComponents = faceGradientComponents[fi][ indices[n] ];
					for( unsigned int d=0 ; d<Dim ; d++ )
						if( evenParity ) gComponents[d] += _faceGradientComponents(n,f)[d].template operator()< Dim >(A);
						else             gComponents[d] -= _faceGradientComponents(n,f)[d].template operator()< Dim >(A);
				}
			}
		}
	}

	for( unsigned int f=0 ; f<faceGradientComponents.size() ; f++ )
		for( const auto & [ i1 , g1 ]: faceGradientComponents[f] ) for( const auto & [ i2 , g2 ] : faceGradientComponents[f] )
		{
			double value = 0;
			for( unsigned int d=0 ; d<Dim ; d++ ) value += RightSimplex< Dim-1 >::Integral( g1[d] * g2[d] , faceMetrics[f] );
			entries.push_back( Eigen::Triplet< double >( i1 , i2 , value ) );
		}

	Eigen::SparseMatrix< double > E( nodes() , nodes() );
	E.setFromTriplets( entries.begin() , entries.end() );
	return E;
}

template< unsigned int Dim , unsigned int Degree >
template< unsigned int FaceDim >
Eigen::SparseMatrix< double > SimplexMesh< Dim , Degree >::_crossFaceGradientEnergy( void ) const
{
	static_assert( FaceDim<Dim , "[ERROR] Face dimension must be smaller than simplex dimension" );
	std::vector< Eigen::Triplet< double > > entries;
	// Ballpark, there are:
	//	Choose(Dim+1,FaceDim) faces on each simplex
	//	two simplices incident on each face,
	//	NodesPerSimplex nodes associated with each simplex (some of which will be double-counted)
	entries.reserve( _simplices.size() * (Dim+1) * ( 2 * NodesPerSimplex ) * ( 2 * NodesPerSimplex ) );

	// Functionality for returning the multi-index associated with a face of a simplex
	auto GetFaceMultiIndex = [&]( unsigned int s , unsigned int f , Permutation< Dim > &p )
	{
		SimplexIndex< Dim-1 , unsigned int > faceIndex = RightSimplex< Dim >::Face( f );
		p = Permutation< Dim >( [&]( unsigned int i , unsigned int j ){ return _simplices[s][ faceIndex[i] ]<_simplices[s][ faceIndex[j] ]; } );
		for( unsigned int d=0 ; d<Dim ; d++ ) faceIndex[d] = _simplices[s][ faceIndex[d] ];
		return FaceMultiIndex( &faceIndex[0] );
	};

	// A mapping from simplex faces to indices
	typename FaceMultiIndex::map faceMap;
	{
		Permutation< Dim > p;
		for( unsigned int s=0 ; s<_simplices.size() ; s++ )	for( unsigned int f=0 ; f<=Dim ; f++ ) faceMap[ GetFaceMultiIndex( s , f , p ) ] = 0;

		unsigned int faceCount = 0;
		for( auto & [ mIdx , idx ] : faceMap ) idx = faceCount++;
	}

	// The per face gradient components of the node functions
#ifdef USE_UNORDERED_SET_MAP
	std::vector< std::unordered_map< unsigned int , Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > > > faceGradientComponents( faceMap.size() );
#else // !USE_UNORDERED_SET_MAP
	std::vector< std::map< unsigned int , Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > > > faceGradientComponents( faceMap.size() );
#endif // USE_UNORDERED_SET_MAP
	// The per face metrics
	std::vector< SquareMatrix< double , Dim-1 > > faceMetrics( faceMap.size() );

	for( unsigned int s=0 ; s<_simplices.size() ; s++ )
	{
		// Compute the gradient components and metrics for the faces
		Matrix< Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > _faceGradientComponents;
		Point< SquareMatrix< double , Dim-1 > , Dim+1 > _faceMetrics;

		_faceGradientComponents = SimplexElements< Dim , Degree >::FaceGradientComponents( _g[s] , &_simplices[s][0] );
		_faceMetrics = SimplexElements< Dim , Degree >::FaceMetrics( _g[s] , &_simplices[s][0] );

		for( unsigned int f=0 ; f<=Dim ; f++ )
		{
			Permutation< Dim > p;
			FaceMultiIndex fmi = GetFaceMultiIndex( s , f , p );
			if( useFaceFunctor(fmi) )
			{
				bool evenParity = ( p.parity()&1 )==0;
				unsigned int fi = faceMap[fmi];
				Matrix< double , Dim , Dim-1 > A = RightSimplex< Dim-1 >::AffineTransform( p );
				// Accumulate the metrices and (signed) gradient components
				faceMetrics[fi] += _faceMetrics[f]/2;
				for( unsigned int n=0 ; n<SimplexElements< Dim , Degree >::NodeNum ; n++ )
				{
					Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > &gComponents = faceGradientComponents[fi][ nodeIndex( s , n ) ];
					for( unsigned int d=0 ; d<Dim ; d++ )
						if( evenParity ) gComponents[d] += _faceGradientComponents(n,f)[d].template operator()< Dim >(A);
						else             gComponents[d] -= _faceGradientComponents(n,f)[d].template operator()< Dim >(A);
				}
			}
		}
	}

	for( unsigned int f=0 ; f<faceGradientComponents.size() ; f++ )
		for( const auto & [ i1 , g1 ]: faceGradientComponents[f] ) for( const auto & [ i2 , g2 ] : faceGradientComponents[f] )
		{
			double value = 0;
			for( unsigned int d=0 ; d<Dim ; d++ ) value += RightSimplex< Dim-1 >::Integral( g1[d] * g2[d] , faceMetrics[f] );
			entries.push_back( Eigen::Triplet< double >( i1 , i2 , value ) );
		}

	Eigen::SparseMatrix< double > E( nodes() , nodes() );
	E.setFromTriplets( entries.begin() , entries.end() );
	return E;
}

template< unsigned int Dim , unsigned int Degree >
void SimplexMesh< Dim , Degree >::hashLocalToGlobalNodeIndex( void )
{
	_localToGlobalNodeIndex.resize( simplices() * NodesPerSimplex );
#pragma omp parallel for
	for( int s=0 ; s<(int)simplices() ; s++ ) for( unsigned int n=0 , i=s*NodesPerSimplex ; n<NodesPerSimplex ; n++ , i++ )
		_localToGlobalNodeIndex[i] = nodeIndex( nodeMultiIndex( s , n ) );
}

template< unsigned int Dim , unsigned int Degree >
Eigen::SparseMatrix< double > SimplexMesh< Dim , Degree >::evaluationMatrix( const std::vector< typename SimplexMesh< Dim >::Sample > &samples ) const
{
	std::vector< Eigen::Triplet< double > > entries( samples.size()*NodesPerSimplex );

	Polynomial::Polynomial< Dim , Degree , double > elements[ NodesPerSimplex ];
	SimplexElements< Dim , Degree >::SetElements( elements );

#pragma omp parallel for
	for( int i=0 ; i<(int)samples.size() ; i++ )
	{
		unsigned int e = i * NodesPerSimplex;
		Point< double , Dim > p;
		for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = samples[i].bcCoordinates[d+1];
		for( unsigned int n=0 ; n<NodesPerSimplex ; n++ )
		{
			double v = elements[n]( p );
			unsigned int idx = nodeIndex( samples[i].sIdx , n );
			entries[e++] = Eigen::Triplet< double >( i , idx , v );
		}
	}
	Eigen::SparseMatrix< double > E( samples.size() , nodes() );
	E.setFromTriplets( entries.begin() , entries.end() );
	return E;
}

template< unsigned int Dim , unsigned int Degree >
double SimplexMesh< Dim , Degree >::evaluate( const Eigen::VectorXd &coefficients , typename SimplexMesh< Dim >::Sample sample ) const
{
	double value = 0;
	static Polynomial::Polynomial< Dim , Degree , double > elements[ NodesPerSimplex ];
	static bool firstTime = true;
	if( firstTime )
	{
		SimplexElements< Dim , Degree >::SetElements( elements );
		firstTime = false;
	}

	Point< double , Dim > p;
	for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = sample.bcCoordinates[d+1];
	for( unsigned int n=0 ; n<NodesPerSimplex ; n++ ) value += elements[n]( p ) * coefficients[ nodeIndex( sample.sIdx , n ) ];

	return value;
}

template< unsigned int Dim , unsigned int Degree >
Polynomial::Polynomial< Dim , Degree , double > SimplexMesh< Dim , Degree >::evaluate( const Eigen::VectorXd &coefficients , unsigned int simplexIndex ) const
{
	static Polynomial::Polynomial< Dim , Degree , double > elements[ NodesPerSimplex ];
	static bool firstTime = true;
	if( firstTime )
	{
		SimplexElements< Dim , Degree >::SetElements( elements );
		firstTime = false;
	}

	Polynomial::Polynomial< Dim , Degree , double > poly;
	for( unsigned int n=0 ; n<NodesPerSimplex ; n++ ) poly += elements[n] * coefficients[ nodeIndex( simplexIndex , n ) ];
	return poly;
}

template< unsigned int Dim , unsigned int Degree >
double SimplexMesh< Dim , Degree >::volume( void ) const
{
	double v = 0;
	for( unsigned int i=0 ; i<_g.size() ; i++ ) v += sqrt( fabs( _g[i].determinant() ) );
	for( unsigned int i=2 ; i<=Dim ; i++ ) v /= i;
	return v;
}

template< unsigned int Dim , unsigned int Degree >
void SimplexMesh< Dim , Degree >::makeUnitVolume( void )
{
	double scale = pow( volume() , -2./Dim );
	for( unsigned int i=0 ; i<_g.size() ; i++ ) _g[i] *= scale;
}

template< unsigned int Dim , unsigned int Degree >
void SimplexMesh< Dim , Degree >::setMetric( std::function< SquareMatrix< double , Dim > (unsigned int) > metricFunction )
{
	for( unsigned int i=0 ; i<_g.size() ; i++ ) _g[i] = metricFunction(i);
}
