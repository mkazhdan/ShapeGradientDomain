#ifndef GRADIENT_DOMAIN_INCLUDED
#define GRADIENT_DOMAIN_INCLUDED

#include <Eigen/Dense>
#include <Misha/MultiThreading.h>
#include <Misha/Geometry.h>
#include <Misha/FEM.h>

namespace GradientDomain
{
	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexVertex( const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyVertexFunctor && High );

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyEdgeFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexEdge( const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyEdgeFunctor && High );

	////////////////////
	// Implementation //
	////////////////////

	template< typename T >
	struct Array{ static const unsigned int Dim=0; };

	template<>
	struct Array< float >
	{
		using T = float;
		using Real = float;
		static const unsigned int Dim=1;
		static const Real & Entry( const T &v , unsigned int ){ return v; }
		static       Real & Entry(       T &v , unsigned int ){ return v; }
		template< unsigned int D > static constexpr const Real & Entry( const T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return v; }
		template< unsigned int D > static constexpr       Real & Entry(       T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return v; }
	};

	template<>
	struct Array< double >
	{
		using T = double;
		using Real = double;
		static const unsigned int Dim=1;
		static const Real & Entry( const T &v , unsigned int ){ return v; }
		static       Real & Entry(       T &v , unsigned int ){ return v; }
		template< unsigned int D > static constexpr const Real & Entry( const T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return v; }
		template< unsigned int D > static constexpr       Real & Entry(       T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return v; }
	};

	template< typename _T , unsigned int _Dim , typename _Real >
	struct Array< Point< _T , _Dim , _Real > >
	{
		using T = Point< _T , _Dim , _Real >;
		using Real = _Real;
		static const unsigned int Dim = _Dim * Array< _T >::Dim;
		static const Real & Entry( const T &v , unsigned int d ){ return Array< _T >::Entry( v[ d%_Dim ] , d/_Dim ); }
		static       Real & Entry(       T &v , unsigned int d ){ return Array< _T >::Entry( v[ d%_Dim ] , d/_Dim ); }
		template< unsigned int D > static constexpr const Real & Entry( const T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return Array< _T >::template Entry< D/_Dim >( v[ D%_Dim ] ); }
		template< unsigned int D > static constexpr       Real & Entry(       T &v ){ static_assert( D<=Dim , "[ERROR] Bad index" ) ; return Array< _T >::template Entry< D/_Dim >( v[ D%_Dim ] ); }
	};

	template< unsigned int D , typename T , typename Real , typename EigenSolver , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ >
	void _Solve( const EigenSolver &solver , const Eigen::SparseMatrix< Real > &lowM , const Eigen::SparseMatrix< Real > &highM , LowFrequencyVertexFunctor && Low , HighFrequencyVertexFunctor && High , std::vector< T > &out  )
	{
		if constexpr( D < Array< T >::Dim )
		{
			Eigen::Matrix< Real , Eigen::Dynamic , 1 > l , h;
			l.resize( lowM.cols() ) , h.resize( highM.cols() );
			for( unsigned int i=0 ; i< lowM.cols() ; i++ ) l[i] = Array< T >::template Entry< D >(  Low(i) );
			for( unsigned int i=0 ; i<highM.cols() ; i++ ) h[i] = Array< T >::template Entry< D >( High(i) );

#ifdef EIGEN_USE_MKL_ALL
			Eigen::Matrix< Real , Eigen::Dynamic , 1 > b = lowM * l + highM * h;
			Eigen::Matrix< Real , Eigen::Dynamic , 1 > x = solver.solve( b );
#else // !EIGEN_USE_MKL_ALL
			Eigen::Matrix< Real , Eigen::Dynamic , 1 > x = solver.solve( lowM * l + highM * h );
#endif // EIGEN_USE_MKL_ALL

			for( unsigned int i=0 ; i<x.size() ; i++ ) Array< T >::template Entry< D >( out[i] ) = x[i];
			_Solve< D+1 >( solver , lowM , highM , std::forward< LowFrequencyVertexFunctor >( Low ) , std::forward< HighFrequencyVertexFunctor >( High ) , out );
		}
	}

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexVertex( const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyVertexFunctor && High )
	{
		static_assert( Array< T >::Dim , "[ERROR] T is not of array type" );
		static_assert( std::is_convertible_v< LowFrequencyVertexFunctor , std::function< T ( unsigned int ) > > , "[ERROR] LowFrequencyVertexFunctor is poorly formed" );
		static_assert( std::is_convertible_v< HighFrequencyVertexFunctor , std::function< T ( unsigned int ) > > , "[ERROR] HighFrequencyVertexFunctor is poorly formed" );

		Eigen::SparseMatrix< Real > mass = mesh.template massMatrix< FEM::BASIS_0_WHITNEY , true >() * lowWeight;
		Eigen::SparseMatrix< Real > stiffness = mesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY , true >() * highWeight;

		EigenSolver solver( mass + stiffness );

		std::vector< T > out( mesh.vCount() );
		_Solve< 0 >( solver , mass , stiffness , std::forward< LowFrequencyVertexFunctor >( Low ) , std::forward< HighFrequencyVertexFunctor >( High ) , out );
		return out;
	}

	template< typename EigenSolver , typename T , typename Real , typename LowFrequencyVertexFunctor /* = std::function< T ( unsigned int ) > */ , typename HighFrequencyEdgeFunctor /* = std::function< T ( unsigned int ) > */ >
	std::vector< T > ProcessVertexEdge( const FEM::RiemannianMesh< Real > &mesh , Real lowWeight , Real highWeight , LowFrequencyVertexFunctor && Low , HighFrequencyEdgeFunctor && High )
	{
		static_assert( Array< T >::Dim , "[ERROR] T is not of array type" );
		static_assert( std::is_convertible_v< LowFrequencyVertexFunctor , std::function< T ( unsigned int ) > > , "[ERROR] LowFrequencyVertexFunctor is poorly formed" );
		static_assert( std::is_convertible_v< HighFrequencyEdgeFunctor , std::function< T ( unsigned int ) > > , "[ERROR] HighFrequencyEdgeFunctor is poorly formed" );

		Eigen::SparseMatrix< Real > mass = mesh.template massMatrix< FEM::BASIS_0_WHITNEY , true >() * lowWeight;
		Eigen::SparseMatrix< Real > d = mesh.template dMatrix< FEM::BASIS_0_WHITNEY , FEM::BASIS_1_WHITNEY , true >();
		Eigen::SparseMatrix< Real > divergence = d.transpose() * mesh.template massMatrix< FEM::BASIS_1_WHITNEY , true >() * highWeight;

		EigenSolver solver( mass + divergence * d );

		std::vector< T > out( mesh.vCount() );
		_Solve< 0 >( solver , mass , divergence , std::forward< LowFrequencyVertexFunctor >( Low ) , std::forward< HighFrequencyEdgeFunctor >( High ) , out );
		return out;
	}

};
#endif // GRADIENT_DOMAIN_INCLUDED