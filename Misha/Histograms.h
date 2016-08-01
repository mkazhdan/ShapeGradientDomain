#ifndef HISTOGRAMS_INCLUDED
#define HISTOGRAMS_INCLUDED
#include <stdio.h>
#include <Misha/Array.h>

namespace Histogram
{
	template< class Real > struct Sample{ Real weight , value; };
	template< class Real > inline void Print( ConstPointer( Real ) histogram , int bins , int units , Real min , Real max );
	template< class Real > inline void Print( const std::vector<         Real   >& values , int bins , int units , Real xMin , Real xMax );
	template< class Real > inline void Print( const std::vector< Sample< Real > >& values , int bins , int units , Real xMin , Real xMax );
	template< class Real > inline void Print( const std::vector< Real >& values , int bins , int units , bool symmetric=false );
	template< class Real > inline void Print( const std::vector< Sample< Real > >& values , int bins , int units , bool symmetric=false );

	template< class Real >
	inline void Print( ConstPointer( Real ) histogram , int bins , int units , Real xMin , Real xMax )
	{
		Real yMax = 0;
		for( int i=0 ; i<bins ; i++ ) yMax = std::max< Real >( yMax , histogram[i] );
		printf( "[ %f , %f ] x [ 0 , %f ]\n" , xMin , xMax , yMax );
		for( int i=0 ; i<units ; i++ ) 
		{
			Real cutOff = (Real) ( yMax / units * ( units - i ) );
			for( int j=0 ; j<bins ; j++ )
			{
				if( histogram[j]>=cutOff ) printf( " " );
				else                       printf( "*" );
			}
			printf( "\n" );
		}
		for( int j=0 ; j<bins ; j++ )
		{
			Real start = xMin + ( xMax - xMin ) / bins * j , end = xMin + ( xMax - xMin ) / bins * ( j+1 );
			if( start<0 && end>=0 ) printf( "0" );
			else                    printf( "-" );
		}
		printf( "\n" );
	}
	template< class Real >
	inline void Print( const std::vector< Real >& values , int bins , int units , Real xMin , Real xMax )
	{
		Pointer( Real ) histogram = AllocPointer< Real >( bins );
		memset( histogram , 0 , sizeof(Real)* bins );
		for( int i=0 ; i<values.size() ; i++ )
		{
			Real b = ( values[i] - xMin ) / ( xMax - xMin ) * ( bins-1 );
			int ib1 = (int)floor( b ) , ib2 = ib1+1;
			Real d1 = (Real)( 1. - (b-ib1) ) , d2 = (Real)( b - ib1 );
			if( ib1>=0 && ib1<bins ) histogram[ib1] += d1;
			if( ib2>=0 && ib2<bins ) histogram[ib2] += d2;
		}
		for( int i=0 ; i<bins ; i++ ) histogram[i] /= values.size();
		Print( histogram , bins , units , xMin , xMax );
		FreePointer( histogram );
	}
	template< class Real >
	inline void Print( const std::vector< Sample< Real > >& values , int bins , int units , Real xMin , Real xMax )
	{
		Real* histogram = AllocPointer< Real >( bins );
		memset( histogram , 0 , sizeof(Real)* bins );
		Real weightSum = 0;
		for( int i=0 ; i<values.size() ; i++ )
		{
			Real b = ( values[i].value - xMin ) / ( xMax - xMin ) * ( bins-1 );
			int ib1 = (int)floor( b ) , ib2 = ib1+1;
			Real d1 = (Real)( 1. - (b-ib1) ) * values[i].weight , d2 = (Real)( b - ib1 ) * values[i].weight;
			weightSum += values[i].weight;
			if( ib1>=0 && ib1<bins ) histogram[ib1] += d1;
			if( ib2>=0 && ib2<bins ) histogram[ib2] += d2;
		}
		for( int i=0 ; i<bins ; i++ ) histogram[i] /= weightSum;
		Print( histogram , bins , units , xMin , xMax );
		FreePointer( histogram );
	}
	template< class Real >
	inline void Print( const std::vector< Real >& values , int bins , int units , bool symmetric )
	{
		Real xMin , xMax;
		xMin = xMax = values[0];
		for( int i=0 ; i<values.size() ; i++ ) xMin = std::min< Real >( xMin , values[i] ) , xMax = std::max< Real >( xMax , values[i] );
		if( symmetric ) xMin = std::min< Real >( xMin , -xMax ) , xMax = std::max< Real >( xMax , -xMin );
		Print( values , bins , units , xMin , xMax );
	}
	template< class Real >
	inline void Print( const std::vector< Sample< Real > >& values , int bins , int units , bool symmetric )
	{
		Real xMin , xMax;
		xMin = xMax = values[0];
		for( int i=0 ; i<values.size() ; i++ ) xMin = std::min< Real >( xMin , values[i].value ) , xMax = std::max< Real >( xMax , values[i].value );
		if( symmetric ) min = std::min< Real >( xMin , -xMax ) , xMax = std::max< Real >( xMax , -xMin );
		Print( values , bins , units , xMin , xMax );
	}
}

#endif // HISTOGRAMS_INCLUDED