#ifndef TIMER_INCLUDED
#define TIMER_INCLUDED
#include <sys/timeb.h>
#if !defined( WIN32 ) && !defined( _WIN64 )
#include <sys/time.h>
#endif // !WIN32 && !_WIN64
#if defined( WIN32 ) || defined( _WIN64 )
struct Timer
{
	static double Time( void )
	{
		struct _timeb t;
		_ftime( &t );
		return double(t.time)+double(t.millitm)/1000.0;
	}
	struct _timeb t;
	Timer( void ){ _ftime( &t ); }
	void reset( void ){ _ftime( &t ); }
	double elapsed( void ) const
	{
		struct _timeb _t;
		_ftime( &_t );
		return (double)(_t.time-t.time)+(double)(_t.millitm-t.millitm)/1000.;
	}
};
#else // !WIN32 && !_WIN64
struct Timer
{
	static double Time( void )
	{
		struct timeval t;
		gettimeofday(&t,NULL);
		return t.tv_sec+(double)t.tv_usec/1000000;
	}
	struct timeval t;
	Timer( void ){ gettimeofday( &t , NULL ); }
	void reset( void ){ gettimeofday( &t , NULL ); }
	double elapsed( void ) const
	{
		struct timeval _t;
		gettimeofday( &_t , NULL );
		return (double)(_t.tv_sec-t.tv_sec)+(double)(_t.tv_usec-t.tv_usec)/1000000;
	}
};
#endif // WIN32 || WIN64
#endif // TIMER_INCLUDED