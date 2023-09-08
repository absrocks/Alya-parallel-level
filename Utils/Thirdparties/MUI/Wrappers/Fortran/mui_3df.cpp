/*
 * mui_3df.cpp
 *
 *  Created on: Jan 21, 2015
 *      Author: skudo
 */
#include "../mui.h"
#include "../string_f2c.h" 

using namespace mui;

static uniface3d* ptr_class = NULL;

extern "C" {


typedef uniface3d                        mui_uniface3d;
typedef sampler_gauss3d<double>          mui_sampler_gauss3d;
typedef sampler_moving_average3d<double> mui_sampler_moving_average3d;
typedef chrono_sampler_exact3d           mui_chrono_sampler_exact3d;
typedef chrono_sampler_mean3d            mui_chrono_sampler_mean3d;


void mui_create_uniface3d_f_( const char *ftype, int* fcomm) 
{
  MPI_Comm  new_world_comm = MPI_COMM_NULL; 
  new_world_comm = mui::mpi_split_by_app();  
  fcomm[0] = MPI_Comm_c2f( new_world_comm );

  std::string  ctype(ftype); 
//std::string  URI = "mpi://"+ ctype +"/interface";
  std::string  URI = "mpi://"+ ctype +"/ifs";

  ptr_class = new mui_uniface3d( URI );
}


void mui_destroy_uniface3d_f_() 
{
  delete ptr_class;
}


void mui_commit_f_( double* t ) 
{
  std::string sent = ( ptr_class->commit(t[0]) )?("ON"):("OFF");
  std::cout<<"sending:"<< sent <<"\n";
}


void mui_push_f_( const char *attr, double* x, double* y, double* z, double* value )
{
  ptr_class->push( std::string(attr), point3d(*x,*y,*z), *value );
}

// spatial sampler: moving_average
// temporal sampler: mean
void mui_fetch_moving_average_mean_f_(const char *attr, double* x, double* y, double* z, double* t, double* value)
{
/*
  mui_sampler_moving_average3d  *spatial
  mui_chrono_sampler_mean3d     *temporal 
  ptr_class->fetch( std::string(attr), point3d(*x,*y,*z), *t, *spatial, *temporal );
*/
  mui::sampler_gauss3d<double>  spatial( 1.0, 1.0/2 );
  mui::chrono_sampler_exact1d   temporal;
  value[0] = ptr_class->fetch( 
                              std::string(attr), 
                              mui::point3d(x[0],y[0],z[0]), t[0], 
                              spatial, temporal 
                             );
  
/*
  for(int i=0,k=0; i<vdata.size(); i++)
  {
    mui::point3d loc( vdata[i][0], vdata[i][1], vdata[i][2] );
    auto ux = interface.fetch( "ux", loc, 0.0, gauss, exact );
  }
*/ 
}

/*
void mui_create_sampler_3d_f_( mui_sampler_gauss3d** ret, double* r, double* h ) {
	*ret = new mui_sampler_gauss3d( *r, *h );
}

void mui_create_sampler_moving_average3d_f_( mui_sampler_moving_average3d** ret, double* dx, double* dy, double* dz ) {
	*ret = new mui_sampler_moving_average3d( point3d(*dx,*dy,*dz) );
}

void mui_create_chrono_sampler_exact3d_f_(mui_chrono_sampler_exact3d** ret) {
	*ret = new mui_chrono_sampler_exact3d;
}

void mui_create_chrono_sampler_mean3d_f_( mui_chrono_sampler_mean3d** ret, double* past, double* future ) {
	*ret = new mui_chrono_sampler_mean3d( *past, *future );
}

// deallocator
void mui_destroy_uniface3d_f_( mui_uniface3d *uniface ) {
	delete uniface;
}

void mui_destroy_sampler_3d_f_( mui_sampler_gauss3d* sampler ) {
	delete sampler;
}

void mui_destroy_sampler_moving_average3d_f_( mui_sampler_moving_average3d* sampler ) {
	delete sampler;
}

void mui_destroy_chrono_sampler_exact3d_f_( mui_chrono_sampler_exact3d* sampler ) {
	delete sampler;
}

void mui_destroy_chrono_sampler_mean3d_f_( mui_chrono_sampler_mean3d* sampler ) {
	delete sampler;
}

// push
void mui_push_f_( mui_uniface3d* uniface, const char *attr, double* x, double* y, double* z, double* value ) {
	uniface->push( std::string(attr), point3d(*x,*y,*z), *value );
}

// spatial sampler: gaussian
// temporal sampler: exact point
double mui_fetch_gaussian_exact_f_( mui_uniface3d* uniface, const char *attr, double* x, double* y, double* z, double* t, mui_sampler_gauss3d *spatial, mui_chrono_sampler_exact3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(*x,*y,*z), *t, *spatial, *temporal );
}

// spatial sampler: gaussian
// temporal sampler: mean
double mui_fetch_gaussian_mean_f_( mui_uniface3d* uniface, const char *attr, double* x, double* y, double* z, double* t, mui_sampler_gauss3d *spatial, mui_chrono_sampler_mean3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(*x,*y,*z), *t, *spatial, *temporal );
}

// spatial sampler: moving average
// temporal sampler: exact point
double mui_fetch_moving_average_exact_f_( mui_uniface3d* uniface, const char *attr, double* x, double* y, double* z, double* t, mui_sampler_moving_average3d *spatial, mui_chrono_sampler_exact3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(*x,*y,*z), *t, *spatial, *temporal );
}

// spatial sampler: moving_average
// temporal sampler: mean
double mui_fetch_moving_average_mean_f_( mui_uniface3d* uniface, const char *attr, double* x, double* y, double* z, double* t, mui_sampler_moving_average3d *spatial, mui_chrono_sampler_mean3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(*x,*y,*z), *t, *spatial, *temporal );
}

// commit all data in buffer
void mui_commit_f_( mui_uniface3d* uniface, double* t ) {
	uniface->commit( *t );
}

// wait for peers
void mui_barrier_f_( mui_uniface3d* uniface, double* t ) {
	uniface->barrier( *t );
}

// remove obsolete data
void mui_forget_f_( mui_uniface3d* uniface, double* first, double* last ) {
	uniface->forget( *first, *last );
}

// set automatic deletion
void mui_set_memory_f_( mui_uniface3d* uniface, double* length ) {
	return uniface->set_memory( *length );
}

*/ 

}
