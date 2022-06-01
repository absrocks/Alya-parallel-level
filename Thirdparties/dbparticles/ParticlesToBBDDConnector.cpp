#include <iostream>
#include <fstream>
#include <ctime>

#include <thrift/protocol/TBinaryProtocol.h>
#include <thrift/transport/TSocket.h>
#include <thrift/transport/TTransportUtils.h>

#include "ParticlesToBBDD.h"

using namespace std;
using namespace apache::thrift;
using namespace apache::thrift::protocol;
using namespace apache::thrift::transport;

using namespace BBDDParticles;

extern "C" 
{
	void sendparticles(int npart, int *recvcounts, int nvar2, int kfl_posla_pts, double cutim, double *recvbuf_rp);
	void timestamp (int ll);
}

void sendparticles(int npart, int *recvcounts, int nvar2, int kfl_posla_pts, double cutim, double *recvbuf_rp)
{
	ofstream myfile;
	int ipars = -1;
	int jpart = 0;
	InsertData data;
	InsertResult result;
   //cc[ll--] = '\0';  // NULL terminate the string

   myfile.open ("connector.txt", std::ios_base::app );
   printf("From doubleIJK: %d\n",npart);
   
  boost::shared_ptr<TTransport> socket(new TSocket("localhost", 9090));
  boost::shared_ptr<TTransport> transport(new TBufferedTransport(socket));
  boost::shared_ptr<TProtocol> protocol(new TBinaryProtocol(transport));
  ParticlesToBBDDClient client(protocol);

  try {
	/*Recorrido correcto de las particulas y obtención de sus propiedades básico
	for( int ipart = 2; ipart <= npart+1; ipart = ipart + 1 ) {
		jpart = ipart - 1;
		for( int ilagr = 1; ilagr <= recvcounts[ipart - 1] / nvar2; ilagr = ilagr + 1 ) {
			if( kfl_posla_pts == 0 ) {
				myfile << cutim << "," << abs(int(recvbuf_rp[ipars+10]))
				<< "," << recvbuf_rp[ipars+1] << "," << recvbuf_rp[ipars+2] 
				<< "," << recvbuf_rp[ipars+3] << "," << recvbuf_rp[ipars+4]
				<< "," << recvbuf_rp[ipars+5] << "," << recvbuf_rp[ipars+6]
				<< "," << int(recvbuf_rp[ipars+11]) << "," << jpart
				<< endl;
			}
			else if (kfl_posla_pts == 1) {
				myfile << cutim << "," << abs(int(recvbuf_rp[ipars+10]))
				<< "," << recvbuf_rp[ipars+1] << "," << recvbuf_rp[ipars+2] 
				<< "," << recvbuf_rp[ipars+3] << "," << recvbuf_rp[ipars+4]
				<< "," << recvbuf_rp[ipars+5] << "," << recvbuf_rp[ipars+6]
				<< "," << int(recvbuf_rp[ipars+11]) << "," << jpart
				<< "," << recvbuf_rp[ipars+7] << "," << recvbuf_rp[ipars+8]
				<< "," << recvbuf_rp[ipars+9] << "," << recvbuf_rp[ipars+15]
				<< "," << recvbuf_rp[ipars+16] << endl;
			}
			else {
				myfile << "not programed kfl_posla_pts: " << kfl_posla_pts << endl;
			}
			ipars = ipars + nvar2;
		}
	}					
	myfile.close();
	*/
	//Go over all the particles and create the list to send to de database
	
	for( int ipart = 2; ipart <= npart+1; ipart = ipart + 1 ) {
		jpart = ipart - 1;
		for( int ilagr = 1; ilagr <= recvcounts[ipart - 1] / nvar2; ilagr = ilagr + 1 ) {
			if( kfl_posla_pts == 0 ) {
				myfile << cutim << "," << abs(int(recvbuf_rp[ipars+10]))
				<< "," << recvbuf_rp[ipars+1] << "," << recvbuf_rp[ipars+2] 
				<< "," << recvbuf_rp[ipars+3] << "," << recvbuf_rp[ipars+4]
				<< "," << recvbuf_rp[ipars+5] << "," << recvbuf_rp[ipars+6]
				<< "," << int(recvbuf_rp[ipars+11]) << "," << jpart
				<< endl;
			}
			else if (kfl_posla_pts == 1) {
				Particle particle;
				std::map<std::string, double> mapDouble;
				std::map<std::string, int32_t> mapInt;
				
				particle.time = cutim;

				particle.partId = (abs(int(recvbuf_rp[ipars+10])));

				particle.x = recvbuf_rp[ipars+1];

				particle.y = recvbuf_rp[ipars+2];

				particle.z = recvbuf_rp[ipars+3];
								
                mapDouble["vx"] = recvbuf_rp[ipars+4];
                mapDouble["vy"] = recvbuf_rp[ipars+5];
                mapDouble["vz"] = recvbuf_rp[ipars+6];
                mapInt["Type number"] = int(recvbuf_rp[ipars+11]);
                mapInt["Subdomain"] = jpart;
                mapDouble["ax"] = recvbuf_rp[ipars+7];
                mapDouble["ay"] = recvbuf_rp[ipars+8];
                mapDouble["az"] = recvbuf_rp[ipars+9];
                mapDouble["Cd"] = recvbuf_rp[ipars+15];
                mapDouble["Re"] = recvbuf_rp[ipars+16];
                
                particle.doubleProperties = mapDouble;
                particle.intProperties = mapInt;
                
                data.data.push_back (particle);
			}
			else {
				myfile << "not programed kfl_posla_pts: " << kfl_posla_pts << endl;
			}
			ipars = ipars + nvar2;
		}
	}
					
    transport->open();
    client.ping();
    client.insertParticles(result, data);
    cout << result.sucess << endl;

    transport->close();
  } catch (TException& tx) {
    cout << "ERROR: " << tx.what() << endl;
  }   

   return;
}


void timestamp (int ll)

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}


