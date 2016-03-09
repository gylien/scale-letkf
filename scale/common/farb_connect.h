#ifndef _CONNECTLIB_H
#define _CONNECTLIB_H

#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#ifndef ENOERR
#define ENOERR 0
#endif
#include "farb_util.h"
#include "farb_conf_parse.h"

#include <semaphore.h>
#include <pthread.h>

//#define DEFAULT_NODE_BUF_SIZE (8*1024*1024)


struct server_connection{
  char service_name[1024];
  char port_name[MPI_MAX_PORT_NAME];
  MPI_Comm client;
  int server_size;		/* size of the server MPI_COMM_WORLD */
  int server_rank; 		/* rank in the server MPI_COMM_WORLD */
  int client_size;		/* size of the client MPI_COMM_WORLD */
  int mpi_tag;
  char* spare_buffer;
  size_t spare_buffer_size;
};

struct client_connection{
  char service_name[1024];
  char port_name[MPI_MAX_PORT_NAME];
  MPI_Comm server;
  int server_size;		/* size of the server MPI_COMM_WORLD */
  int client_size;		/* size of the client MPI_COMM_WORLD */
  int client_rank; 		/* rank in the client MPI_COMM_WORLD */
  int mpi_tag;
  char* spare_buffer;
  size_t spare_buffer_size;
};

enum req_status {
  initialized,
  being_sent,
  sent
};

#if 0
  char direction[8];    /* direction of send & recv */
  char alias_name[1024];    /* alias name for the file */
  int client_server_flag;
#endif
struct async_request {
  void* buffer;
  int count;
  MPI_Datatype datatype;
  int client_server_flag;
  int tag;
  MPI_Comm comm;

  enum req_status status;
  struct async_request*next;
};

#if 0
int save_to_buffer_node(struct file_buffer *fbuf, \
                int offset, int size, char * data);

int read_from_buffer_node(struct file_buffer *fbuf, \
                int offset, int size, char * data);

int flush_data_to_disk(char *file_name, struct file_buffer *fbuf);

int copy_buffer_file(struct file_buffer *src_fbuf, struct file_buffer *dst_fbuf);
#endif

/* initialize things */
void conn_init_();

/* finalize things */
void conn_finalize_();

/* wait for a connection */
struct server_connection* conn_serv_publish_port(const char * service_name, 
					char *port_name, MPI_Comm *client);

/* close a connection */
void conn_serv_unpublish_port(struct server_connection* connection);

/* connect to a server */
struct client_connection* conn_client_connect(const char* service_name, \
					MPI_Comm*server);

/* disconnect a client */
void conn_client_disconnect(struct client_connection*connection);


/* receive data from a server
 * (Currently for SYNC between server and client) 
 * should be replaced with netCDF related function 
 */
void conn_client_recv_data(struct client_connection*connection,
			   void* buffer,
			   int len,
			   int tag,
			   int src);

/* send data to a client 
 * (Currently for SYNC between server and client) 
 * should be replaced with netCDF related function 
 */
void conn_server_send_data(struct server_connection*connection,
			   void* buffer,
			   int len,
			   int tag,
			   int dest);

void pub_server_send_data2(char *buff, long size, int tag);

/* receive data from a client
 *  * (Currently for SYNC between server and client)
 *   * should be replaced with netCDF related function
 *    */
void conn_server_recv_data(struct server_connection*connection,
                           void* buffer,
                           int len,
                           int tag,
                           int src);

/* send data to a server
 *  * (Currently for SYNC between server and client)
 *   * should be replaced with netCDF related function
 *    */
void conn_client_send_data(struct client_connection*connection,
                           void* buffer,
                           int len,
                           int tag,
                           int dest);

void pub_client_send_data2(char *buff, long size, int tag);

/* --------------  functions used for SCALE-LETKF ----------------- */

void pub_server_connect_(char *service_name);

void pub_client_connect_(char *service_name);

int pub_send_data_(char *service_name);

int pub_recv_data_(char *service_name);

int pub_send_data1_(char *service_name);

int pub_recv_data1_(char *service_name);

void pub_server_send_sz(long size);
void pub_server_send_data(char *buff, long size);
long pub_server_recv_sz(void);
void pub_server_recv_data(char *buff, long);


void pub_netcdf_unpublish(void);

void pub_client_send_sz(long size);
void pub_client_send_data(char *buff, long size);
long pub_client_recv_sz(void);
void pub_client_recv_data(char *buff, long);
void pub_netcdf_disconnect(void);

struct file_buffer * pub_check_in_list(char *file_name);

struct async_request* async_mpi_send(const void *buf,
                     int count,
                     MPI_Datatype datatype,
                     int dest,
                     int tag,
                     MPI_Comm comm);

void* communication_thread(void*arg);


/* all to all diffusion between SIM and DA processes */
void conn_client_alltoall_data(struct client_connection*connection,
			       void* recv_buffer,
			       int send_count);

void conn_server_alltoall_data(struct server_connection*connection,
			       void* send_buffer,
			       int send_count);

/* all to all diffusion between SIM and DA processes that uses MPI_Scatter
 * to reduce the number of messages
 */
void conn_client_allscatter_data(struct client_connection*connection,
				  void* recv_buffer,
				  int send_count);

void conn_server_allscatter_data(struct server_connection*connection,
				  void* send_buffer,
				  int send_count);

/* Similar to allscatter_data, but uses MPI_Iscatter */
void conn_client_allscatter_nb_data(struct client_connection*connection,
				    void* recv_buffer,
				    int send_count);

void conn_server_allscatter_nb_data(struct server_connection*connection,
				    void* send_buffer,
				    int send_count);

#endif
