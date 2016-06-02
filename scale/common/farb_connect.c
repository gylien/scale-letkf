/*
 * Copyright (C) 2015, Advanced Institute for Computational Science, RIKEN
 * Author: Jianwei Liao(liaotoad@gmail.com)
 */

#include <string.h>
#include <assert.h>
#include "farb_connect.h"
#include <sys/time.h>
#include <dlfcn.h>

static struct server_connection *server_conn = NULL;
static struct client_connection *client_conn = NULL;

static void parse_config();

struct timeval tv[10];

#define DEFAULT_TAG 42

#if 0
int save_to_buffer_node(struct file_buffer *fbuf, \
                int offset, int size, char * data)
{
    int remain = size;
    int avail_buf;
    int tocpy;
    struct buffer_node *node;
    while (remain)
    {
        avail_buf = DEFAULT_NODE_BUF_SIZE - offset % DEFAULT_NODE_BUF_SIZE;
        tocpy = remain > avail_buf ? avail_buf : remain;
        node = search_buffer_node(fbuf, offset/DEFAULT_NODE_BUF_SIZE, NULL);
        if (node != NULL){
            memcpy(node->node_ptr + offset % DEFAULT_NODE_BUF_SIZE, data, tocpy);
            node->dirty_flag = 1;
        }
        else{
            char *nbuf = malloc(DEFAULT_NODE_BUF_SIZE);
            if (nbuf == NULL)
                return ENOMEM;
            memset(nbuf, 0, DEFAULT_NODE_BUF_SIZE);
            memcpy(nbuf + offset % DEFAULT_NODE_BUF_SIZE, data, tocpy);
            add_to_buffer_node_list(fbuf, offset/DEFAULT_NODE_BUF_SIZE, nbuf, 1);
        }
        data += tocpy;
        remain -= tocpy;
        offset += tocpy;
    }
    return 0;
}

int read_from_buffer_node(struct file_buffer *fbuf, \
                int offset, int size, char * data)
{

    int remain = size;
    int avail_buf;
    int tocpy;
    struct buffer_node *node;
    while (remain)
    {
        avail_buf = DEFAULT_NODE_BUF_SIZE - offset % DEFAULT_NODE_BUF_SIZE;
        tocpy = remain > avail_buf ? avail_buf : remain;
        node = search_buffer_node(fbuf, offset/DEFAULT_NODE_BUF_SIZE, NULL);
        if (node != NULL){
            memcpy(data, node->node_ptr + offset % DEFAULT_NODE_BUF_SIZE, tocpy);
        }
        else{
            data = NULL;
        }
        data += tocpy;
        remain -= tocpy;
        offset += tocpy;
    }
    return 0;
}

int flush_data_to_disk(char *file_name, struct file_buffer *fbuf){
    printf("flush data to the disk %s\n", file_name);
    FILE *f = fopen(file_name, "wb");
    struct buffer_node *node = (struct buffer_node*)fbuf->buffer_list;
    int ret, i = 0;

    fseek(f, 0, SEEK_END);
    fseek(f, 0, SEEK_SET);
    while (node != NULL){
        if (node->next == NULL && fbuf->buffer_size % DEFAULT_NODE_BUF_SIZE > 0){
            ret = fwrite(node->node_ptr, 1, fbuf->buffer_size % DEFAULT_NODE_BUF_SIZE, f);
            if(ret != fbuf->buffer_size % DEFAULT_NODE_BUF_SIZE){
                printf("\n fwrite() failed\n");
                fclose(f);
                return -1;
            }
        }
        else{
            ret = fwrite(node->node_ptr, 1, DEFAULT_NODE_BUF_SIZE, f);
            if(ret != DEFAULT_NODE_BUF_SIZE){
                printf("\n fwrite() failed\n");
                fclose(f);
                return -1;
            }

            if(0 != fseek(f, DEFAULT_NODE_BUF_SIZE, SEEK_CUR)){
                printf("\n fseek() failed\n");
                fclose(f);
                return -1;
            }
        }
        node = node->next;
        i++;
    }
    fclose(f);
    return 0;
}


int copy_buffer_file(struct file_buffer *src_fbuf, struct file_buffer *dst_fbuf){
    struct buffer_node *node = src_fbuf->buffer_list;
    char *nbuf;
    dst_fbuf->buffer_pointer = NULL;
    dst_fbuf->buffer_size = src_fbuf->buffer_size;
    int i = 0;
    while (node != NULL){
        if (node->next == NULL && src_fbuf->buffer_size % DEFAULT_NODE_BUF_SIZE > 0){
            nbuf = malloc(src_fbuf->buffer_size % DEFAULT_NODE_BUF_SIZE);
            if (nbuf==NULL)
                return ENOMEM;

            memcpy(nbuf, node->node_ptr, src_fbuf->buffer_size % DEFAULT_NODE_BUF_SIZE);
        }
        else{
            nbuf = malloc(DEFAULT_NODE_BUF_SIZE);
            if (nbuf==NULL)
                return ENOMEM;

            memcpy(nbuf, node->node_ptr, DEFAULT_NODE_BUF_SIZE);
        }
        add_to_buffer_node_list(dst_fbuf, i, nbuf, 1);
        node = node->next;
        i++;
    }
	return 0;
}
#endif


/*-------------------------------------------------------------------------*/
/**
  @brief	Parsing the config file for each component
  @param	ini_name		the name of the config file
  @param	service_name	the serivce name, remmended to use the name component
  @return 	void
 */
/*--------------------------------------------------------------------------*/

static void parse_config(char *ini_name, char *service_name){
  dictionary * ini ;
  char       * s; 
  char       para_name[1024];
  int        i, j;
  ini = iniparser_load(ini_name);
  iniparser_dump(ini, NULL);
  int comp_number =  iniparser_getint(ini, "COMPONENT:number", 0);
  int max_file_number = iniparser_getint(ini, "FILE:number", 0);

  for (j = 0; j < comp_number; j++){
  	sprintf(para_name, "COMP%d:service_name0", j);
  	s = iniparser_getstring(ini, para_name, NULL);
  	if (s == NULL || strcmp(s, service_name)){
	  continue;
  	}
  	for (i = 0; i < max_file_number; i++){
   	  sprintf(para_name, "COMP%d:file_name%d", j, i);
      s = iniparser_getstring(ini, para_name, NULL);
      if (s == NULL){
	  	break;
      }

   	add_to_list(s, 1);
   	struct file_buffer *fbuf = search_in_list(s, NULL);
   	sprintf(para_name, "COMP%d:service_name%d", j, i);
   	s = iniparser_getstring(ini, para_name, NULL);
   	strcpy(fbuf->service_name, s);

   	sprintf(para_name, "COMP%d:alias_name%d", j, i);
   	s = iniparser_getstring(ini, para_name, NULL);
   	strcpy(fbuf->alias_name, s);

   	sprintf(para_name, "COMP%d:direction%d", j, i);
   	s = iniparser_getstring(ini, para_name, NULL);
   	strcpy(fbuf->direction, s);

   	sprintf(para_name, "COMP%d:client_server_flag%d", j, i);
   	fbuf->client_server_flag = iniparser_getint(ini, para_name, 0);
	
	fbuf->buffer_list = NULL;
   	}
  }
   //print_all_list();
   iniparser_freedict(ini);
}

struct server_connection *conn_serv_publish_port(const char *service_name,
						char *port_name,
						MPI_Comm *client){
  int ierr;

  struct server_connection *connection = malloc(sizeof(struct server_connection));
  strncpy(connection->service_name, service_name, 1024);

  MPI_Comm_size(MPI_COMM_WORLD, &connection->server_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &connection->server_rank);

  connection->spare_buffer = NULL;
  connection->spare_buffer_size = 0;

  /* The root server published the sevice name */
  if(connection->server_rank==0) {
    strncpy(connection->port_name, port_name, MPI_MAX_PORT_NAME);

    ierr = MPI_Open_port(MPI_INFO_NULL, connection->port_name);
    assert(ierr== MPI_SUCCESS);

    ierr = MPI_Publish_name(connection->service_name, \
		MPI_INFO_NULL, connection->port_name);
    assert(ierr== MPI_SUCCESS);
  }

  MPI_Bcast(connection->port_name, MPI_MAX_PORT_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);

  ierr = MPI_Comm_accept(connection->port_name, MPI_INFO_NULL, 0, \
		MPI_COMM_WORLD, &connection->client);

  assert(ierr == MPI_SUCCESS);
  MPI_Comm_remote_size(connection->client, &connection->client_size);

  server_conn = connection;
  return connection;
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Release the opened port after application's execution
  @param	connection		the connection instance on the server side.
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_serv_unpublish_port(struct server_connection *connection){
  MPI_Comm_disconnect( &connection->client );
  if(connection->server_rank == 0) {
    MPI_Unpublish_name(connection->service_name, MPI_INFO_NULL, \
		connection->port_name);
    MPI_Close_port(connection->port_name);
  }
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Connect to the constructed inter-connection
  @param	service_name	name of the connection service
  @return 	server			record the server information after the connection
 */
/*--------------------------------------------------------------------------*/
struct client_connection* conn_client_connect(const char* service_name, 
					     MPI_Comm*server) {
  struct client_connection* connection = malloc(sizeof(struct client_connection));
  MPI_Comm_rank(MPI_COMM_WORLD, &connection->client_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &connection->client_size);

  connection->spare_buffer = NULL;
  connection->spare_buffer_size = 0;

  strncpy(connection->service_name, service_name, 1024);
  int ierr = MPI_Lookup_name(connection->service_name, MPI_INFO_NULL, \
		connection->port_name);
  assert(ierr == MPI_SUCCESS);

  ierr = MPI_Comm_connect( connection->port_name, MPI_INFO_NULL, 0, \
		MPI_COMM_WORLD, &connection->server );
  assert(ierr == MPI_SUCCESS);

  ierr = MPI_Comm_remote_size(connection->server, &connection->server_size);
  assert(ierr == MPI_SUCCESS);
  client_conn = connection;
  return connection;
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Disconect inter-connection from client side
  @param	connection		the connection instance on the client side.
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_client_disconnect(struct client_connection *connection) {
  MPI_Comm_disconnect(&connection->server);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Sending data in the inter-connection environment on the server side
  @param	connection		the connection instance on the server side
  @param	buffer			data for sending
  @param	len				length for sending
  @param	tag				tag of the data
  @param	dest			the rank of the destination process
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_server_send_data(struct server_connection*connection,
			   void* buffer,
			   int len,
			   int tag,
			   int dest) {
  assert(dest < connection->client_size);
  int ierr = MPI_Send(buffer, len, MPI_BYTE, dest, tag, connection->client);
  assert(ierr == MPI_SUCCESS);
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Receving data in the inter-connection environment on the client side
  @param	connection		the connection instance on the client side
  @param	buffer			buffer for storing the received data
  @param	len				length for receiving data
  @param	src				the rank of the source process
  @return	void
 */
/*--------------------------------------------------------------------------*/
void conn_client_recv_data(struct client_connection*connection,
                           void* buffer,
                           int len,
                           int tag,
                           int src) {

  assert(src < connection->server_size);
  MPI_Status status;
  int ierr = MPI_Recv(buffer, len, MPI_BYTE, src, MPI_ANY_TAG, \
                connection->server, &status);
  connection->mpi_tag = status.MPI_TAG;
  assert(ierr == MPI_SUCCESS);
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Receving data in the inter-connection environment on the server side
  @param	connection		the connection instance on the server side
  @param	buffer			buffer for storing the received data
  @param	len				length for receiving data
  @param	src				the rank of the source process
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_server_recv_data(struct server_connection*connection,
                           void* buffer,
                           int len,
                           int tag,
                           int src){
  assert(src < connection->server_size);
  MPI_Status status;
  int ierr = MPI_Recv(buffer, len, MPI_BYTE, src, MPI_ANY_TAG, \
                connection->client, &status);
  connection->mpi_tag = status.MPI_TAG;
  assert(ierr == MPI_SUCCESS);
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Sending data in the inter-connection environment on the client side
  @param	connection		the connection instance on the client side
  @param	buffer			data for sending
  @param	len				length for sending
  @param	tag				tag of the data
  @param	dest			the rank of the destination process
  @return	void
 */
/*--------------------------------------------------------------------------*/
void conn_client_send_data(struct client_connection*connection,
                           void* buffer,
                           int len,
                           int tag,
                           int dest){
  assert(dest < connection->client_size);
  int ierr = MPI_Send(buffer, len, MPI_BYTE, dest, tag, \
		connection->server);
  
  assert(ierr == MPI_SUCCESS);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Alltoall Comm. in the inter-connection environment on the server side
  @param	connection		the connection instance on the server side
  @param	send_buffer		the buffer for storing the data requiring to be sent
  @param	send_count		length for each process in alltoall
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_server_alltoall_data(struct server_connection*connection,
			       void* send_buffer,
			       int send_count) {
  size_t buffer_size = send_count*connection->client_size;
    if(connection->spare_buffer_size < buffer_size) {
    connection->spare_buffer = realloc(connection->spare_buffer, buffer_size);
    assert(connection->spare_buffer);
    connection->spare_buffer_size = send_count*buffer_size;
  }

  int ierr = MPI_Alltoall(send_buffer, send_count, MPI_BYTE,
			  connection->spare_buffer, send_count, MPI_BYTE,
			  connection->client);
  assert(ierr == MPI_SUCCESS);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Alltoall Comm. in the inter-connection environment on the client side
  @param	connection		the connection instance on the client side
  @param	recv_buffer		the buffer for storing the received data in alltoall
  @param	send_count		length for each process in alltoall
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_client_alltoall_data(struct client_connection*connection,
			       void* recv_buffer,
			       int send_count) {
  size_t buffer_size = send_count * connection->server_size;
  if(connection->spare_buffer_size < buffer_size) {
    connection->spare_buffer = realloc(connection->spare_buffer, buffer_size);
    assert(connection->spare_buffer);
    connection->spare_buffer_size = buffer_size;
  }

  int ierr = MPI_Alltoall(connection->spare_buffer, send_count, MPI_BYTE,
			  recv_buffer, send_count, MPI_BYTE,
			  connection->server);
  assert(ierr == MPI_SUCCESS);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Scattering (send) the data on the server side
  @param	connection		the connection instance on the server side
  @param	send_buffer		the buffer for the data requiring to be sent
  @param	send_count		length for each process in alltoall
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_server_allscatter_data(struct server_connection*connection,
				 void* send_buffer,
				 int send_count) {
  int i;
  for(i=0; i<connection->server_size; i++) {
    int root = i == connection->server_rank ? MPI_ROOT : MPI_PROC_NULL;
    int ierr = MPI_Scatter(send_buffer, send_count, MPI_BYTE,
			   NULL, 0, MPI_BYTE,
			   root,
			   connection->client);
    assert(ierr == MPI_SUCCESS);
  }
}



/*-------------------------------------------------------------------------*/
/**
  @brief	Scattering (recv) the data on the client side (MPI_Scatter Impl.)
  @param	connection		the connection instance on the client side
  @param	recv_buffer		the buffer for storing the received data
  @param	send_count		length for each process in alltoall
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_client_allscatter_data(struct client_connection*connection,
				  void* recv_buffer,
				  int send_count) {
  int i;
  for(i=0; i<connection->server_size; i++) {
    int ierr = MPI_Scatter(NULL, 0, MPI_BYTE,
			   &((char*)recv_buffer)[i*send_count], send_count, MPI_BYTE,
			   i,
			   connection->server);
    assert(ierr == MPI_SUCCESS);
  }
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Scattering (send) the data on the server side (MPI_Iscatter Impl.)
  @param	connection		the connection instance on the server side
  @param	send_buffer		the buffer for the data requiring to be sent
  @param	send_count		length for each process in alltoall
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_server_allscatter_nb_data(struct server_connection*connection,
				  void* send_buffer,
				  int send_count) {
#if HAVE_MPI3
  MPI_Request reqs[connection->server_size];
  int i;
  for(i=0; i<connection->server_size; i++) {
    int root = i == connection->server_rank ? MPI_ROOT : MPI_PROC_NULL;
    int ierr = MPI_Iscatter(send_buffer, send_count, MPI_BYTE,
			   NULL, 0, MPI_BYTE,
			   root,
			    connection->client,
			    &reqs[i]);
    assert(ierr == MPI_SUCCESS);
  }
  MPI_Waitall(i, reqs, MPI_STATUSES_IGNORE);
#else
  conn_server_allscatter_data(connection, send_buffer, send_count);
#endif
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Scattering (recv) the data on the client side (MPI_Iscatter Impl.)
  @param	connection		the connection instance on the client side
  @param	recv_buffer		the buffer for storing the received data
  @param	send_count		length for each process in alltoall
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_client_allscatter_nb_data(struct client_connection*connection,
				  void* recv_buffer,
				  int send_count) {
#if HAVE_MPI3
  MPI_Request reqs[connection->server_size];
  int i;
  for(i=0; i<connection->server_size; i++) {
    int ierr = MPI_Iscatter(NULL, 0, MPI_BYTE,
			   &((char*)recv_buffer)[i*send_count], send_count, MPI_BYTE,
			   i,
			    connection->server,
			    &reqs[i]);
    assert(ierr == MPI_SUCCESS);
  }
  MPI_Waitall(i, reqs, MPI_STATUSES_IGNORE);
#else
  conn_client_allscatter_data(connection, recv_buffer, send_count);
#endif
}

/******************************************************************************/
/* 
 *  functions used for scale-letkf											  *  
*/ 
/******************************************************************************/


/*-------------------------------------------------------------------------*/
/**
  @brief	alltoallv communication (recv) targeting the 3D array on the client side
  @param	connection		the connection instance on the client side
  @param	recv_buffer		the buffer for storing the received data
  @param	datatype		datatype of each send buffer element (handle)
  @param	buffer_size		Each process may send a different amount of data (alltoallv)
  @param	recv_count		length for each process in alltoallv
  @param	axis			(Important) 0 for normal case, and 2 for SCALE-LETKF 
  @return	void
 */
/*--------------------------------------------------------------------------*/

void pub_client_allscatter_3d_range_data(struct client_connection*connexion,
				       int nb_groups,
				       int* recv_buffer,
				       MPI_Datatype datatype,
				       size_t buffer_size[3],
				       int recv_count,
				       int axis) {
  int my_range = connexion->client_rank % nb_groups;
  int i;
  MPI_Datatype row;
  int ierr;
  switch(axis) {
  case 0: // blue
    ierr = MPI_Type_vector(recv_count, // nb blocks
			   buffer_size[0]*buffer_size[1], //block_length
			   0, // stride
			   datatype,
			   &row);
    break;
  case 1: // green
    ierr = MPI_Type_vector(recv_count*buffer_size[2], // nb blocks
			   buffer_size[0], //block_length
			   buffer_size[0]*buffer_size[1], // stride
			   //buffer_size[1], // stride
			   datatype,
			   &row);
    break;
  case 2:
    ierr = MPI_Type_vector(recv_count*buffer_size[1]*buffer_size[2], // nb blocks
			   1, //block_length
			   buffer_size[0], // stride
			   datatype,
			   &row);
    break;
  default:
    fprintf(stderr, "%s: Incorrect axis (%d)\n", __FUNCTION__, axis);
    abort();
  }
  assert(ierr == MPI_SUCCESS);
  ierr = MPI_Type_commit(&row);
  assert(ierr == MPI_SUCCESS);

  int nreqs = connexion->server_size/nb_groups;
  MPI_Request reqs[nreqs];
  for(i=0; i<connexion->server_size/nb_groups; i++) {
    int src = (nb_groups*i)+my_range;
    int tag = DEFAULT_TAG;
    void* start_addr = NULL;
    switch(axis) {
    case 0: // blue
      start_addr =&recv_buffer[buffer_size[0]*buffer_size[1]*i*recv_count];
      break;
    case 1:
      start_addr =&recv_buffer[buffer_size[0]*i*recv_count];
      break;
    case 2:
      start_addr =&recv_buffer[i*recv_count];
      break;
    default:
      abort();
    }

    int ierr = MPI_Irecv(start_addr, 1, row, src, tag, connexion->server, &reqs[i]);
    assert(ierr == MPI_SUCCESS);
  }
  MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	alltoallv communication (send) targeting the 3D array on the client side
  @param	connection		the connection instance on the server side
  @param	send_buffer		the buffer for storing the sent data
  @param	datatype		datatype of each send buffer element (handle)
  @param	buffer_size		Each process may send a different amount of data (alltoallv)
  @param	send_count		length for each process in alltoallv
  @param	axis			(Important) 0 for normal case, and 2 for SCALE-LETKF 
  @return	void
 */
/*--------------------------------------------------------------------------*/

void pub_server_allscatter_3d_range_data(struct server_connection*connexion,
					 int nb_groups,
					 int* send_buffer,
					 MPI_Datatype datatype,
					 size_t buffer_size[3],
					 int send_count,
					 int axis) {
  int my_range = connexion->server_rank % nb_groups;
  int i;
  MPI_Datatype row;
  int ierr;

  switch(axis) {
  case 0: // blue
    ierr = MPI_Type_vector(send_count, // nb blocks
			   buffer_size[0]*buffer_size[1], //block_length
			   0, // stride
			   datatype,
			   &row);
    break;
  case 1: // green
    ierr = MPI_Type_vector(send_count*buffer_size[2], // nb blocks
			   buffer_size[0], //block_length
			   buffer_size[0]*buffer_size[1], // stride
			   datatype,
			   &row);
    break;
  case 2:
    ierr = MPI_Type_vector(send_count*buffer_size[1]*buffer_size[2], // nb blocks
			   1, //block_length
			   buffer_size[0], // stride
			   datatype,
			   &row);
    break;
  default:
    fprintf(stderr, "%s: Incorrect axis (%d)\n", __FUNCTION__, axis);
    abort();
  }
  assert(ierr == MPI_SUCCESS);
  ierr = MPI_Type_commit(&row);
  assert(ierr == MPI_SUCCESS);

  int nreqs = connexion->client_size/nb_groups;
  MPI_Request reqs[nreqs];

  for(i=0; i<connexion->client_size/nb_groups; i++) {
    int dest = (nb_groups*i)+my_range;
    int tag = DEFAULT_TAG;
    void* start_addr = NULL;
    switch(axis) {
    case 0: // blue
      start_addr =&send_buffer[buffer_size[0]*buffer_size[1]*i*send_count];
      break;
    case 1:
      start_addr =&send_buffer[buffer_size[0]*i*send_count];
      break;
    case 2:
      start_addr=&send_buffer[i*send_count];
      break;
    default:
      abort();
    }
    int ierr = MPI_Isend(start_addr, 1, row, dest, tag, connexion->client, &reqs[i]);
    assert(ierr == MPI_SUCCESS);
  }
  MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	alltoallv communication (recv) targeting the 2D array on the client side
  @param	connection		the connection instance on the client side
  @param	recv_buffer		the buffer for storing the received data
  @param	datatype		datatype of each send buffer element (handle)
  @param	buffer_size		Each process may send a different amount of data (alltoallv)
  @param	recv_count		length for each process in alltoallv
  @param	axis			(Important) 0 for normal case, and 1 for SCALE-LETKF 
  @return	void
 */
/*--------------------------------------------------------------------------*/

void pub_client_allscatter_2d_range_data(struct client_connection*connexion,
				       int nb_groups,
				       int* recv_buffer,
				       MPI_Datatype datatype,
				       size_t buffer_size[2],
				       int recv_count,
				       int axis) {
  switch(axis) {
  case 0:
    pub_client_allscatter_2d_row_range_data(connexion, nb_groups, recv_buffer,
					   datatype, buffer_size, recv_count);
    break;
  case 1:
    pub_client_allscatter_2d_col_range_data(connexion, nb_groups, recv_buffer,
					  datatype, buffer_size, recv_count);
    break;
  default:
    fprintf(stderr, "%s: Incorrect axis (%d)\n", __FUNCTION__, axis);
    abort();
  }
}

/*-------------------------------------------------------------------------*/
/**
  @brief	alltoallv communication (send) targeting the 2D array on the client side
  @param	connection		the connection instance on the server side
  @param	send_buffer		the buffer for storing the sent data
  @param	datatype		datatype of each send buffer element (handle)
  @param	buffer_size		Each process may send a different amount of data (alltoallv)
  @param	send_count		length for each process in alltoallv
  @param	axis			(Important) 0 for normal case, and 1 for SCALE-LETKF 
  @return	void
 */
/*--------------------------------------------------------------------------*/


void pub_server_allscatter_2d_range_data(struct server_connection*connexion,
					 int nb_groups,
					 int* send_buffer,
					 MPI_Datatype datatype,
					 size_t buffer_size[2],
					 int send_count,
					 int axis) {
  switch(axis) {
  case 0:
    pub_server_allscatter_2d_row_range_data(connexion, nb_groups, send_buffer,
					    datatype, buffer_size, send_count);
    break;
  case 1:
    pub_server_allscatter_2d_col_range_data(connexion, nb_groups, send_buffer,
					    datatype, buffer_size, send_count);
    break;
  default:
    fprintf(stderr, "%s: Incorrect axis (%d)\n", __FUNCTION__, axis);
    abort();
  }
}


/* scatter the columns of a 2D array  */
void pub_server_allscatter_2d_row_range_data(struct server_connection*connexion,
					     int nb_groups,
					     int* send_buffer,
					     MPI_Datatype datatype,
					     size_t buffer_size[2],
					     int row_send_count) {

  int my_range = connexion->server_rank % nb_groups;
  int i;

  MPI_Datatype row;
  int ierr = MPI_Type_vector(row_send_count,
			     buffer_size[0],
			     0,
			     datatype,
			     &row);
  assert(ierr == MPI_SUCCESS);
  ierr = MPI_Type_commit(&row);
  assert(ierr == MPI_SUCCESS);

  int nreqs = connexion->client_size/nb_groups;
  MPI_Request reqs[nreqs];

  for(i=0; i<connexion->client_size/nb_groups; i++) {
    int dest = (nb_groups*i)+my_range;
    int tag = DEFAULT_TAG;
    int ierr = MPI_Isend(&send_buffer[buffer_size[0]*i*row_send_count], 1, row, dest, \
									tag, connexion->client, &reqs[i]);
    assert(ierr == MPI_SUCCESS);
  }
  MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
}


void pub_client_allscatter_2d_row_range_data(struct client_connection*connexion,
					   int nb_groups,
					   int* recv_buffer,
					   MPI_Datatype datatype,
					   size_t buffer_size[2],
					   int row_recv_count) {
  int my_range = connexion->client_rank % nb_groups;
  int i;
  MPI_Datatype row;
  int ierr = MPI_Type_vector(row_recv_count,
			     buffer_size[0],
			     0,
			     datatype,
			     &row);
  assert(ierr == MPI_SUCCESS);
  ierr = MPI_Type_commit(&row);
  assert(ierr == MPI_SUCCESS);

  int nreqs = connexion->server_size/nb_groups;
  MPI_Request reqs[nreqs];
  for(i=0; i<connexion->server_size/nb_groups; i++) {
    int src = (nb_groups*i)+my_range;
    int tag = DEFAULT_TAG;
    int ierr = MPI_Irecv(&recv_buffer[buffer_size[0]*row_recv_count*i], 1, row, src, \
									tag, connexion->server, &reqs[i]);
    assert(ierr == MPI_SUCCESS);
  }
  MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
}

/* scatter the columns of a 2D array  */
void pub_server_allscatter_2d_col_range_data(struct server_connection*connexion,
					     int nb_groups,
					     int* send_buffer,
					     MPI_Datatype datatype,
					     size_t buffer_size[2],
					     int col_send_count) {

  int my_range = connexion->server_rank % nb_groups;
  int i;

  MPI_Datatype col;
  int ierr = MPI_Type_vector(buffer_size[0],
			     col_send_count,
			     buffer_size[1],
			     datatype,
			     &col);
  assert(ierr == MPI_SUCCESS);
  ierr = MPI_Type_commit(&col);
  assert(ierr == MPI_SUCCESS);

  int nreqs = connexion->client_size/nb_groups;
  MPI_Request reqs[nreqs];

  for(i=0; i<connexion->client_size/nb_groups; i++) {
    int dest = (nb_groups*i)+my_range;
    int tag = DEFAULT_TAG;
    int ierr = MPI_Isend(&send_buffer[i*col_send_count], 1, col, dest, \
							tag, connexion->client, &reqs[i]);
    assert(ierr == MPI_SUCCESS);
  }
  MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
}


void pub_client_allscatter_2d_col_range_data(struct client_connection*connexion,
					   int nb_groups,
					   int* recv_buffer,
					   MPI_Datatype datatype,
					   size_t buffer_size[2],
					   int col_recv_count) {
  int my_range = connexion->client_rank % nb_groups;
  int i;
  MPI_Datatype col;
  int ierr = MPI_Type_vector(buffer_size[0],
			     col_recv_count,
			     buffer_size[1],
			     datatype,
			     &col);
  assert(ierr == MPI_SUCCESS);
  ierr = MPI_Type_commit(&col);
  assert(ierr == MPI_SUCCESS);

  int nreqs = connexion->server_size/nb_groups;
  MPI_Request reqs[nreqs];
  for(i=0; i<connexion->server_size/nb_groups; i++) {
    int src = (nb_groups*i)+my_range;
    int tag = DEFAULT_TAG;
    int ierr = MPI_Irecv(&recv_buffer[col_recv_count*i], 1, col, src, \
						tag, connexion->server, &reqs[i]);
    assert(ierr == MPI_SUCCESS);
  }
  MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
}

void pub_server_allscatter_range_data(struct server_connection*connexion,
				      int nb_groups,
				      char* send_buffer,
				      int send_count) {
  int my_range = connexion->server_rank % nb_groups;
  int i;
  MPI_Request reqs[connexion->client_size/nb_groups];

  for(i=0; i<connexion->client_size/nb_groups; i++) {
    int dest = (nb_groups*i)+my_range;
    int tag = DEFAULT_TAG;
    int ierr = MPI_Isend(&send_buffer[send_count*i], send_count, MPI_BYTE, dest,\
						tag, connexion->client, &reqs[i]);
    assert(ierr == MPI_SUCCESS);
  }
  MPI_Waitall(connexion->client_size/nb_groups, reqs, MPI_STATUSES_IGNORE);
}

void pub_client_allscatter_range_data(struct client_connection*connexion,
				  int nb_groups,
				  char* recv_buffer,
				  int recv_count) {

  int my_range = connexion->client_rank % nb_groups;
  int i;
  MPI_Request reqs[connexion->server_size/nb_groups];
  for(i=0; i<connexion->server_size/nb_groups; i++) {
    int src = (nb_groups*i)+my_range;
    int tag = DEFAULT_TAG;
    int ierr = MPI_Irecv(&recv_buffer[recv_count*i], recv_count, MPI_BYTE, src, \
						tag, connexion->server, &reqs[i]);
    assert(ierr == MPI_SUCCESS);
  }
  MPI_Waitall(connexion->server_size/nb_groups, reqs, MPI_STATUSES_IGNORE);
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Sending an int about the length of data on the server side
  @param 	size		the integer of size
  @return 	void
 */
/*--------------------------------------------------------------------------*/

void pub_server_send_sz(size_t size){
  assert(server_conn);
  int k = server_conn->server_rank;
  conn_server_send_data(server_conn, &size, sizeof(size), 99, k);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Sending the data on the server side (used for Sync Transfer)
  @param	buff		the buffer for storing the sent data
  @param	size		maximum number of elements in the buffer
  @return 	void
 */
/*--------------------------------------------------------------------------*/

void pub_server_send_data(void*buff, size_t size){
  assert(server_conn);
  int k = server_conn->server_rank;
  conn_server_send_data(server_conn, buff, size, 100, k);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Sending the data on the server side (used for Async Transfer)
  @param	buff		the buffer for storing the sent data
  @param	size		maximum number of elements in the buffer
  @param	tag			the tag of this message
  @return 	void
 */
/*--------------------------------------------------------------------------*/

void pub_server_send_data2(void *buff, size_t size, int tag){
  assert(server_conn);
  int k = server_conn->server_rank;
  conn_server_send_data(server_conn, buff, size, tag, k);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Receiving an int about the length of data on the server side
  @return 	size
 */
/*--------------------------------------------------------------------------*/

size_t pub_server_recv_sz(){
  assert(server_conn);
  size_t s;
  pub_server_recv_data(&s, sizeof(s));
  return s;
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Receiving the data on the server side
  @param	buff		the buffer for storing the received data
  @param	size		maximum number of elements in the receive buffer
  @return 	void
 */
/*--------------------------------------------------------------------------*/

void pub_server_recv_data(void *buff, size_t size){
  assert(server_conn);
  conn_server_recv_data(server_conn, buff, size, 100, server_conn->server_rank);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Sending an int about the length of data on the client side
  @param 	size		the integer of size
  @return 	void
 */
/*--------------------------------------------------------------------------*/

void pub_client_send_sz(size_t size){
  assert(client_conn);
  pub_client_send_data(&size, sizeof(size));
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Sending the data on the client side (used for Async Transfer)
  @param	buff		the buffer for storing the sent data
  @param	size		maximum number of elements in the buffer
  @param	tag			the tag of this message
  @return 	void
 */
/*--------------------------------------------------------------------------*/

void pub_client_send_data2(void *buff, size_t size, int tag){
  assert(client_conn);
  int k = client_conn->client_rank;
  conn_client_send_data(client_conn, buff, size, tag, k);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Sending the data on the client side (used for Sync Transfer)
  @param	buff		the buffer for storing the sent data
  @param	size		maximum number of elements in the buffer
  @return 	void
 */
/*--------------------------------------------------------------------------*/

void pub_client_send_data(void *buff, size_t size){
  assert(client_conn);
  int k = client_conn->client_rank;
  conn_client_send_data(client_conn, buff, size, 100, k);
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Receiving an int about the length of data on the client side
  @return 	size
 */
/*--------------------------------------------------------------------------*/

size_t pub_client_recv_sz(void){
  assert(client_conn);
  size_t k = -1;
  conn_client_recv_data(client_conn, &k, sizeof(k), 99, client_conn->client_rank);
  return k;
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Receiving the data on the client side
  @param	buff		the buffer for storing the received data
  @param	size		maximum number of elements in the receive buffer
  @return 	void
 */
/*--------------------------------------------------------------------------*/

void pub_client_recv_data(void *buff, size_t size){
  assert(client_conn);
  conn_client_recv_data(client_conn, buff, size, 100, client_conn->client_rank);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Wrapper of "conn_serv_unpublish_port"
 */
/*--------------------------------------------------------------------------*/

void pub_netcdf_unpublish(void){
  conn_serv_unpublish_port(server_conn);
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Wrapper of "conn_client_disconnect)"
 */
/*--------------------------------------------------------------------------*/

void pub_netcdf_disconnect(void){
  conn_client_disconnect(client_conn);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Wrapper of "match_in_list(char * file_name)"
  @param	file_name		the name of the target netcdf file
 */
/*--------------------------------------------------------------------------*/

struct file_buffer * pub_check_in_list(char *file_name){
  return  match_in_list(file_name);
}


/******************************************************************************/
/* 
 *  functions used for scale-letkf, and called by Fortran					  *  
*/ 
/******************************************************************************/

/*-------------------------------------------------------------------------*/
/**
  @brief	Initialize component (module) by parsing the config file.
  @param	file_name		the name of configure file.
  @param	module_name		the name of module, should be same in the config file.
  @return	void
 */
/*--------------------------------------------------------------------------*/
void conn_init_(char *file_name, char *module_name) {
  //MPI_Init( NULL, NULL );
  parse_config(file_name, module_name);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Finalize the component, and delete the list
  @return	void
 */
/*--------------------------------------------------------------------------*/

void conn_finalize_() {
  delete_list();
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Publish the port by rank 0 server process
  @param	service_name	name of the connection service
  @return	port_name 		name of the opened port
  @return 	client			record the client information after the connection
 */
/*--------------------------------------------------------------------------*/


void pub_server_connect_(char *service_name){
  char port_name[MPI_MAX_PORT_NAME];
  MPI_Comm client;
  server_conn = conn_serv_publish_port(service_name, port_name, &client);
}

void pub_client_connect_(char *service_name){
  MPI_Comm server;
  client_conn = conn_client_connect(service_name, &server);
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Receiving the data relevant to the component, and is called by Fortran
  @param	server_name		the name of the component
  @return 	void
 */
/*--------------------------------------------------------------------------*/


int pub_recv_data_(char *service_name){
  int i, j, file_len, loops;
  int list_size = list_get_size();
  struct file_buffer *fbuf;
  char *nbuf;
  for (i = 0; i < list_size; i++){
    fbuf = search_one_node(i);
    if (fbuf && !strncmp(fbuf->direction, "rd", 2) && (fbuf->client_server_flag == 0)){
	  	file_len = pub_client_recv_sz();
		if (file_len % DEFAULT_NODE_BUF_SIZE > 0)
			loops = file_len/DEFAULT_NODE_BUF_SIZE + 1;
		else loops = file_len/DEFAULT_NODE_BUF_SIZE;

		for (j = 0; j < loops; j++)
		{
          nbuf = malloc(DEFAULT_NODE_BUF_SIZE);
          if (nbuf == NULL)
            return ENOMEM;
	
      	  pub_client_recv_data(nbuf, DEFAULT_NODE_BUF_SIZE);
		  add_to_buffer_node_list(fbuf, client_conn->mpi_tag/DEFAULT_NODE_BUF_SIZE, nbuf, 1);
	  	} 
	    fbuf->buffer_pointer = NULL;
        fbuf->buffer_size =file_len;
      }
      else if (!strncmp(fbuf->direction, "rd", 2)){
        file_len = pub_server_recv_sz();
        if (file_len % DEFAULT_NODE_BUF_SIZE > 0)
            loops = file_len/DEFAULT_NODE_BUF_SIZE + 1;
        else loops = file_len/DEFAULT_NODE_BUF_SIZE;

        for (j = 0; j < loops; j++)
        {
          nbuf = malloc(DEFAULT_NODE_BUF_SIZE);
          if (nbuf == NULL)
            return ENOMEM;

          pub_server_recv_data(nbuf, DEFAULT_NODE_BUF_SIZE);
          add_to_buffer_node_list(fbuf, client_conn->mpi_tag/DEFAULT_NODE_BUF_SIZE, nbuf, 1);
        }
        fbuf->buffer_pointer = NULL;
        fbuf->buffer_size =file_len;	
      }
  }
  return ENOERR;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Sending the data relevant to the component, and is called by Fortran
  @param	server_name		the name of the component
  @return 	void
 */
/*--------------------------------------------------------------------------*/

int pub_send_data_(char *service_name){
  int i, j, node_size;
  int size = list_get_size();
  struct file_buffer *fbuf;
  //gettimeofday(&tv[0], NULL);
  for (i = 0; i < size; i++){
    fbuf = search_one_node(i);
    struct buffer_node *node = fbuf->buffer_list;

    if (fbuf && !strcmp(fbuf->direction, "wr")){
	  node_size = get_buffer_list_size(fbuf);

      if (fbuf->client_server_flag == 0){
		pub_client_send_sz(fbuf->buffer_size);
		for (j = 0; j < node_size; j++){ 
		  if (node && node->dirty_flag){
		    pub_client_send_data2(node->node_ptr, DEFAULT_NODE_BUF_SIZE, j*DEFAULT_NODE_BUF_SIZE);
	        //printf("client send, node %d is dirty\n", j);
		  }
		  node = node->next;
		}
      }
      else{
		pub_server_send_sz(fbuf->buffer_size);
        for (j = 0; j < node_size; j++){
          if (node && node->dirty_flag){
            pub_server_send_data2(node->node_ptr, DEFAULT_NODE_BUF_SIZE, j*DEFAULT_NODE_BUF_SIZE);
			//printf("server_send, node %d is dirty\n", j);
		  }
          node = node->next;
		}
	  }
    }
  }
  //gettimeofday(&tv[1], NULL);
  //printf("SEND TIME: %ld\n", (tv[1].tv_sec - tv[0].tv_sec) * 1000000 + \
		 tv[1].tv_usec - tv[0].tv_usec);
  return ENOERR; 
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Receiving the data relevant to the component, and is called by Fortran
  			(used for receiving the obsda data: OBSOPE -> LETKF)
  @param	server_name		the name of the component
  @return 	void
 */
/*--------------------------------------------------------------------------*/

int pub_recv_data1_(char *service_name){

  void *handle;
  int size;
  void  (*func_set_address)(char *) = NULL;
  void (*func_set_size)(int) = NULL;
  char *local_buffer = NULL;

  size = pub_client_recv_sz();
  
  if (size != 0){
  	local_buffer =  malloc(size);
  	if (local_buffer == NULL)
  		return ENOMEM;
  	pub_client_recv_data(local_buffer, size);
  }
   
  handle = dlopen("../../hack.so", RTLD_LAZY);
  if (!handle) {
	fprintf(stderr, "%s\n", dlerror());
	exit(EXIT_FAILURE);
  }
  func_set_address = (void(*)(char *)) dlsym(handle, "set_address"); 
  func_set_address(local_buffer);
  func_set_size = (void(*)(int)) dlsym(handle, "set_size");
  func_set_size(size);
  dlclose(handle);
  return size;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Sending the data relevant to the component, and is called by Fortran
  			(used for sending the obsda data: OBSOPE -> LETKF)
  @param	server_name		the name of the component
  @return 	void
 */
/*--------------------------------------------------------------------------*/

int pub_send_data1_(char *service_name){

  void *handle;
  char *my_buf = NULL;
  int size;
  char* (*func_get_address)(void) = NULL;
  int (*func_get_size)(void) = NULL;
  
  handle = dlopen("../../hack.so", RTLD_LAZY); 
  if (!handle) {
	fprintf(stderr, "%s\n", dlerror());
	exit(EXIT_FAILURE);
  }

  func_get_address = (char* (*)(void)) dlsym(handle, "get_address"); 
  my_buf = func_get_address();

  func_get_size = (int(*)(void)) dlsym(handle, "get_size");
  size = func_get_size();
  
  pub_server_send_sz(size);
  if (size != 0)
  	pub_server_send_data(my_buf, size);
  dlclose(handle);
  
  /* free the buffer of the send processes */
  
  if (my_buf != NULL && size !=0)
  	free(my_buf);
  return size;
}

