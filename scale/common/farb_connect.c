/*
 * Copyright (C) 2015, Advanced Institute for Computational Science, RIKEN
 * Author: Jianwei Liao(liaotoad@gmail.com)
 */

#include <string.h>
#include <assert.h>
#include "farb_connect.h"
#include <sys/time.h>
#include <dlfcn.h>

static struct server_connection *server_conn;
static struct client_connection *client_conn;

struct timeval tv[10];

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

void conn_init_(char *file_name, char *module_name) {
  //MPI_Init( NULL, NULL );
  parse_config("../../sample.ini", module_name);
}

void conn_finalize_() {
  delete_list();
}

struct server_connection *conn_serv_publish_port(const char *service_name,
						char *port_name,
						MPI_Comm *client){
  int ierr;

  struct server_connection *connection = malloc(sizeof(struct server_connection));
  strncpy(connection->service_name, service_name, 1024);

  MPI_Comm_size(MPI_COMM_WORLD, &connection->server_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &connection->server_rank);

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
  return connection;
}


void conn_serv_unpublish_port(struct server_connection *connection){
  MPI_Comm_disconnect( &connection->client );
  if(connection->server_rank == 0) {
    MPI_Unpublish_name(connection->service_name, MPI_INFO_NULL, \
		connection->port_name);
    MPI_Close_port(connection->port_name);
  }
}

struct client_connection* conn_client_connect(const char* service_name, 
					     MPI_Comm*server) {
  struct client_connection* connection = malloc(sizeof(struct client_connection));
  MPI_Comm_rank(MPI_COMM_WORLD, &connection->client_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &connection->client_size);

  strncpy(connection->service_name, service_name, 1024);
  int ierr = MPI_Lookup_name(connection->service_name, MPI_INFO_NULL, \
		connection->port_name);
  assert(ierr == MPI_SUCCESS);

  ierr = MPI_Comm_connect( connection->port_name, MPI_INFO_NULL, 0, \
		MPI_COMM_WORLD, &connection->server );
  assert(ierr == MPI_SUCCESS);

  ierr = MPI_Comm_remote_size(connection->server, &connection->server_size);
  assert(ierr == MPI_SUCCESS);
  return connection;
}

void conn_client_disconnect(struct client_connection *connection) {
  MPI_Comm_disconnect(&connection->server);
}

void conn_server_send_data(struct server_connection*connection,
			   void* buffer,
			   int len,
			   int tag,
			   int dest) {
  assert(dest < connection->client_size);
  int ierr = MPI_Send(buffer, len, MPI_BYTE, dest, tag, connection->client);
  assert(ierr == MPI_SUCCESS);
}

void conn_server_recv_data(struct server_connection*connection,
                           void* buffer,
                           int len,
                           int tag,
                           int src){
  assert(src < connection->server_size);
  MPI_Status status;
  MPI_Request request;
  int ierr = MPI_Recv(buffer, len, MPI_BYTE, src, MPI_ANY_TAG, \
                connection->client, &status);
  connection->mpi_tag = status.MPI_TAG;
  assert(ierr == MPI_SUCCESS);
}

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

/* ------------------- functions used for scale-letkf  -----------------------*/

void pub_server_connect_(char *service_name){
  char port_name[MPI_MAX_PORT_NAME];
  MPI_Comm client;
  server_conn = conn_serv_publish_port(service_name, port_name, &client);
}

void pub_client_connect_(char *service_name){
  MPI_Comm server;
  client_conn = conn_client_connect(service_name, &server);
}

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

int pub_send_data1_(char *service_name){

  void *handle;
  char *my_buf = NULL;
  int size;
  void * (*func_get_address)(void) = NULL;
  int (*func_get_size)(void) = NULL;
  
  handle = dlopen("../../hack.so", RTLD_LAZY); 
  if (!handle) {
	fprintf(stderr, "%s\n", dlerror());
	exit(EXIT_FAILURE);
  }

  func_get_address = (char (*)(void)) dlsym(handle, "get_address"); 
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

int pub_recv_data_(char *service_name){
  int i, j, file_len, loops;
  int list_size = list_get_size();
  struct file_buffer *fbuf;
  char *nbuf;
  for (i = 0; i < list_size; i++){
    fbuf = search_one_node(i);
	struct buffer_node *node;
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
		  //printf("client recv: tag %d\n", client_conn->mpi_tag);
		  add_to_buffer_node_list(fbuf, client_conn->mpi_tag/DEFAULT_NODE_BUF_SIZE, nbuf, 1);
	  	} 
	    printf("client recv: node %d, data size %d\n", get_buffer_list_size(fbuf), file_len);
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
          //printf("server recv: tag %d\n", client_conn->mpi_tag);
          add_to_buffer_node_list(fbuf, client_conn->mpi_tag/DEFAULT_NODE_BUF_SIZE, nbuf, 1);
        }
        printf("server recv: node %d, data size %d\n", get_buffer_list_size(fbuf), file_len);
        fbuf->buffer_pointer = NULL;
        fbuf->buffer_size =file_len;	
      }
  }
  return ENOERR;
}

int pub_send_data_(char *service_name){
  int i, j, node_size;
  int size = list_get_size();
  struct file_buffer *fbuf;
  gettimeofday(&tv[0], NULL);
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
  gettimeofday(&tv[1], NULL);
  printf("SEND TIME: %ld\n", (tv[1].tv_sec - tv[0].tv_sec) * 1000000 + \
		 tv[1].tv_usec - tv[0].tv_usec);
  return ENOERR; 
}

void pub_server_send_sz(long size){
  int k = server_conn->server_rank;
  conn_server_send_data(server_conn, &size, sizeof(long), 99, k);
}

void pub_server_send_data(char *buff, long size){
  int k = server_conn->server_rank;
  conn_server_send_data(server_conn, buff, size, 100, k);
}

void pub_server_send_data2(char *buff, long size, int tag){
  int k = server_conn->server_rank;
  conn_server_send_data(server_conn, buff, size, tag, k);
}

long pub_server_recv_sz(){
  char buff[16];
  pub_server_recv_data(buff, 16);
  return atol(buff);

}

void pub_server_recv_data(char *buff, long size){
    conn_server_recv_data(server_conn, buff, size, 100, server_conn->server_rank);
}

void pub_client_send_sz(long size){

  char buff[16];
  sprintf(buff, "%d", size);
  pub_client_send_data(buff, 16);

}

void pub_client_send_data2(char *buff, long size, int tag){
  int k = client_conn->client_rank;
  conn_client_send_data(client_conn, buff, size, tag, k);
}

void pub_client_send_data(char *buff, long size){
  int k = client_conn->client_rank;
  conn_client_send_data(client_conn, buff, size, 100, k);
}

long pub_client_recv_sz(void){
  long k = -1;
  conn_client_recv_data(client_conn, &k, sizeof(long), 99, client_conn->client_rank);
  return k;
}

void pub_client_recv_data(char *buff, long size){
  conn_client_recv_data(client_conn, buff, size, 100, client_conn->client_rank);
}

void pub_netcdf_unpublish(void){
  conn_serv_unpublish_port(server_conn);
}

void pub_netcdf_disconnect(void){
  conn_client_disconnect(client_conn);
}

struct file_buffer * pub_check_in_list(char *file_name){
  return  match_in_list(file_name);
}

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


