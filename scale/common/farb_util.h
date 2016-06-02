/*
 * Copyright (C) 2015, Advanced Institute for Computational Science, RIKEN
 * Author: Jianwei Liao(liaotoad@gmail.com)
 */

#ifndef LINKED_LIST_H
#define LINKED_LIST_H

#define DEFAULT_NODE_BUF_SIZE (4*1024*1024)

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

typedef struct buffer_node{
	int offset_start;
	char *node_ptr;
  	int dirty_flag;
	struct buffer_node *next;
}node_buffer_t;

struct file_buffer{
  char file_path[1024];         /* path of the file */
  char *buffer_pointer;         /* pointer of the buffer */
  struct buffer_node *buffer_list; /* buffer node list */
  int buffer_size;              /* size of the buffer */
  char service_name[1024];		/* service name */
  char direction[8];			/* direction of send & recv */
  char alias_name[1024];		/* alias name for the file */
  int client_server_flag;
  int files_cnt;
  struct file_buffer *next;     /* pointer to the next record */
};

int save_to_buffer_node(struct file_buffer *fbuf, \
                int offset, int size, char * data);

int read_from_buffer_node(struct file_buffer *fbuf, \
                int offset, int size, char * data);

int flush_data_to_disk(char *file_name, struct file_buffer *fbuf);

int copy_buffer_file(struct file_buffer *src_fbuf, struct file_buffer *dst_fbuf);


/*----------------- linked list (buffer node)--------------------------*/
static struct buffer_node *buffer_node_curr __attribute__((unused));

struct buffer_node* create_buffer_node_list(struct file_buffer *filebuf, int offset, char *ptr);

struct buffer_node* add_to_buffer_node_list(struct file_buffer *filebuf, int offset, char *ptr, int add_to_end);

int delete_from_buffer_node_list(struct file_buffer *filebuf, int pos);

int delete_buffer_node_list(struct file_buffer *filebuf);

int get_buffer_list_size(struct file_buffer *filebuf);

struct buffer_node *search_buffer_node(struct file_buffer *filebuf, int pos, struct buffer_node **prev);

/*----------------- linked list (file buffer)--------------------------*/
static struct file_buffer *head __attribute__((unused));
static struct file_buffer *curr __attribute__((unused));

struct file_buffer* create_list(char *file_path);

struct file_buffer* add_to_list(char *file_path, int add_to_end);

struct file_buffer* search_in_list(char *file_path, struct file_buffer **prev);

int delete_from_list(char *file_path);

int delete_list();

int list_get_size();

struct file_buffer *search_one_node(int pos);
void print_all_list();

struct file_buffer* match_in_list(char *file_path);
#endif
