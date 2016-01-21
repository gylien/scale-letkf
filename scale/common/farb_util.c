/*
 * Copyright (C) 2015, Advanced Institute for Computational Science, RIKEN
 * Author: Jianwei Liao(liaotoad@gmail.com)
 */

#include "farb_util.h"

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

struct buffer_node* create_buffer_node_list(struct file_buffer *filebuf, int offset, char *bptr)
{
	struct buffer_node *ptr = (struct buffer_node*)malloc(sizeof(struct buffer_node));
    if(NULL == ptr)
	{
		printf("\n Node creation failed \n");
		return NULL;
	}
	ptr->offset_start = offset;
	ptr->node_ptr = bptr;
	ptr->dirty_flag = 0;
	ptr->next = NULL;
	filebuf->buffer_list = ptr; 
	buffer_node_curr = ptr;
	return ptr;
}

struct buffer_node* add_to_buffer_node_list(\
			struct file_buffer *filebuf, int offset, char *bptr, int add_to_end)
{
    fflush(stdout);
    if(NULL == filebuf->buffer_list)

	{
		return (create_buffer_node_list(filebuf, offset, bptr));
	}
	struct buffer_node *ptr = (struct buffer_node*)malloc(sizeof(struct buffer_node));
	if(NULL == ptr)
	{
		printf("\n Node creation failed \n");
		return NULL;
	}
	ptr->offset_start = offset;
	ptr->node_ptr = bptr;
	ptr->dirty_flag = 0;
	ptr->next = NULL;

	buffer_node_curr = filebuf->buffer_list;
	while (buffer_node_curr->next != NULL)
		buffer_node_curr = buffer_node_curr->next;
	
	if(add_to_end)
	{
		buffer_node_curr->next = ptr;
		//buffer_node_curr = ptr;
	}
    else
    {
        ptr->next = filebuf->buffer_list;
        filebuf->buffer_list = ptr;
    }
	return ptr;
}

int delete_from_buffer_node_list(struct file_buffer *filebuf, int pos)
{
	struct buffer_node *prev = NULL;
	struct buffer_node *ptr= search_buffer_node(filebuf, pos, &prev);
    if (ptr == NULL)
    {
		printf("deleted node not found!\n");
		return -1;
    }
    else
    {
        if (prev != NULL)
		{
            prev->next = ptr->next;
		}
        else if (ptr == filebuf->buffer_list)
        {
            filebuf->buffer_list = ptr->next;
        }

        if (ptr == buffer_node_curr)
        {
            buffer_node_curr = prev;
        }
    }
    if (ptr->node_ptr != NULL)
    	free(ptr->node_ptr);

    free(ptr);
	return 0;
}

int delete_buffer_node_list(struct file_buffer *filebuf){
    struct buffer_node *ptr = (struct buffer_node *)filebuf->buffer_list;
	while(ptr != NULL)
	{
		buffer_node_curr = ptr->next;
		if (ptr->node_ptr != NULL)
			free(ptr->node_ptr);

		free(ptr);
		ptr = buffer_node_curr;			   
	}
	filebuf->buffer_list = NULL;
	return 0;
}

int get_buffer_list_size(struct file_buffer *filebuf){
    int i = 0;
	buffer_node_curr = (struct buffer_node *)filebuf->buffer_list;
	while (buffer_node_curr != NULL){
		i++;
		buffer_node_curr = buffer_node_curr -> next;
	}
	return i;
}

struct buffer_node *search_buffer_node(struct file_buffer *filebuf, int pos, struct buffer_node **prev)
{
    int i = 0;
    struct buffer_node *node = (struct buffer_node *)filebuf->buffer_list;
    while (i < pos && node != NULL){
        if (prev != NULL)
			*prev = node;
        i++;
        node = node -> next;
    }
    if (i == pos){
        return node;
	}
    else{
        //*prev = NULL;
        return NULL;
    }
}

/* for file buffer */

struct file_buffer* create_list(char *file_path)
{
    struct file_buffer *ptr = (struct file_buffer*)malloc(sizeof(struct file_buffer));
    if(NULL == ptr)
    {
        printf("\n Node creation failed \n");
        return NULL;
    }
    strncpy(ptr->file_path, file_path, 1024);
    ptr->buffer_pointer = NULL;
	ptr->buffer_list = NULL;
    ptr->buffer_size = 0;
	ptr->files_cnt = 0;
    ptr->next = NULL;

    head = curr = ptr;
    return ptr;
}

struct file_buffer* add_to_list(char *file_path, int add_to_end)
{
    if(NULL == head)
    {
        return (create_list(file_path));
    }

    struct file_buffer *ptr = (struct file_buffer*)malloc(sizeof(struct file_buffer));
    if(NULL == ptr)
    {
        printf("\n Node creation failed \n");
        return NULL;
    }
    strncpy(ptr->file_path, file_path, 1024);
    ptr->buffer_pointer = NULL;
	ptr->buffer_list = NULL;
    ptr->buffer_size = 0;
	ptr->files_cnt = 0;
    ptr->next = NULL;

    if(add_to_end)
    {
        curr->next = ptr;
        curr = ptr;
    }
    else
    {
        ptr->next = head;
        head = ptr;
    }
    return ptr;
}

struct file_buffer* search_in_list(char * file_path, struct file_buffer **prev)
{
    struct file_buffer *ptr = head;
    struct file_buffer *tmp = NULL;
    int found = 0;

    while(ptr != NULL)
    {
        if(strcmp(file_path, ptr->file_path) == 0)
        {
            found = 1;
            break;
        }
        else
        {
            tmp = ptr;
            ptr = ptr->next;
        }
    }
    
    if(found)
    {
        if(prev)
            *prev = tmp;
        return ptr;
    }
    else
    {
        return NULL;
    }
}

int delete_from_list(char *file_path)
{
    struct file_buffer *prev = NULL;
    struct file_buffer *del = NULL;

    del = search_in_list(file_path,&prev);
    if(del == NULL)
    {
        return -1;
    }
    else
    {
        if(prev != NULL)
            prev->next = del->next;
        else if(del == head)
        {
            head = del->next;
        }
        if(del == curr)
        {
            curr = prev;
        }
    }
    if (del->buffer_pointer != NULL)
		free(del->buffer_pointer);

    free(del);
    del = NULL;

    return 0;
}

void print_all_list(){
        curr = head;
        while (curr != NULL){
                printf("file_name [%s]\n", curr->file_path);
                printf("service_name [%s]\n", curr->service_name);
                printf("alias_name [%s]\n", curr->alias_name);
                printf("direction [%s]\n", curr->direction);
                curr = curr -> next;
        }
}

int list_get_size(){
   	int i = 0;
        curr = head;
        while (curr != NULL){
        	i++;
           	curr = curr -> next;
        }
        return i;
}

struct file_buffer * search_one_node(int pos){
    int i = 0;
    curr = head;
    while (i < pos && curr != NULL){
        i++;
        curr = curr -> next;
    }
    if (i == pos)
	return curr;
    else return NULL;
}

int delete_list(void)
{
    struct file_buffer *ptr = head;
    while(ptr != NULL)
    {
        curr = ptr->next;	
		delete_buffer_node_list(ptr);
		
		if (ptr->buffer_pointer != NULL)
		  free(ptr->buffer_pointer);
		free(ptr);
		ptr = curr;
    }
    return 0;
}

struct file_buffer* match_in_list(char * file_path)
{
    struct file_buffer *ptr = head;
    struct file_buffer *tmp = NULL;
    int found = 0;

    while(ptr != NULL)
    {
        if(strstr(file_path, ptr->file_path) != NULL )
        {
            found = 1;
            break;
        }
        else if (strstr(file_path, ptr->alias_name) != NULL && !strcmp(ptr->direction, "rdwr"))
        {
	    found =  1;
	    break;
	}
	else
	{
            tmp = ptr;
            ptr = ptr->next;
        }
    }
    if(found)
    {
        return ptr;
    }
    else
    {
        return NULL;
    }
}
