/*  transaction.h -- ref-cache http transactions

    Copyright (C) 2025 Genome Research Ltd.

    Author: Rob Davies <rmd@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef TRANSACTION_H_INCLUDED
#define TRANSACTION_H_INCLUDED

#include <stdint.h>
#include <stddef.h>

#include "../htslib/hts_defs.h"
#include "types.h"

typedef enum {
    WRITE_BLOCKED,
    WRITE_BLOCKED_UPSTREAM,
    WRITE_MORE,
    WRITE_COMPLETE,
    WRITE_EOF,
    WRITE_ERROR
} Write_result;

Transaction *new_transaction(Client *client, Http_Parser *parser);
void free_transaction(Transaction *transact);
void free_transaction_list(Transaction *head);
Transaction * switch_to_next_transaction(Transaction *transact);
Client * transaction_get_client(Transaction *transact);
int transaction_get_keep_alive(Transaction *transact);
void transaction_set_ref(Transaction *transact, RefFile *ref);
void transaction_set_req_str(Transaction *transact, const char *requested)
HTS_ACCESS(read_only, 2);
void set_error_response(Transaction *transact, unsigned int code);

void set_message_response(Transaction *transact, const char *content_type,
                          const char *message, size_t len)
HTS_ACCESS(read_only, 2)
HTS_ACCESS(read_only, 3, 4);

Write_result transaction_send_data(Transaction *transact, int fd);
int transaction_have_content(Transaction *transact);
void transaction_set_next(Transaction *transact, Transaction *next);
int transaction_has_data_to_send(Transaction *transact);
void set_transaction_file_range(Transaction *transact, int64_t size,
                                int ref_data_available);

void got_download_started(unsigned int id, int64_t val, int fd,
                          Client **write_stack);
void got_download_part(unsigned int id, int64_t val, Client **write_stack);
void got_download_clen(unsigned int id, int64_t val, Client **write_stack);
void got_download_result(unsigned int id, int64_t val, Client **write_stack);


size_t make_log_message(Transaction *transact, char *buffer, size_t size);

#endif /* TRANSACTION_H_INCLUDED */
