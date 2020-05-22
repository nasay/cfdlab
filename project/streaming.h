#ifndef _STREAMING_H_
#define _STREAMING_H_

#include "utils.h"
/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming (Fields &fields,
                    int * length,           /*Interior length of the cavity */
                    int n_threads,          /*Number of threads to be used */
                    float exchange);        /*Exchange factor, this can make more splashing effects but can introduce instabilities*/

#endif
